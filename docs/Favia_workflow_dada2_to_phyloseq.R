#This workflow provides a step by step guide to process demultiplexed sequences 
#of the Favia fragum project through the DADA2 pipeline 
#for analysis of 16s amplicon libraries.

# The inputs for this workflow are 
# 1. Demultiplexed read 1 and read 2 files from NCBI Sequence Read Archive
# Submission # SUB13878240
# 2. metadata file available in DRYAD and Github
# https://marineinvert.github.io/microbiome/resources.html
# 3. SILVA database available at:
# https://benjjneb.github.io/dada2/training.html


#Install the packages needed

#For additional installation help, please see:
#  https://benjjneb.github.io/dada2/dada-installation.html

#For additional information about using DADA2, please see:
#  https://benjjneb.github.io/dada2/tutorial.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.20")

# For information about DECIPHER, please see: 
# https://www.bioconductor.org/packages/release/bioc/html/DECIPHER.html
BiocManager::install("DECIPHER")

# For information about phangorn, please see:
# https://github.com/KlausVigo/phangorn
install.packages("phangorn")

# For information about phyloseq, please see:
# https://bioconductor.org/packages/release/bioc/html/phyloseq.html
BiocManager::install("phyloseq")

# For information about Biostrings, please see:
# https://bioconductor.org/packages/release/bioc/html/Biostrings.html
BiocManager::install("Biostrings")

path <- "" #include the path of the directory containing the read files
list.files(path). #check the files in the directory

#sort read 1 and read 2 files to make sure they are listed in the same order
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

#extract sample names from the filenames
#the assumption is that filenames have the format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#begin the same steps as the dada2 tutorial
#for details see: https://benjjneb.github.io/dada2/tutorial.html

library(dada2)#load package
plotQualityProfile(fnFs, aggregate = TRUE)
plotQualityProfile(fnRs, aggregate = TRUE)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filter and trim
#trim at length 240 for forward reads and 160 for reverse reads to ensure enough overlap
#between paired reads

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
print(out)

#error rate calculation
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#error correction
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#merge reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#output statistics
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
print(track)

#Assign taxonomy using SILVA database
#the SILVA database is available at: https://benjjneb.github.io/dada2/training.html

taxa <- assignTaxonomy(seqtab.nochim, "###EDIT YOUR PATH/silva_nr99_v138.2_wSpecies_train_set.fasta", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#create distance matrix and gene tree for 'tree_glom' agglomeration
sequences<-getSequences(seqtab.nochim)
names(sequences)<-sequences
alignment <- DECIPHER::AlignSeqs(DNAStringSet(sequences), anchor=NA)
length(alignment)

dm<-DistanceMatrix(alignment,includeTerminalGaps = FALSE)
dm.dataframe<-as.data.frame(dm)
is.na(dm.dataframe)<-sapply(dm.dataframe, is.infinite)
dm.dataframe[is.na(dm.dataframe)]<-1
dm.matrix<-as.matrix(dm.dataframe)

my.upgma<-phangorn::upgma(dm2)

# create phyloseq object
library(phyloseq)
library(Biostrings)

#metadata file available in this dryad submission 
favia <- read.csv("#INSERT PATH TO METADATA FILE HERE/favia_metadata_2.csv") 
samples.out <- rownames(seqtab.nochim)
rownames(favia) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(favia), 
               tax_table(taxa),phy_tree(my.upgma))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)

#change ASV names to numerical ordered ASVs
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


