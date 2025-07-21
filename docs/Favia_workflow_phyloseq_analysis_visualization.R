# This workflow provides a step by step guide to analyze and visualize 
# the phyloseq object of the Favia fragum project.

# The input for this workflow is the hyloseq object available on DRYAD, 
# or the output of the previous workflow



# For information about phyloseq, please see:
# https://bioconductor.org/packages/release/bioc/html/phyloseq.html
BiocManager::install("phyloseq")
library(phyloseq)

# 1. Input phyloseq object
Favia.ps <- readRDS("###EDIT YOUR PATH/Favia.ps.RDS")

# 2. Remove Chloroplasts and Mitochondria
Favia.ps.2 <- subset_taxa(Favia.ps, (Family!="Mitochondria" | is.na(Family)))
Favia.ps.2
Favia.ps.3 <- subset_taxa(Favia.ps.2, (Order!="Chloroplast" | is.na(Order)))
Favia.ps.3

# 3. Remove sample "CC10" because it has really low numbers of reads 
# and it's likely not a good representation of its microbial community.
Favia.ps.4 <- subset_samples(Favia.ps.3, sample_names(Favia.ps.3)!="CC10")
Favia.ps.4

# 4. Use tree glom to agglomerate ASVs. This will reduce our data set.
# For information about speedyseq, please see:
# https://github.com/mikemc/speedyseq
BiocManager::install("speedyseq")
library(speedyseq)

glom.Favia.ps.4<-tree_glom(Favia.ps.4,resolution=0.05)
glom.Favia.ps.4

# 5. Generate beta-diversity ordination using Bray-Curtis distances
ord_Favia_bc <- phyloseq::ordinate(glom.Favia.ps.4, "PCoA", "bray")
library(ggplot2)
glom.ps4.ordination.plot = 
  plot_ordination(glom.Favia.ps.4, ord_Favia_bc, 
                  type = "samples", 
                  color="habitat", 
                  label="sample", 
                  shape = "site") +
  geom_point(size=2) +
  stat_ellipse(aes(group = habitat), linetype = 2, linewidth = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_color_manual(values = c(mangrove = "#a6611a", reef = "#018571"))

glom.ps4.ordination.plot

# 6. Perform PERMANOVA to test for differences in sites and habitats using bray-curtis distances.
# For information on vegan, please see:
# https://github.com/vegandevs/vegan
BiocManager::install("vegan")
library(vegan)

glom.Favia.ps.4.bray <-phyloseq::distance(glom.Favia.ps.4, method="bray")
sampledf <- data.frame(sample_data(glom.Favia.ps.4))
vegan::adonis2(glom.Favia.ps.4.bray~habitat*site, data=sampledf)

# 7. Test for differences in dispersion (B-C Dissimilarity) between groups.
dispr2 <- vegan::betadisper(glom.Favia.ps.4.bray, phyloseq::sample_data(glom.Favia.ps.4)$habitat)
dispr2
p21 = plot(dispr2, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
p22 = boxplot(dispr2, main = "", xlab = "")
vegan::permutest(dispr2)

# 8. Determine "biomarkers", as in Segata (2011)"
# For information on microbiomeMarker, please see:
# https://github.com/yiluheihei/microbiomeMarker

BiocManager::install("microbiomeMarker")
library(microbiomeMarker)

#create a new phyloseq to be able to edit and add taxa names to tax table. 
glom.Favia.ps.5 <- glom.Favia.ps.4
phyloseq::tax_table(glom.Favia.ps.5)[,7] <-phyloseq::taxa_names(glom.Favia.ps.5)

#Linear Discriminant Analysis Effect Size
favia_lefse<-run_lefse(glom.Favia.ps.5, wilcoxon_cutoff = 0.01, group="habitat",
                       kw_cutoff = 0.01, multigrp_strat = TRUE, lda_cutoff = 4)

# the ef_lda measures the relative size of enrichment
# the padj is the p-value adjusted for multiple comparisons
print(marker_table(favia_lefse))
data.frame(marker_table(favia_lefse))

#now we can plot it 
plot_ef_bar(favia_lefse) + 
  scale_fill_manual(values = c(mangrove = "#a6611a", reef = "#018571"))

# 9. Visualize an unrooted tree of Bacteria. 
# The color of the branches indicate the log 2 fold change in median proportion of abundance between habitats. 
# More brown = more enriched in mangroves, more teal = more enriched in reefs.

# For more information on metacoder, please see:
# https://grunwaldlab.github.io/metacoder_documentation/

install.packages("metacoder")
library(metacoder)

#make phyloseq object have the same orientation for all inputs
glom.Favia.ps.4 <- speedyseq::orient_taxa(glom.Favia.ps.4, "rows")

#create a taxmap object
favia.taxmap <- parse_phyloseq(glom.Favia.ps.4)

#caculate abundance of each taxon
favia.taxmap$data$tax_abund <- calc_taxon_abund(favia.taxmap, "otu_table")

#calculate differences in abundance of each taxon between habitats
favia.taxmap$data$diff_table <- compare_groups(favia.taxmap,
                                               dataset = "tax_abund",
                                               cols = sample_data(glom.Favia.ps.4)$sample, 
                                               groups = sample_data(glom.Favia.ps.4)$habitat)


set.seed(070473) #the location of the tree tips are randomized, so setting a seed will allow replicability of the same topology
favia.taxmap %>%
  metacoder::filter_taxa(taxon_names == "Bacteria", subtaxa = TRUE) %>% 
  #metacoder::filter_taxa(taxon_ranks == "Order", supertaxa = TRUE) %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = log2_median_ratio,
            node_color_interval = c(-5,5),
            node_color_range = c("#a6611a","gray90","#018571"),
            node_size_axis = NULL,
            title = "Log 2 ratio of median proportions",
            title_size = 0.05,
            node_color_axis_label = "reef        mangrove",
            node_size_axis_label = "no. taxa within node",
            layout = "da", # The primary layout algorithm
            initial_layout = "fr", 
            overlap_avoidance = 1)

