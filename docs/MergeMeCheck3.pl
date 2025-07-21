#!perl -w

# want to merge two files by merging each line

# get input on command line

# usage: perl MergeMeCheck3.pl inputfilename1 inputfilename2

my $file1 = $ARGV[0];
chomp $file1;
open (IN, "< $file1") || die "\n can't open file: $! \n";

my $file2 = $ARGV[1];
chomp $file2;
open (INTWO, "< $file2") || die "\n can't open file: $! \n";

my $output = "merged.fastq";
open (OUT, "> $output") || die "\n can't open file: $! \n";

print "Welcome to MergeMeCheck3\n";
print "First: $file1\n";
print "Second: $file2\n";

# note: this script is meant for use with fastq files
# do in batches of 4 lines or 2 lines because of pattern

while (my $line = <IN>) {
	my $lineTwo = <INTWO>;
	my $modLine = substr($line,0,-8);
	my $modLineTwo = substr ($lineTwo,0,-8);
	if ($modLine eq $modLineTwo) {
		print OUT $lineTwo; #just keep one of the definition lines, no chomp because want return
		$line = <IN>; #get the sequence itself
		$lineTwo = <INTWO>;
		chomp($line);
		my $newline = $line . $lineTwo;
		print OUT $newline;
	} 
	# by putting these steps inside the if statement, should skip to next line if the definition lines do not match
	# for "no check" removed this part
	# why: to combine an index1 and an index 2, version 2 wants to skip everything
	# can combine I1/I2 into a single 24bp golay barcode
	# to have a check, need to remove final 8 characters of definition line
	# now have this code in version 3
}
close OUT;
print "Done!\n";
