#!usr/bin/perl
use 5.010;

# Files
$input = 'Output/aai/aai_summary.tsv';
$strains = 'intermediaryFiles/genomeList.txt';

# Get the strains
open($inStrains, '<', $strains);
while(<$inStrains>) {
	chomp;
	@line = split("\t", $_);
	push(@strains, @line[0]);
}
close($inStrains);

# Prepare the first row
print('Strains');
foreach $i (@strains) {
	$i =~ s/Genomes\///;
	$i =~ s/.fna//;
	print("\t$i");
}

# Make the matrix
foreach $i (@strains) {
	print("\n$i");
	foreach $j (@strains) {
		if($i eq $j) {
			print("\t100");
		}
		else {
			open($in, '<', $input);
			while(<$in>) {
				if(/$i/) {
					if(/$j/) {
						@line = split("\t", $_);
						print("\t@line[5]");
					}
				}
			}
			close($in);
		}
	}
}
