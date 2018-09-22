#!usr/bin/perl
use 5.010;

# Files
$input = 'Output/fastani_output_twoWay.txt';
$strains = 'intermediaryFiles/genomePaths.txt';

# Get the strains
open($inStrains, '<', $strains);
while(<$inStrains>) {
	chomp;
	push(@strains, $_);
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
						@line = split(' ', $_);
						$ani = (@line[2] + @line[3]) / 2;
						print("\t$ani");
						last;
					}
				}
			}
			close($in);
		}
	}
}
