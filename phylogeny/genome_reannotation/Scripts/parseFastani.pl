#!usr/bin/perl
use 5.010;

$genomeList = '../one_way_fastani/intermediaryFiles/genomeList.txt';

while(<>) {
	@line = split(' ', $_);
	if(@line[2] >= 85) {
		@line[1] =~ s/Genomes\///;
		@line[1] =~ s/.fna//;
		open($genomes, '<', $genomeList);
		while(<$genomes>) {
			if(/@line[1]/) {
				print("$_");
			}
		}
		close($genomes);
	}
}
