#!usr/bin/perl
use 5.010;

while(<>) {
	@line = split("\t", $_);
	system("cp ../phylogeny/genome_reannotation/Output/@line[0]/@line[0].faa Proteomes/");
}
