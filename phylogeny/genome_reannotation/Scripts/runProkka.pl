#!usr/bin/perl
use 5.010;

while(<>) {
	@line = split("\t", $_);
	system("prokka --outdir Output/@line[0] --cpus 20 --prefix @line[0] Genomes/@line[0].fna");
}
