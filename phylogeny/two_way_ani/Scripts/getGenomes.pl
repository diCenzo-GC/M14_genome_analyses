#!usr/bin/perl
use 5.010;

while(<>) {
	@line = split("\t", $_);
	$file = '../one_way_fastani/Genomes/' . @line[0] . '.fna';
	system("cp $file Genomes/");
}
