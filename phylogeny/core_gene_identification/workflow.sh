## Prepare necessary directories
mkdir Genomes/
mkdir intermediaryFiles/

## Get the genomes list
cp ../genome_reannotation/intermediaryFiles/genomeList.txt intermediaryFiles # Get the list of genomes
perl Scripts/getGenomes.pl intermediaryFiles/genomeList.txt # Get the gff files

## Run roary to get the core genome
roary -p 20 -f Output -e -i 80 -g 100000 Genomes/*.gff # Run roary

