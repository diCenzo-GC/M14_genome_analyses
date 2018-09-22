## Prepare necessary directories
mkdir Genomes/
mkdir Output/
mkdir intermediaryFiles/

## Get the genomes list
grep 'Complete' ../one_way_fastani/intermediaryFiles/refseqGenomeInformation.txt > intermediaryFiles/refseqGenomeInformation.txt # Find complete genomes
grep 'Chromosome' ../one_way_fastani/intermediaryFiles/refseqGenomeInformation.txt >> intermediaryFiles/refseqGenomeInformation.txt # Find chromosome level assemblies
grep -v 'RAC02' intermediaryFiles/refseqGenomeInformation.txt > temp.txt # Remove the one odd strain
mv temp.txt intermediaryFiles/refseqGenomeInformation.txt # Move the file
perl Scripts/parseGenomeList.pl intermediaryFiles/refseqGenomeInformation.txt # Parse the NCBI Genome information
mv genomeList.txt intermediaryFiles/ # Move the file
perl Scripts/parseFastani.pl ../one_way_fastani/Output/output_one_way_fastani.txt >> intermediaryFiles/genomeList.txt # Get the unfinished genomes
sort -u intermediaryFiles/genomeList.txt > temp.txt # Remove duplicates
mv temp.txt intermediaryFiles/genomeList.txt # Move the file
perl Scripts/getGenomes.pl intermediaryFiles/genomeList.txt # Get the genomes

## Reannotate the genomes
perl Scripts/runProkka.pl intermediaryFiles/genomeList.txt # Run prokka to annotate the genomes


