## Prepare necessary directories
mkdir Genomes/
mkdir Output/

## Download genomes
perl Scripts/parseGenomeList.pl intermediaryFiles/refseqGenomeInformation.txt # parse the NCBI Genome information
echo "Sinorhizobium_sp_M14\tXXX" >> genomeList.txt # add M14 to the genome list
sed -i 's/\\t/\t/' genomeList.txt # add a proper tab to the file
mv genomeList.txt intermediaryFiles/ # move the file
perl Scripts/downloadGenomes.pl intermediaryFiles/genomeList.txt # Download the genomes
cp /home/georged/M14_analyses/assembly_annotation/prokka/Sinorhizobium_sp_M14.fna Genomes/ # Get the M14 genome

## Perform the fastANI analysis
find Genomes/*.fna > intermediaryFiles/genomePaths.txt # Get the genome paths
fastANI -q Genomes/Sinorhizobium_sp_M14.fna --rl intermediaryFiles/genomePaths.txt -o Output/output_one_way_fastani.txt # Run fastANI
