## Prepare necessary directories
mkdir Genomes/
mkdir Output/
mkdir intermediaryFiles/

## Get the genomes
cp ../genome_reannotation/intermediaryFiles/genomeList.txt intermediaryFiles # Get the list of genomes
perl Scripts/getGenomes.pl intermediaryFiles/genomeList.txt # Get the genomes

## Run fastANI
find Genomes/*.fna > intermediaryFiles/genomePaths.txt # Get the genome paths
fastANI --ql intermediaryFiles/genomePaths.txt --rl intermediaryFiles/genomePaths.txt -o Output/fastani_output.txt # Run fastANI

## Parse fastANI output
sort -k1,1 -k2,2 Output/fastani_output.txt > Output/fastani_output_sorted_1.txt # Sort the file by first column then by second column
sort -k2,2 -k1,1 Output/fastani_output.txt > Output/fastani_output_sorted_2.txt # Sort the file by second column then by first column
cut -d ' ' -f 1,2,3 Output/fastani_output_sorted_1.txt > temp.txt # Get the relevant columns of the first sorted file
cut -d ' ' -f 3 Output/fastani_output_sorted_2.txt > temp2.txt # Get the relevant columns of the second sorted file
paste -d ' ' temp.txt temp2.txt > Output/fastani_output_twoWay.txt # Combine the relevant columns
rm temp* # remove the temporary files
perl Scripts/prepareANImatrix.pl > Output/ANI_matrix.txt # make a two-way ANI matrix from the fastANI output

