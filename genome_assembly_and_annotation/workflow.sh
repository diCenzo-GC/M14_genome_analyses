## Prepare necessary directories
mkdir Results/
mkdir Results/M14/
mkdir Results/M14/referenceGenomes/

## Perform workflow for sample AM2
spades.py -o Results/M14/ --pe1-1 fastq_files/M14_1_R1.fastq.gz --pe1-2 fastq_files/M14_1_R2.fastq.gz --pe2-1 fastq_files/M14_2_R1.fastq.gz --pe2-2 fastq_files/M14_2_R2.fastq.gz -t 32 # Assemble the reads
perl removeShort.pl Results/M14/scaffolds.fasta > Results/M14/scaffolds_shortRemoved.fasta # Remove scaffolds below 200 bp
perl removeLowCov.pl Results/M14/scaffolds_shortRemoved.fasta > Results/M14/scaffolds_shortRemoved_lowRemoved.fasta # Remove scaffolds below 50x coverage
fastANI -q Results/M14/scaffolds_shortRemoved_lowRemoved.fasta --rl genomePaths.txt -o Results/M14/scaffolds_fastANI_output.txt # Run fastANI to find most related genomes
perl parseFastANI.pl Results/M14/scaffolds_fastANI_output.txt # Collect 10 most related genomes (first update script to be correct path)
java -jar /home/medusa/medusa/medusa.jar -scriptPath /home/medusa/medusa/medusa_scripts/ -v -o Results/M14/scaffolds_medusa.fasta -i Results/M14/scaffolds_shortRemoved_lowRemoved.fasta -f Results/AM2/referenceGenomes/ # Run MeDuSa to try and combine scaffolds
prokka --outdir Results/M14/prokka/ --cpus 32 --rfam --prefix Sinorhizobium_sp_M14 Results/M14/scaffolds_medusa.fasta # Annotate the assembly
