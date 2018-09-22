# Remember to first prepare a HMM_file directory containing the seed proteins to build the hidden Markov models, and add hmm_list.txt file to the main directory.

# Set up directories
mkdir Proteomes
mkdir ProteomesHMM
mkdir HMMsearch
mkdir HMMsearchParsed
mkdir HMMsearchHits
mkdir HMMscan
mkdir HMMscanParsed
mkdir HMMscanTop
mkdir OutputFiles
mkdir HMMscanProteins
mkdir intermediaryFiles
mkdir hmmDatabaseFiles
mkdir HMMscanTopLists

# Prepare protein files
cp ../phylogeny/two_way_aai/intermediaryFiles/genomeList.txt intermediaryFiles # Get the list of genomes
mv hmm_list.txt intermediaryFiles # Move this file
perl Scripts/getGenomes.pl intermediaryFiles/genomeList.txt # Get the proteomes
perl Scripts/switchNames.pl # Switch the names of the proteins in the faa files
cat ProteomesHMM/*.faa > intermediaryFiles/combined_proteomes_HMM.faa # Combine the faa files into one file
perl Scripts/modifyFasta.pl intermediaryFiles/combined_proteomes_HMM.faa > intermediaryFiles/combined_proteomes_HMM_modified.faa # Modify the fasta file for easy extraction

# Perform the HMMsearch screens
perl Scripts/performHMMsearch.pl # A short script to repeat for all HMM files, the build, hmmsearch, parsing, and hit extraction

# Perform the HMM scan screens
cp /datadisk1/georged/HMM_Database_Files/* hmmDatabaseFiles/ # get the HMM database
gunzip hmmDatabaseFiles/*.gz # unzip the database files
perl Scripts/performHMMscan.pl # a short script to repeat for all the HMM search output files, to perform hmmscan, parse, and hit extraction
rm hmmDatabaseFiles/* # delete the HMM database files

# Determine strains with each protein
perl Scripts/determineProteinPresence.pl > OutputFiles/proteinPresenceAbsence.txt # determine which of the six proteins are in each of the strains
