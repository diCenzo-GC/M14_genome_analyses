## Prepare necessary directories
mkdir Proteomes/
mkdir Output/
mkdir intermediaryFiles/

## Get the proteomes
cp ../genome_reannotation/intermediaryFiles/genomeList.txt intermediaryFiles # Get the list of proteomes
perl Scripts/getGenomes.pl intermediaryFiles/genomeList.txt # Get the proteomes

## Perform AAI analysis
comparem aai_wf -e 1e-12 -p 40.0 -a 70.0 --proteins --sensitive -x faa --tmp_dir /datadisk3/georged/ -c 5 Proteomes/ Output/ # Run compareM to do AAI calculation
perl Scripts/prepareAAImatrix.pl > Output/AAI_matrix.txt # Make a two-way AAI matrix from the compareM output
sed -i 's/_FSM_MA//g' Output/AAI_matrix.txt # Shorten name of one strain
