cp ../assembly_annotation/prokka/Sinorhizobium_sp_M14.gbf organism_directory/ # Get the genbank file
python PhiSpy/genbank_to_seed.py organism_directory/Sinorhizobium_sp_M14.gbf organism_directory/ # Convert genbank file into SEED format
PhiSpy/./PhiSpy.py -i organism_directory -o output_directory -c # Run PhiSpy, and chose number 14 for training set
PhiSpy/./PhiSpy.py -i organism_directory -o output_directory_2 -c # Run PhiSpy, and chose number 0 for training set
PhiSpy/./PhiSpy.py -i organism_directory -o output_directory_3 -c # Run PhiSpy, and chose number 4 for training set

