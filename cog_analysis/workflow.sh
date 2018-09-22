emapper.py -i Proteomes/Sinorhizobium_sp_M14.faa --output Output/Sinorhizobium_sp_M14 -d bact --usemem --cpu 1
emapper.py -i Proteomes/Sinorhizobium_sp_A49.faa --output Output/Sinorhizobium_sp_A49 -d bact --usemem --cpu 1
emapper.py -i Proteomes/Ensifer_adhaerens_OV14.faa --output Output/Ensifer_adhaerens_OV14 -d bact --usemem --cpu 1
emapper.py -i Proteomes/Ensifer_adhaerens_Casida_A.faa --output Output/Ensifer_adhaerens_Casida_A -d bact --usemem --cpu 1

mkdir COG_categories/
mv Output/ emapperOutput/
perl COGanalysis.pl
