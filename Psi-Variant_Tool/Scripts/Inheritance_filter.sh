#!/bin/bash
#$ -S /bin/bash
#$ -cwd

source activate Python3_Apurba # Activate the conda environment (e.g., Python3_Apurba) after installing Python version 3 and above by using conda
##################################################################################################################
# Script: Inheritance.py
#            -i: input file (the output (-o) of the first script)   !-- update in each run --!
#            -o: Results file   !-- update in each run --!
###################################################################################################################
    
python Inheritance.py \ # Add the location of the file
    -i vcf_table_rare_LoF_Missense_Conseq_Revised.table \ # Add the location of the file
    -o vcf_table_rare_LoF_Missense_Inheritance.table # Add the location of the file
	
source deactivate