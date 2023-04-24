#!/bin/bash
#$ -S /bin/bash
#$ -pe shared 40 
#$ -q bioinfo.q
#$ -cwd 

source activate GATK

# ----The following line will remove all those variant sites with GATK's 'VQSR' and 'ExcessHet' tags in the FILTER column in the vcf

bcftools view -f PASS \
dummy_data.vcf.gz \ # Add the location of the file
> dummy_data_filtered.vcf # Add the location of the file

source deactivate

#---4. Ultra rare SNVs annotated with VEP----
source activate ensembl-vep #Install first the Ensembl's VEP under a conda environement (e.g., ensembl-vep in this case)
export PERL5LIB=''
vep \
	--offline  \
	--format vcf \
	--assembly GRCh38 \
	--dir /v100/ \
	--force_overwrite  \
  --everything \
	--vcf  \
  --per_gene \
	--species homo_sapiens \
	--fork 40 \
  --dir_plugins /VEP_plugins-release-100/ \ # Add the location of the directory
  --plugin dbNSFP,/dbNSFP4/dbNSFP4.1a_grch38.gz,MetaSVM_score,MetaLR_score,REVEL_score,CADD_phred,MPC_score,M-CAP_score,FATHMM_score,MutationAssessor_score,VEST4_score,SIFT_score,Polyphen2_HVAR_score,clinvar_clnsig \ # Add the location of the database
  --plugin LoFtool,LoFtool_scores.txt \ # Add the location of the file
	-i dummy_data_filtered.vcf -o stdout | \ # Add the location of the input file
   filter_vep --filter "MAX_AF < 0.01 or not MAX_AF" \
	-o dummy_data_filtered_ultra_rare_VEP.vcf # Add the location of the output file
	
source deactivate

#----VCF to rare LOF+missense SNVs extraction for all samples

source activate Python3_Apurba # Activate the conda environment (e.g., Python3_Apurba) after installing Python version 3 and above by using conda
##################################################################################################################
# Script: vcf to table generation python script
#           -g: file with genes to included (only SNVs in the provided gene list will be included)    !-- chack in each run --!
#           -v: VEP annotated VCF file
#           -p: Pedigree file for the VCF
#           -c: file with Consequence to exluded    !-- check in each run --!
#           -f: Max allel frequncy to include (exlude SNVs with higher population allele frequncy) deafult value 0.01
#           -o: Results file   !-- update in each run --!
###################################################################################################################


python vcf_to_table.py \ # Add the location of the file
    -g Biomart_ensemble_hg37_Gene_Symbols.txt \ # Add the location of the file
    -v dummy_data_filtered_ultra_rare_VEP.vcf \ # Add the location of the file
    -p Family_Pedigree.table \ # Add the location of the file
    -c Func_conseq_to_exclude.txt \ # Add the location of the file
    -o vcf_table_rare_LoF_Missense_Conseq.table # Add the location of the file
	
source deactivate
