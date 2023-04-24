import pandas as pd
import numpy as np
import argparse

#filter by genotype (homozygots, conpund hetrozygots and denovo mut) 
parser = argparse.ArgumentParser(description="This program filter by genotype (homozygots, conpund hetrozygots and denovo mut) by getting the results of filter_AutismGenes_to_table_new.py (VCF at table format with parents genotype) and return and filter information about the inheritance mutation behavior (denovo mutation, compound hetrozygot, homo, etc.)")
parser.add_argument('--input', '-i' ,help = "VCF file at table format with parents genotype columns", type = str)
parser.add_argument('--output',  '-o' ,help = 'output file', type = str)

args = parser.parse_args()

my_data = pd.read_csv(args.input, sep = "\t")

my_child_data = my_data[my_data['mother_GT'].notnull() | my_data['father_GT'].notnull()]

my_child_data['Mut_type'] = '1'

homo_mut = my_child_data[(my_child_data['genotype'].str.contains('[1-9]/[1-9]', regex = True)) & (my_child_data['father_GT'].str.contains('0/[1-9]', regex = True))& (my_child_data['mother_GT'].str.contains('0/[1-9]', regex = True)) ]
my_child_data.loc[homo_mut.index,'Mut_type'] = "Homo"
denovo_mut = my_child_data[(my_child_data['genotype'].str.contains('0/[1-9]', regex = True)) & (my_child_data['father_GT'].str.contains('0/0'))& (my_child_data['mother_GT'].str.contains('0/0')) ]
my_child_data.loc[denovo_mut.index,'Mut_type'] ="deNovo"
homo_by_denovo_mut = my_child_data[(my_child_data['genotype'].str.contains('[1-9]/[1-9]')) & (((my_child_data['father_GT'].str.contains('0/[1-9]', regex = True))& (my_child_data['mother_GT'].str.contains('0/0'))) | ((my_child_data['father_GT'].str.contains('0/0'))& (my_child_data['mother_GT'].str.contains('0/[1-9]', regex = True)))) ]
my_child_data.loc[homo_by_denovo_mut.index,'Mut_type'] = "Homo_by_deNovo"
homo_by_duble_denovo_mut = my_child_data[(my_child_data['genotype'].str.contains('[1-9]/[1-9]')) & (my_child_data['father_GT'].str.contains('0/0'))& (my_child_data['mother_GT'].str.contains('0/0')) ]
my_child_data.loc[homo_by_duble_denovo_mut.index, 'Mut_type'] = "Homo_by_duble_deNovo"

compound_het = my_child_data[my_child_data['genotype'].str.contains("0/[1-9]" ,regex = True)]
compound_het_step1 = compound_het[compound_het['father_GT'] != compound_het['mother_GT']]
compound_het_step2 = compound_het_step1[~(compound_het_step1['father_GT'].str.contains('[1-9]/[1-9]', regex = True) | compound_het_step1['mother_GT'].str.contains('[1-9]/[1-9]', regex = True) )]
compound_het_step3 = compound_het_step2.groupby(['SampleID','Symbol','mother_GT','father_GT']).genotype.nunique()
compound_het_step4 = compound_het_step3[compound_het_step3.reset_index().groupby(["SampleID","Symbol"]).size() >1]
compound_het_list = compound_het_step4.reset_index()[['SampleID','Symbol']].drop_duplicates()
compound = pd.merge(my_child_data,compound_het_list, how = 'inner', on = ["SampleID","Symbol"])
compound['Mut_type'] = 'compound'

final_table = pd.merge(my_child_data,compound[["SampleID","Symbol","CHROM","POS","ALT", "Mut_type" ]], how = 'left', on = ["SampleID","Symbol","CHROM","POS","ALT" ])
final_table.loc[final_table.Mut_type_y=='compound', "Mut_type_x"] = "compound"
final_table = final_table.drop("Mut_type_y", axis=1)

final_table.to_csv(args.output, sep = "\t", index = False)



