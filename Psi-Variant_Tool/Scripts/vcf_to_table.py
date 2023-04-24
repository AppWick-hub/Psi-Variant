import vcf
import argparse

parser = argparse.ArgumentParser(description="This program gets the following parameters: 1. vcf with VEP annotation, 2. list of genes (symbol gene at each line), 3. list of Consequence not to includ, 4. MAX AF cutoff to include and 5. Pedigree file. OUTPUT: Filter variants at table format, with parents genotyping ")
parser.add_argument('--vcf', '-v' ,help = "VCF file", type = str)
parser.add_argument('--gene_list',  '-g' ,help = 'file with gene list (symbol), one gene at line - *filter in*', type = str)
parser.add_argument('--ConsequenceExcluded',  '-c' ,help = 'file with consequence terms - *filter out*', type = str, default = "")
parser.add_argument('--output',  '-o' ,help = 'output file', type = str)
parser.add_argument('--AF_cutoff', '-f',help = 'max Allel frequncy cutoff to *filter in*', type = float, default = 0.01)
parser.add_argument('-ped', '-p' ,help = 'ped file at plink format', type = str)

args = parser.parse_args()


my_vcf = vcf.Reader(open(args.vcf),"r")
Autism_genes = open(args.gene_list,"r")
output = open(args.output,"w")
my_ped = open(args.ped, "r")

if (args.ConsequenceExcluded != ""):
    my_ConsequenceNotInclud = open(args.ConsequenceExcluded,"r")
    ConsequenceNotInclud_set = set(my_ConsequenceNotInclud.read().splitlines())
else:
    ConsequenceNotInclud_set = set()

my_ped_lines = my_ped.read().splitlines()

my_ped_dict = dict()
for line in my_ped_lines:
    sp_line = line.split("\t")
    my_ped_dict[sp_line[1]] = {'father' : sp_line[2], 'mother' : sp_line[3], 'line' : line}

my_samples = my_vcf.samples
VEP = my_vcf.infos["CSQ"].desc.split(":")[1].split("|")
Autism_genes_set = set(Autism_genes.read().splitlines())


output.write("SampleID\tCHROM\tPOS\tREF\tALT\tQUAL\tMQ\tSOR\tQD\tFS\tMQRankSum\tReadPosRankSum\tFILTER\tgenotype\tAD_ref\tAD_alt\tDP\tGQ\tSymbol\tcDNA_position\tProtein_position\tExisting_variation\tConsequence\tSIFT\tPolyPhen\tREVEL_score\tMPC_score\tM_CAP_score\tCADD_phred\tMetaLR_score\tMetaSVM_score\tLoFtool\tGENE_PHENO\tCLIN_SIG\tMAX_AF\tMAX_AF_POPS\tAmino_acids\tCodons\tfather_GT\tfather_AD_REF\tfather_AD_ALT\tfather_GQ\tmother_GT\tmother_AD_REF\tmother_AD_ALT\tmother_GQ\n")

for record in my_vcf: 
    if 'QD' not in record.INFO.keys():
        record.INFO['QD']='NA'
    if 'MQRankSum' not in record.INFO.keys():
        record.INFO['MQRankSum']='NA'
    if 'ReadPosRankSum' not in record.INFO.keys():
        record.INFO['ReadPosRankSum']='NA'
    if 'MQ' not in record.INFO.keys():
        record.INFO['MQ']='NA'
    if 'SOR' not in record.INFO.keys():
        record.INFO['SOR']='NA'
    if 'FS' not in record.INFO.keys():
        record.INFO['FS']='NA'
    for temp_my_vep in record.INFO['CSQ']:
        my_vep = temp_my_vep.split("|")
        my_symbol = my_vep[VEP.index("SYMBOL")]
        if my_symbol in Autism_genes_set:
            MAX_AF = my_vep[VEP.index("MAX_AF")]
            Consequence = my_vep[VEP.index("Consequence")]
            if (MAX_AF == "" or float(MAX_AF) < args.AF_cutoff) and (Consequence not in ConsequenceNotInclud_set):
                MAX_AF_POPS = my_vep[VEP.index("MAX_AF_POPS")]
                CLIN_SIG = my_vep[VEP.index("CLIN_SIG")]
                cDNA_position = my_vep[VEP.index("cDNA_position")]
                Protein_position = my_vep[VEP.index("Protein_position")]
                Amino_acids = my_vep[VEP.index("Amino_acids")]
                Codons = my_vep[VEP.index("Codons")]
                GENE_PHENO = my_vep[VEP.index("GENE_PHENO")]
                SIFT =my_vep[VEP.index("SIFT")]
                PolyPhen =my_vep[VEP.index("PolyPhen")]
                REVEL_score =my_vep[VEP.index("REVEL_score")]
                MPC_score =my_vep[VEP.index("MPC_score")]
                M_CAP_score =my_vep[VEP.index("M-CAP_score")]
                CADD_phred =my_vep[VEP.index("CADD_phred")]
                MetaLR_score =my_vep[VEP.index("MetaLR_score")]
                MetaSVM_score =my_vep[VEP.index("MetaSVM_score")]
                LoFtool =my_vep[VEP.index("LoFtool")]
                Existing_variation = my_vep[VEP.index("Existing_variation")]
                my_genotypes = record.get_hom_alts() + record.get_hets()
                for genotype in my_genotypes:
                    if genotype['GQ'] > 20:
                        try:
                            my_father = my_ped_dict[genotype.sample]['father']
                            my_mother = my_ped_dict[genotype.sample]['mother']
                        except:
                            print ("something wrong with ped dict")
                            my_father = 0
                            my_mother = 0
                        if (my_father in my_samples):
                            father_GT = record.genotype(my_father)['GT']
                            father_AD = record.genotype(my_father)['AD']
                            if father_AD==None:
                                    father_AD=['NA','NA']
                            father_GQ = record.genotype(my_father)['GQ']
                            if father_GQ==None:
                                    father_GQ='NA'
                            if (my_mother in my_samples):
                                mother_GT = record.genotype(my_mother)['GT']
                                mother_AD = record.genotype(my_mother)['AD']
                                if mother_AD==None:
                                    mother_AD=['NA','NA']
                                mother_GQ = record.genotype(my_mother)['GQ']
                                if mother_GQ==None:
                                    mother_GQ='NA' 
                                output.write(f"{genotype.sample}\t{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT}\t{record.QUAL}\t{record.INFO['MQ']}\t{record.INFO['SOR']}\t{record.INFO['QD']}\t{record.INFO['FS']}\t{record.INFO['MQRankSum']}\t{record.INFO['ReadPosRankSum']}\t{record.FILTER}\t{genotype['GT']}\t{genotype['AD'][0]}\t{genotype['AD'][1]}\t{genotype['DP']}\t{genotype['GQ']}\t{my_symbol}\t{cDNA_position}\t{Protein_position}\t{Existing_variation}\t{Consequence}\t{SIFT}\t{PolyPhen}\t{REVEL_score}\t{MPC_score}\t{M_CAP_score}\t{CADD_phred}\t{MetaLR_score}\t{MetaSVM_score}\t{LoFtool}\t{GENE_PHENO}\t{CLIN_SIG}\t{MAX_AF}\t{MAX_AF_POPS}\t{Amino_acids}\t{Codons}\t{father_GT}\t{father_AD[0]}\t{father_AD[1]}\t{father_GQ}\t{mother_GT}\t{mother_AD[0]}\t{mother_AD[1]}\t{mother_GQ}\n")   
                        else:
                            output.write(f"{genotype.sample}\t{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT}\t{record.QUAL}\t{record.INFO['MQ']}\t{record.INFO['SOR']}\t{record.INFO['QD']}\t{record.INFO['FS']}\t{record.INFO['MQRankSum']}\t{record.INFO['ReadPosRankSum']}\t{record.FILTER}\t{genotype['GT']}\t{genotype['AD'][0]}\t{genotype['AD'][1]}\t{genotype['DP']}\t{genotype['GQ']}\t{my_symbol}\t{cDNA_position}\t{Protein_position}\t{Existing_variation}\t{Consequence}\t{SIFT}\t{PolyPhen}\t{REVEL_score}\t{MPC_score}\t{M_CAP_score}\t{CADD_phred}\t{MetaLR_score}\t{MetaSVM_score}\t{LoFtool}\t{GENE_PHENO}\t{CLIN_SIG}\t{MAX_AF}\t{MAX_AF_POPS}\t{Amino_acids}\t{Codons}\t\t\t\t\t\t\t\t\n")                