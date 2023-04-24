# Psi-Variant
This is an integrative in-house tool to detect ultra rare Likely gene disrupting (LGD) SNVs in Whole Exome Sequencing (WES) data from simplex/multiplex ASD affected families. 

After cleaning the vcf file based on DP, GQ, GATK's False positive filters, denovo false positive filter, and an in-house machine learning filtering algorithm, the Psi-Variant workflow utilizes Ensembl’s VEP to annotate the functional consequences for each variant in a multi- sample vcf file. 

Then, all LoF (frameshift indels, stop gain, splice acceptor/donor, etc.) SNVs will be prioritized by LoFtool. 

Further, three out of six different in-silico tool based scores (SIFT, PolyPhen-2, CADD, REVEL, M_CAP and MPC) will be used to prioritizes missense SNVs (as per the recommended cut-off the in-silico scores). However, one can further opt for one out of six or two out of six scores, and so on to prioritize missense SNVs. Also its possible to relax the cut-off of these scores if required. These scores were extracted by utilizing the dbNSFP database using VEP. 

At present, the tool only includes the (a) denovo, (b) autosomal recessive and (c) x-linked type of SNVs, but one can extend this approach for (d) dominant (e) compound heterozygote, etc. type of SNVs as well. A dedicated script (Inheritance.py) was developed in Python for this purpose and is available on the Scripts directory of Psi-Variant.

A population allele frequency cut-off of < 1% has been used, one can flexibly extend (the allele frequency cut off) further by updating the script. 

Psi-Variant also integrates with the ACMG/AMP SNV interpretation tools e.g., InterVar and TAPES. 

Interestingly, Psi-Variant detects ultra-rare LGD SNVs in any genes (or any gene-panel) having any functional consequences (e.g., LOFs/non-LoFs). 

All the scripts/input files/outputs related to the paper can be found on the "Paper_Related" directory. However, all the R, Python, bash scripts and dummy vcf data can be accessed on the "Psi-Variant" directory.

A brief tutorial (to implement Psi-Variant on a dummy WES dataset) will be described very soon on this page.

For any further query, please contact: apurba.shil316@gmail.com (https://twitter.com/ApurbaShil19) or idanmen@bgu.ac.il (https://www.idanme.com/).


