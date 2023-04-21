# Psi-Variant
This is an integrative in-house tool to detect ultra rare Likely gene disrupting (LGD) SNVs in Whole Exome Sequencing (WES) data from simplex/multiplex ASD affected families. 

After cleaning the vcf file based on DP, GQ, GATK's False positive filters, denovo false positive filter, and an in-house machine learning filtering algorithm, the Psi-Variant workflow utilizes Ensembl’s VEP to annotate the functional consequences for each variant in a multi- sample vcf file. 

Then, all LoF (frameshift indels, nonsense and splice acceptor/donor) SNVs will be prioritized by LoFtool. 

Further, six different in-silico tool based scores (SIFT, PolyPhen-2, CADD, REVEL, M_CAP and MPC) were used to prioritizes non-LoF SNVs. These scores were extracted by utilizing the dbNSFP database. 

At present, the tool only includes the (a) denovo, (b) autosomal recessive and (c) x-linked type of SNVs, but one can extend this approach for (d) dominant (e) compound heterozygote, etc. type of SNVs as well. A dedicated script (Inheritance.py) was developed in Python for this purpose.

A population allele frequency cut-off of < 1% has been used, one can flexibly extend further by updating the script. 

Psi-Variant also integrates with the ACMG/AMP SNV interpretation tool InterVar and TAPES. 

Interestingly, Psi-Variant detects LGD SNVs in any genes (or any gene-panel) having any functional consequences (e.g., LOFs/non-LoFs). 

A brief tutorial (to implement Psi-Variant on a dummy WES dataset) will be described very soon on this page.

A manuscript has been submitted already to a journal and the preprint will be available soon (the preprint link will be added here). 

For any query, please contact mea at apurba.shil316@gmail.com.


