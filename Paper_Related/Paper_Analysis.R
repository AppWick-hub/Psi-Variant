#-----Effectiveness and concordance of InterVar, TAPES and Psi-Variant-Paper analysis-----
#-----------Author: Apurba Shil----------------------

#---------Import InterVar output------------------
InterVar=read.delim("InterVar_Output.table") #Add the location of the file

library(tidyverse)
I=InterVar %>%
  select(Unique_ID, SampleID, X.Chr, Start, Ref, Alt, Ref.Gene, ExonicFunc.refGene, genotype, 
         father_GT, mother_GT, Mut_type_x, ACMG_InterVar, Prediction_ACMG_Revised) %>%
  rename(ACMG_Rev_I=Prediction_ACMG_Revised) %>%
  rename(ACMG_Old_I=ACMG_InterVar) %>%
  mutate(Det_I="InterVar") %>%
  rename(SampleID_I=SampleID) %>%
  rename(CHROM_I=X.Chr) %>%
  rename(POS_I=Start) %>%
  rename(REF_I=Ref) %>%
  rename(ALT_I=Alt) %>%
  rename(Symbol_I=Ref.Gene) %>%
  rename(Consequence_I=ExonicFunc.refGene) %>%
  rename(genotype_I=genotype) %>%
  rename(father_GT_I=father_GT) %>%
  rename(mother_GT_I=mother_GT) %>%
  rename(Mut_type_x_I=Mut_type_x)

library(readxl)
SFARI_scores=read_excel("SFARI_Genes_11012022.xlsx", sheet="Sheet 1") #Add the location of the file

library(tidyverse)
InterVar=left_join(InterVar, SFARI_scores, by="Symbol") %>%
  mutate(SFARI_Score=replace(Score, is.na(Score), 0))

#---------Import TAPES Output-----------------
TAPES=read.delim("TAPES_Output.table") #Add the location of the file

T=TAPES%>%
  select(Unique_ID, SampleID, CHROM, POS, Ref, ALT, Gene.refGene, ExonicFunc.refGene, genotype, 
         father_GT, mother_GT, Mut_type_x, Prediction_ACMG_tapes, Prediction_ACMG_Revised) %>%
  rename(ACMG_Rev_T=Prediction_ACMG_Revised) %>%
  rename(ACMG_Old_T=Prediction_ACMG_tapes) %>%
  mutate(Det_T="TAPES") %>%
  rename(SampleID_T=SampleID) %>%
  rename(CHROM_T=CHROM) %>%
  rename(POS_T=POS) %>%
  rename(REF_T=Ref) %>%
  rename(ALT_T=ALT) %>%
  rename(Symbol_T=Gene.refGene) %>%
  rename(Consequence_T=ExonicFunc.refGene) %>%
  rename(genotype_T=genotype) %>%
  rename(father_GT_T=father_GT) %>%
  rename(mother_GT_T=mother_GT) %>%
  rename(Mut_type_x_T=Mut_type_x)

library(tidyverse)
TAPES=left_join(TAPES, SFARI_scores, by="Symbol") %>%
  mutate(SFARI_Score=replace(Score, is.na(Score), 0))

#---------Import Psi-Variant Output-----------------
Psi_Variant=read.delim("Psi-Variant_Output.table") #Add the location of the file

P=Psi_Variant %>%
  mutate(Det_P="Psi-Variant") %>%
  rename(SampleID_P=SampleID) %>%
  rename(CHROM_P=CHROM) %>%
  rename(POS_P=POS) %>%
  rename(REF_P=REF) %>%
  rename(ALT_P=ALT) %>%
  rename(Symbol_P=Symbol) %>%
  rename(Consequence_P=Consequence) %>%
  rename(genotype_P=genotype) %>%
  rename(father_GT_P=father_GT) %>%
  rename(mother_GT_P=mother_GT) %>%
  rename(Mut_type_x_P=Mut_type_x) %>%
  select(Det_P, Unique_ID, SampleID_P, CHROM_P, POS_P, REF_P, ALT_P, Symbol_P, Consequence_P,
         genotype_P, father_GT_P, mother_GT_P, Mut_type_x_P)

library(tidyverse)
Psi_Variant=left_join(Psi_Variant, SFARI_scores, by="Symbol") %>%
  mutate(SFARI_Score=replace(Score, is.na(Score), 0))

#----------Merge I and T-----------------
I_T=full_join(I, T, by="Unique_ID", keep=FALSE) %>%
  mutate(Det_I=ifelse(is.na(Det_I)==TRUE, 0, Det_I)) %>%
  mutate(Det_T=ifelse(is.na(Det_T)==TRUE, 0, Det_T)) %>%
  mutate(Det=ifelse(Det_I=="InterVar" & Det_T=="TAPES", "InterVar-TAPES", 
                    ifelse(Det_I==0, "TAPES", 
                           ifelse(Det_T==0, "InterVar", NA)))) 

#-------------Merge I-T and P_3-----------
I_T_P=full_join(I_T, P, by="Unique_ID", keep=FALSE) %>%
  mutate(Det_I=ifelse(is.na(Det_I)==TRUE, 0, Det_I)) %>%
  mutate(Det_T=ifelse(is.na(Det_T)==TRUE, 0, Det_T)) %>%
  mutate(Det_P=ifelse(is.na(Det_P)==TRUE, 0, Det_P)) %>%
  mutate(Detected=paste(Det_I, Det_T, Det_P, sep='-')) %>%
  mutate(Detected=gsub("0", "", Detected)) %>%
  mutate(Detected=gsub("-", "", Detected)) %>%
  mutate(Detected=ifelse(Detected=="InterVarPsiVariant", "InterVar-Psi-Variant",
                         ifelse(Detected=="InterVarTAPES", "InterVar-TAPES",
                                ifelse(Detected=="InterVarTAPESPsiVariant", "InterVar-TAPES-Psi-Variant",
                                       ifelse(Detected=="TAPESPsiVariant", "TAPES-Psi-Variant", Detected))))) 

#Exporting the dataset
library(xlsx)
write.xlsx(I_T_P, "I_T_P_Combined_Output.xlsx", #Add the location to save the file
           sheetName = "Combined", col.names = TRUE, row.names = TRUE, append = TRUE)


#Removing all the files generated untill now after exporting the above dataset I_T_P
rm("I",
   "I_T",
   "I_T_P",
   "InterVar",
   "P_3",
   "Psi_Variant_3",
   "T",
   "TAPES")

#Import the Mother file now

library(readxl)
M_3=read_excel("I_T_P_Combined_Output.xlsx", #Add the location of the file
               sheet="Combined")

#------------Annotate SFARI Gene scores to the symbol----------------

library(readxl)
SFARI_scores=read_excel("SFARI_Genes_11012022.xlsx", #Add the location of the file
sheet="Sheet 1") #the database was accessed on 11th January 2022 from the SFARI gene database

library(tidyverse)
Mother_File=left_join(M_3, SFARI_scores, by="Symbol") %>%
  mutate(SFARI_Score=replace(Score, is.na(Score), 0)) %>%
  select(-c(Score, syndromic))

#-------------11 Group wise Analysis--------------------
#-----------------------------1. InterVar-Whole-----------
df1=Mother_File %>%
  filter(grepl("InterVar", Detected)) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl("InterVar", Detected)) %>%
  count(SampleID, sort = TRUE) #To know about the number of probands detected

df1=Mother_File %>%
  filter(grepl("InterVar", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl("InterVar", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(SampleID, sort = TRUE)#To know about the number of probands detected

df1=Mother_File %>%
  filter(grepl("InterVar", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl("InterVar", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(SampleID, sort = TRUE)#To know about the number of probands detected

#------------2. TAPES-Whole-----------------
df1=Mother_File %>%
  filter(grepl("TAPES", Detected)) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl("TAPES", Detected)) %>%
  count(SampleID, sort = TRUE) #To know about the number of probands detected

df1=Mother_File %>%
  filter(grepl("TAPES", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl("TAPES", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(SampleID, sort = TRUE)#To know about the number of probands detected

df1=Mother_File %>%
  filter(grepl("TAPES", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl("TAPES", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(SampleID, sort = TRUE)#To know about the number of probands detected

#-------------3. Psi-Variant-Whole--------------
df1=Mother_File %>%
  filter(grepl("Psi", Detected)) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl("Psi", Detected)) %>%
  count(SampleID, sort = TRUE) #To know about the number of probands detected

df1=Mother_File %>%
  filter(grepl("Psi", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl("Psi", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(SampleID, sort = TRUE)#To know about the number of probands detected

df1=Mother_File %>%
  filter(grepl("Psi", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl("Psi", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(SampleID, sort = TRUE)#To know about the number of probands detected

#---------------4. InterVar+TAPES (Union)---------------
df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('TAPES', Detected)) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('TAPES', Detected)) %>%
  count(SampleID, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('TAPES', Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('TAPES', Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(SampleID, sort = TRUE)

df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('TAPES', Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('TAPES', Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(SampleID, sort = TRUE)

#---------------------5. InterVar+Psi-Variant (Union)---------------------
df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('Psi', Detected)) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('Psi', Detected)) %>%
  count(SampleID, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('Psi', Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('Psi', Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(SampleID, sort = TRUE)

df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('Psi', Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl('InterVar', Detected) |
           grepl('Psi', Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(SampleID, sort = TRUE) 

#----------------------6. TAPES+Psi-Variant (Union)-------------------------
df1=Mother_File %>%
  filter(grepl('TAPES', Detected) |
           grepl('Psi', Detected)) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl('TAPES', Detected) |
           grepl('Psi', Detected)) %>%
  count(SampleID, sort = TRUE)

df1=Mother_File %>%
  filter(grepl('TAPES', Detected) |
           grepl('Psi', Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl('TAPES', Detected) |
           grepl('Psi', Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(SampleID, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl('TAPES', Detected) |
           grepl('Psi', Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl('TAPES', Detected) |
           grepl('Psi', Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(SampleID, sort = TRUE) 

#-----------------------------7. InterVar-TAPES (Intersection/Overlap)-------------------
df1=Mother_File %>%
  filter(grepl("InterVar-TAPES", Detected)) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("InterVar-TAPES", Detected)) %>%
  count(SampleID, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("InterVar-TAPES", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("InterVar-TAPES", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(SampleID, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("InterVar-TAPES", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(grepl("InterVar-TAPES", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(SampleID, sort = TRUE)

#-----------------8. InterVar-Psi-Variant (Intersection/Overlap)------------------------
df1=Mother_File %>%
  filter(grepl("InterVar-Psi-Var", Detected) | 
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("InterVar-Psi-Var", Detected) | 
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  count(SampleID, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("InterVar-Psi-Var", Detected) | 
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("InterVar-Psi-Var", Detected) | 
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(SampleID, sort = TRUE)

df1=Mother_File %>%
  filter(grepl("InterVar-Psi-Var", Detected) | 
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("InterVar-Psi-Var", Detected) | 
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(SampleID, sort = TRUE) 

#--------------9. TAPES-Psi-Variant (Intersection/Overlap)-------------------
df1=Mother_File %>%
  filter(grepl("TAPES-Psi-Var", Detected) |
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("TAPES-Psi-Var", Detected) |
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  count(SampleID, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("TAPES-Psi-Var", Detected) |
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("TAPES-Psi-Var", Detected) |
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  filter(SFARI_Score==1) %>%
  count(SampleID, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("TAPES-Psi-Var", Detected) |
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(grepl("TAPES-Psi-Var", Detected) |
           grepl("InterVar-TAPES-Psi-Var", Detected)) %>%
  filter(SFARI_Score !=0) %>%
  count(SampleID, sort = TRUE) 

#---------------10. All_Gene-InterVar-TAPES (Intersection/Overlap)-------------------
df1=Mother_File %>%
  filter(Detected=="InterVar-TAPES-Psi-Variant") %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  filter(Detected=="InterVar-TAPES-Psi-Variant") %>%
  count(SampleID, sort = TRUE) 

df1=Mother_File %>%
  filter(Detected=="InterVar-TAPES-Psi-Variant") %>%
  filter(SFARI_Score==1) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(Detected=="InterVar-TAPES-Psi-Variant") %>%
  filter(SFARI_Score==1) %>%
  count(SampleID, sort = TRUE)

df1=Mother_File %>%
  filter(Detected=="InterVar-TAPES-Psi-Variant") %>%
  filter(SFARI_Score !=0) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(Detected=="InterVar-TAPES-Psi-Variant") %>%
  filter(SFARI_Score !=0) %>%
  count(SampleID, sort = TRUE)

#-----------------11. Psi-Variant+InterVar+TAPES (Union)-----------------------
df1=Mother_File %>%
  count(Symbol, sort = TRUE) 

df1=Mother_File %>%
  count(SampleID, sort = TRUE)

df1=Mother_File %>%
  filter(SFARI_Score==1) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(SFARI_Score==1) %>%
  count(SampleID, sort = TRUE)

df1=Mother_File %>%
  filter(SFARI_Score !=0) %>%
  count(Symbol, sort = TRUE)

df1=Mother_File %>%
  filter(SFARI_Score !=0) %>%
  count(SampleID, sort = TRUE)

#-----------------OR and chi square calculations for each of the 11 groups------------------

#---------------------InterVar---------------------------
dat2 <- matrix(c(15, 176, 191, 19618),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("SFARI 1","Other Genes")
#rownames(dat2) <- c("InterVar", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

dat2 <- matrix(c(34, 157, 997, 18812),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("All SFARI","Other Genes")
#rownames(dat2) <- c("InterVar", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)
#--------------------TAPES----------------
dat2 <- matrix(c(12, 159, 194, 19635),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("SFARI 1","Other Genes")
#rownames(dat2) <- c("TAPES", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

dat2 <- matrix(c(27, 144, 1004, 18825),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("All SFARI","Other Genes")
#rownames(dat2) <- c("TAPES", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)
#-------------------------Psi-Variant-----------------
dat2 <- matrix(c(18, 394, 188, 19400),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("SFARI 1","Other Genes")
#rownames(dat2) <- c("Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

dat2 <- matrix(c(59, 353, 972, 18616),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("All SFARI","Other Genes")
#rownames(dat2) <- c("Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

#-------------------------InterVar+TAPES------------------
dat2 <- matrix(c(16, 202, 190, 19592),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("SFARI 1","Other Genes")
#rownames(dat2) <- c("InterVar+TAPES", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

dat2 <- matrix(c(36, 182, 995, 18787),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("All SFARI","Other Genes")
#rownames(dat2) <- c("InterVar+TAPES", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

#------------------------------InterVar+Psi-Variant------------------------
dat2 <- matrix(c(24, 472, 182, 19322),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("SFARI 1","Other Genes")
#rownames(dat2) <- c("InterVar+Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

dat2 <- matrix(c(70, 426, 961, 18543),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("All SFARI","Other Genes")
#rownames(dat2) <- c("InterVar+Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

#------------------------------TAPES+Psi-Variant--------------------
dat2 <- matrix(c(22, 443, 184, 19351),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("SFARI 1","Other Genes")
#rownames(dat2) <- c("TAPES+Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

dat2 <- matrix(c(65, 400, 966, 18569),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("All SFARI","Other Genes")
#rownames(dat2) <- c("TAPES+Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

#------------------------------InterVar-TAPES----------------------
dat2 <- matrix(c(11, 132, 195, 19662),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("SFARI 1","Other Genes")
#rownames(dat2) <- c("InterVar-TAPES", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

dat2 <- matrix(c(24, 119, 1007, 18850),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("All SFARI","Other Genes")
#rownames(dat2) <- c("InterVar-TAPES", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

#------------------------------InterVar-Psi-Variant------------------------
dat2 <- matrix(c(9, 90, 197, 19704),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("SFARI 1","Other Genes")
#rownames(dat2) <- c("InterVar-Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

dat2 <- matrix(c(22, 77, 1009, 18892),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("All SFARI","Other Genes")
#rownames(dat2) <- c("InterVar-Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

#------------------------------TAPES-Psi-Variant----------------------
dat2 <- matrix(c(8, 106, 198, 19688),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("SFARI 1","Other Genes")
#rownames(dat2) <- c("TAPES-Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

dat2 <- matrix(c(20, 94, 1011, 18875),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("All SFARI","Other Genes")
#rownames(dat2) <- c("TAPES-Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

#------------------------------InterVar-TAPES-Psi-Variant---------------------------
dat2 <- matrix(c(7, 81, 199, 19713),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("SFARI 1","Other Genes")
#rownames(dat2) <- c("InterVar-TAPES-Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

dat2 <- matrix(c(17, 71, 1014, 18898),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("All SFARI","Other Genes")
#rownames(dat2) <- c("InterVar-TAPES-Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

#------------------------------InterVar+TAPES+Psi-Variant----------------------------
dat2 <- matrix(c(24, 475, 182, 19319),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("SFARI 1","Other Genes")
#rownames(dat2) <- c("InterVar+TAPES+Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

dat2 <- matrix(c(70, 429, 961, 18540),
               ncol=2, byrow=TRUE)
#colnames(dat2) <- c("All SFARI","Other Genes")
#rownames(dat2) <- c("InterVar+TAPES+Psi-Variant", "Raw VCF")
dat2 <- as.table(dat2)
dat2

summary(dat2)

#Generate Odds ratio and relative risk
#install.packages("epitools")
library(epitools)
library(caret)
oddsratio(dat2)
riskratio(dat2)
sensitivity(dat2)
specificity(dat2)
posPredValue(dat2)
negPredValue(dat2)

#----------------Comparing pipelines for all variants-----------------------------
#Generate the venn diagram
#---------------------------
#install.packages("VennDiagram")                       # Install VennDiagram package
library(VennDiagram) 
library(tidyverse)

setwd("~/Outputs/") #Add the location of the graph to save
png("Venn_Diagram-All.png", width = 8, height = 6.8, units = 'in', res = 900)

grid.newpage()                                        # Move to new plotting page
g=draw.triple.venn(area1 = 220, #InterVar
                   area2 = 199, #TAPES
                   area3 = 483, #Psi-Variant
                   n12 = 164, #InterVar-TAPES
                   n23 = 121, #TAPES-Psi-Variant
                   n13 = 102, #Psi-Variant-InterVar
                   n123 = 90, #Psi-Variant-InterVar-TAPES
                   col=c('black', 'black', 'black'),
                   #fill = c(alpha('green',0.8), alpha('red',0.8), alpha('blue',0.8)),
                   lwd=1.8,
                   lty=1,
                   cex=1.50,
                   fontface = "italic",
                   fontfamily = "sans",
                   height = 680 , 
                   width = 480 , 
                   resolution = 1000,
                   compression = "lzw",
                   cat.cex = 1.22,
                   cat.default.pos = "outer",
                   cat.fontfamily = "sans",
                   rotation = 1, print.mode = "raw", sigdigs=3,
                   #print.mode = c("raw","percent"),
                   category = c("InterVar (220)", "TAPES (199)", "Psi-Variant (483)"))

dev.off()

#---------------Stacked CI Plot (All SFARI+SFARI 1)-----------------
library(ggplot2)
library(ggthemes)
#install.packages("ggthemes")

Outcome_order <- c("I",
                   "T",
                   "P",
                   "I∪T",
                   "I∪P",
                   "T∪P",
                   "I∪T∪P",
                   "I∩T",
                   "I∩P",
                   "T∩P",
                   "I∩T∩P")

#this is the first dataset you have
df1 <- data.frame(Outcome=c("I",
                            "T",
                            "P",
                            "I∪T",
                            "I∪P",
                            "T∪P",
                            "I∪T∪P",
                            "I∩T",
                            "I∩P",
                            "T∩P",
                            "I∩T∩P"),
                  OR=c(8.83, 7.73, 4.75, 8.24, 5.43, 5.25, 5.39, 8.51, 10.15, 7.64, 8.73),
                  Lower=c(4.90, 4.00, 2.80, 4.67, 3.42, 3.25, 3.40, 4.26, 4.68, 3.36, 3.60),
                  Upper=c(14.77,13.57,7.57,13.56,8.22,8.08,8.17,15.31,19.40,14.95,17.88))
# add a group column
df1$ASD_Genes <- "SFARI 1"

# create a second dataset, similar format to first
df2 <- data.frame(Outcome=c("I",
                            "T",
                            "P",
                            "I∪T",
                            "I∪P",
                            "T∪P",
                            "I∪T∪P",
                            "I∩T",
                            "I∩P",
                            "T∩P",
                            "I∩T∩P"),
                  OR=c(4.10, 3.53, 3.21, 3.75, 3.18, 3.13, 3.15, 3.79, 5.38, 4.00, 4.49),
                  Lower=c(2.77, 2.28, 2.39, 2.57, 2.43, 2.37, 2.41, 2.38, 3.25, 2.39, 2.55),
                  Upper=c(5.90, 5.27, 4.22, 5.32, 4.10, 4.07, 4.07, 5.80, 8.53, 6.37, 7.48))

# different group
df2$ASD_Genes = "All SFARI"

# combine the two datasets                      
df = rbind(df1, df2)
# you can do the factoring here
df$Outcome = factor (df$Outcome, level=Outcome_order)

#define colours for dots and bars
dotCOLS = c("darkgray", "black")
barCOLS = c("darkgray","black")

setwd("~/Outputs/") #Add the location of the graph to save
png("Stacked CI plot-ORs.png", width = 15, height = 8, units = 'in', res = 1000)

ggplot(df, aes(x=Outcome, y=OR, ymin=Lower, ymax=Upper,
               col=ASD_Genes, fill=ASD_Genes)) +
  #specify position here
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  #geom_hline(yintercept=1, lty=4, lwd=0.4, color="gray") +
  #specify position here too
  geom_point(size=5, shape=21, col="white", stroke = 0.5,
             position=position_dodge(width = 0.5)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  scale_x_discrete(name="Pipeline Combinations") +
  scale_y_continuous(name="Odds Ratio (95% C.I.)", limits = c(0, 20)) +
  #guides(fill=guide_legend(title="ASD Genes"))+
  theme_par(base_size = 18)+
  #coord_flip()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), 
                                    size=20),
        axis.text.y=element_text(size=18, face="bold"))+
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0),
                                    size=20), 
        axis.text.x=element_text(size=18, face="bold",
                                 angle = 15))+
  theme(title=element_text(face=NULL), legend.position = "top") #+
  #labs(caption = "Ref Category: Raw VCF")

dev.off()
#------------------------------

#count the number of SNVs and probands among I, T and P
M_3 %>%
  count(SampleID)

library(tidyverse)
I=M_3 %>%
  filter(grepl("InterVar", Detected)) %>%
  count(SampleID)

T=M_3 %>%
  filter(grepl("TAPES", Detected)) %>%
  count(SampleID)

P=M_3 %>%
  filter(grepl("Psi", Detected)) %>%
  count(SampleID)

#-----Stacked Line diagram----------
library(readxl)
data=read_excel("Graphs_Data.xlsx", #Add the location of the file to import
                sheet = "Line_Graph_Data")

#----Pipeline Combinationwise PPV-------
setwd("~/Outputs/") #Add the location of the graph to save
png("SFARI_1_PPV-pipelinewise.png", width = 10, height = 6, 
    units = 'in', res = 900)

#Drawing first plot-PPV with axis 1
plot(data$PPV_T3, type="b", col="black", pch=20, xaxt="n",
     lty=2, ylab="Positive Predictive Value (PPV)", xlab="Pipeline Combinations", 
     lwd=1.2, cex=2.5, ylim = c(0.03, 0.23), cex.lab=1.3, cex.axis=1.3)

axis(1, at=data$ID, labels=c("I",
                             "T",
                             "P",
                             "I∪T",
                             "I∪P",
                             "T∪P",
                             "I∪T∪P",
                             "I∩T",
                             "I∩P",
                             "T∩P",
                             "I∩T∩P"), cex.axis=1.3)

points(data$PPV_T3_ALL, col="darkgray", pch=20, cex=2.5)
lines(data$PPV_T3_ALL, col="darkgray",lty=2, lwd=1.2, type="b")

legend(x=1, y=0.235, col = c("black", "darkgray"), # Color of lines or symbols
       border = "black", legend=c("SFARI 1", "All SFARI"), # Fill box border color
       lwd=c(1.2, 1.2), lty=c(2, 2), box.lwd=1, box.lty = 1,
       box.col='gray',  bg=NA, pch=c(20, 20))

dev.off()
#-------------------

#two proportions Z test-to know if the proportion is different staistically or not
#-----------------------------Testing whether proportion of frameshift variant is statistically dfferent in one of the pipeline 
#combinations as compared to the vcf 

#1. Raw Output vs InterVar
prop.test(x = c(39, 4232), n = c(220, 1213319))

#prop.test(x = c(39, 4232), n = c(220, 1213319),
#alternative = "greater")

#2. Raw Output vs TAPES
prop.test(x = c(22, 4232), n = c(199, 1213319))

#3. Raw Output vs Psi-Variant
prop.test(x = c(58, 4232), n = c(483, 1213319))

#-----------------------------Testing whether proportion of missense variants is statistically dfferent in one of the pipeline 
#combinations as compared to the vcf 

#1. Raw Output vs InterVar
prop.test(x = c(113, 95919), n = c(220, 1213319))

#2. Raw Output vs TAPES
prop.test(x = c(117, 95919), n = c(199, 1213319))

#3. Raw Output vs Psi-Variant
prop.test(x = c(394, 95919), n = c(483, 1213319))

#-----------------------------Testing whether proportion of stop gain/start loss variants is statistically dfferent in one of the pipeline 
#combinations as compared to the vcf 

#1. Raw Output vs InterVar
prop.test(x = c(16, 2105), n = c(220, 1213319))

#2. Raw Output vs TAPES
prop.test(x = c(13, 2105), n = c(199, 1213319))

#3. Raw Output vs Psi-Variant
prop.test(x = c(16, 2105), n = c(483, 1213319))

#-----------------------------Non-Frameshift----- 
#combinations as compared to the vcf 

#1. Raw Output vs InterVar
prop.test(x = c(42, 4062), n = c(220, 1213319))

#2. Raw Output vs TAPES
prop.test(x = c(43, 4062), n = c(199, 1213319))

#3. Raw Output vs Psi-Variant
prop.test(x = c(2, 4062), n = c(483, 1213319))

#-----------------------------Splice Acceptor/Donor----- 
#combinations as compared to the vcf 

#1. Raw Output vs InterVar
prop.test(x = c(4, 18817), n = c(220, 1213319))

#2. Raw Output vs TAPES
prop.test(x = c(4, 18817), n = c(199, 1213319))

#3. Raw Output vs Psi-Variant
prop.test(x = c(12, 18817), n = c(483, 1213319))

#-----------------------------Denovo----- 
#combinations as compared to the vcf 

#1. Raw Output vs InterVar
prop.test(x = c(202, 43052), n = c(220, 1213319))

#2. Raw Output vs TAPES
prop.test(x = c(188, 43052), n = c(199, 1213319))

#3. Raw Output vs Psi-Variant
prop.test(x = c(169, 43052), n = c(483, 1213319))

#-----------------------------Autosomal Recessive----- 
#combinations as compared to the vcf 

#1. Raw Output vs InterVar
prop.test(x = c(10, 70948), n = c(220, 1213319))

#2. Raw Output vs TAPES
prop.test(x = c(5, 70948), n = c(199, 1213319))

#3. Raw Output vs Psi-Variant
prop.test(x = c(183, 70948), n = c(483, 1213319))

#-----------------------------Sex Linked----- 
#combinations as compared to the vcf 

#1. Raw Output vs InterVar
prop.test(x = c(8, 9103), n = c(220, 1213319))

#2. Raw Output vs TAPES
prop.test(x = c(6, 9103), n = c(199, 1213319))

#3. Raw Output vs Psi-Variant
prop.test(x = c(131, 9103), n = c(483, 1213319))

#------------------Gene-Type wise------
#-----------------------------SFARI 1------------ 
#combinations as compared to the vcf 

#1. Raw Output vs InterVar
prop.test(x = c(15, 19236), n = c(220, 1213319))

#2. Raw Output vs TAPES
prop.test(x = c(12, 19236), n = c(199, 1213319))

#3. Raw Output vs Psi-Variant
prop.test(x = c(21, 19236), n = c(483, 1213319))

#-----------------------------Any SFARI----- 
#combinations as compared to the vcf 

#1. Raw Output vs InterVar
prop.test(x = c(32, 93681), n = c(220, 1213319))

#2. Raw Output vs TAPES
prop.test(x = c(24, 93681), n = c(199, 1213319))

#3. Raw Output vs Psi-Variant
prop.test(x = c(75, 93681), n = c(483, 1213319))

#-----------------------------Other----- 
#combinations as compared to the vcf 

#1. Raw Output vs InterVar
prop.test(x = c(188, 1119638), n = c(220, 1213319))

#2. Raw Output vs TAPES
prop.test(x = c(175, 1119638), n = c(199, 1213319))

#3. Raw Output vs Psi-Variant
prop.test(x = c(408, 1119638), n = c(483, 1213319))

#------------------------------------------------The End-------------------------------------