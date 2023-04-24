#-------------------------------------Applying Psi-Variant-----------------------------------------
#-------------Author: Apurba Shil (Ph.D. Student), Ben Gurion University----------------------------

# Import the .table file generated from the vcf file
data1=read.delim("vcf_table_rare_LoF_Missense_Conseq.table") #Add the location of the file to import

#Create few fields to use further
library(tidyverse)
data1=data1 %>%
  mutate(Unique_ID=as.character(paste0(data1$SampleID, "-", data1$CHROM, "-", data1$POS))) %>%
  mutate(SIFT_New=as.numeric(gsub(pattern="[^0-9.]", "", data1$SIFT))) %>%
  mutate(PolyPhen_New=as.numeric(gsub(pattern="[^0-9.]", "", data1$PolyPhen))) %>%
  mutate(MPC_New=as.numeric(gsub("\\.*&\\.*", "", data1$MPC_score))) #%>%
  #select(Unique_ID, SIFT_New, SIFT, PolyPhen_New, PolyPhen, MPC_New, MPC_score, everything())

#Use the recommended cut-off of different in-silico tools to be used (recode them as 1 or 0)
data2=data1 %>%
  filter(grepl('missense', Consequence) | grepl('start_lost', Consequence)) %>%
  mutate(M_CAP_N=replace(M_CAP_score, MAX_AF>0.01 & is.na(M_CAP_score), 0 )) %>%
  mutate(M_CAP_N=replace(M_CAP_N, MAX_AF<=0.01 & is.na(M_CAP_N), 
                         median(M_CAP_N, na.rm = TRUE))) %>%
  mutate(M_CAP=replace(M_CAP_N, is.na(MAX_AF) & is.na(M_CAP_N), 
                       median(M_CAP_N, na.rm = TRUE))) %>%
  mutate(MPC=replace(MPC_New, is.na(MPC_New), 
                     median(MPC_New, na.rm = TRUE))) %>%
  mutate(MetaLR=replace(MetaLR_score, is.na(MetaLR_score), 
                        median(MetaLR_score, na.rm = TRUE))) %>%
  mutate(MetaSVM=replace(MetaSVM_score, is.na(MetaSVM_score), 
                         median(MetaSVM_score, na.rm = TRUE))) %>%
  mutate(REVEL=replace(REVEL_score, is.na(REVEL_score), 
                       median(REVEL_score, na.rm = TRUE))) %>%
  mutate(CADD_phred=replace(CADD_phred, is.na(CADD_phred), 
                            median(CADD_phred, na.rm = TRUE))) %>%
  mutate(PolyPhen=replace(PolyPhen_New, is.na(PolyPhen_New), 
                          median(PolyPhen_New, na.rm = TRUE))) %>%
  mutate(SIFT=replace(SIFT_New, is.na(SIFT_New), 
                      median(SIFT_New, na.rm = TRUE))) %>%
  mutate(SIFT_cat=ifelse(SIFT < 0.05, 1, 0)) %>%
  mutate(PolyPhen_cat=ifelse(PolyPhen > 0.446, 1, 0)) %>%
  mutate(MetaLR_cat=ifelse(MetaLR > 0.5, 1, 0)) %>%
  mutate(MetaSVM_cat=ifelse(MetaSVM > 0, 1, 0)) %>%
  mutate(CADD_cat=ifelse(CADD_phred > 20, 1, 0)) %>%
  mutate(REVEL_cat=ifelse(REVEL > 0.5, 1, 0)) %>%
  mutate(MPC_cat=ifelse(MPC >= 2, 1, 0)) %>%
  mutate(M_CAP_cat=ifelse(M_CAP > 0.025, 1, 0)) %>%
  mutate(Total_Score=SIFT_cat+PolyPhen_cat+CADD_cat+REVEL_cat+MPC_cat+M_CAP_cat) %>%
  mutate(x_cat=ifelse(Total_Score==0, "0",
                      ifelse(Total_Score==1, "1",
                             ifelse(Total_Score==2, "2",
                                    ifelse(Total_Score==3, "3",
                                           ifelse(Total_Score==4, "4",
                                                  ifelse(Total_Score==5, "5","6"))))))) %>%
  select(Unique_ID, Total_Score)

#Need to run before applying the inheritance filter
#matching based on the unique ID and then merging data2 with data1 (actual consequence data)
library(tidyverse)
dat=left_join(data1, data2, by="Unique_ID") %>%
  distinct(Unique_ID, .keep_all = TRUE) 

#Exporting the refined results 
write.table(dat, file = "vcf_table_rare_LoF_Missense_Conseq_Revised.table", #Add the location of the file to export
            sep = "\t", row.names = TRUE, col.names = TRUE)

rm(dat, data1, data2)

#Apply the inheritance filter to the exported "vcf_All_Sample_rare_LoF_Missense_Conseq_Revised.table" file.
#and generate vcf_All_Sample_rare_LoF_Missense_inheritance.table

#Then import it
data4=read.delim("vcf_table_rare_LoF_Missense_Inheritance.table") #Add the location of the file to import

library(readxl)
LoFtool_scores=read_excel("LoFtool_Scores.xlsx", sheet="Sheet1") #Add the location of the file to import
quantile(LoFtool_scores$ExACtool_percentile)

LoFtool_scores$LoF_cat=cut(LoFtool_scores$ExACtool_percentile, 
                           breaks=c(0.0000000, 0.2032380, 0.4688254, 0.7344127, 1.0000000),
                           labels=c('< 25%','25-50%', '50-75%','> 75%'))

#Now apply some filters to refine the variant list as folow
library(tidyverse)

data2=left_join(data4, LoFtool_scores, by="Symbol") %>%
  filter(Mut_type_x != 1) %>% 
  filter(Mut_type_x != "compound") %>%
  add_count(CHROM, POS, Mut_type_x, sort=T) %>%
  filter(Mut_type_x != "deNovo" | 
           Mut_type_x == "deNovo" & n < 3) %>% #Remove all the denovo snvs which occured in 3 and more individuals
  select(-n) %>%
  add_count(CHROM, POS, sort=T) %>%
  filter(n < 4) %>% #Any SNVs (including denovo)
  select(-n) %>%
  mutate(father_AD_ALT=as.numeric(as.character(father_AD_ALT))) %>%
  mutate(father_AD_REF=as.numeric(as.character(father_AD_REF))) %>%
  mutate(father_DP=father_AD_ALT+father_AD_REF) %>%
  mutate(mother_AD_ALT=as.numeric(as.character(mother_AD_ALT))) %>%
  mutate(mother_AD_REF=as.numeric(as.character(mother_AD_REF))) %>%
  mutate(mother_DP=mother_AD_ALT+mother_AD_REF) %>%
  filter(father_GQ > 50 & mother_GQ > 50 & 
           father_DP > 20 & mother_DP > 20 & DP > 20 & GQ > 50) %>%
  filter(is.na(Total_Score)=="TRUE" | Total_Score >= 3) %>% # 3 out of 6 in-silioc tools
  mutate(n=nchar(as.character(ALT))) %>%
  select(1:5,n,everything()) %>%
  filter(grepl('missense', Consequence) & n < 4 | !grepl('missense', Consequence)) %>%
  filter(grepl('missense', Consequence) | grepl('start', Consequence) |
           grepl('< 25%', LoF_cat) & !grepl('missense', Consequence) |
           is.na(LoF_cat)=="TRUE" & !grepl('missense', Consequence)) %>% #Applying the LoFTool filter to LoFs only
  select(-c(n,Synonymous,ExACtool_percentile)) 
  
#Concatenation for preparing a Unique id for matching 2 files
data2$Unique_ID=paste0(data2$SampleID, "-", data2$CHROM, "-", data2$POS)

#Then apply the machine learning filter as follow to annotate the false positives
#ML-Filter need to apply
#-------------------------------------------------------------------------
#To apply the ML Based filter
#Score_1: If at a site >=3 SNVs occur then tag them as 1; else 0
#-----------------------------------------------------------------
Final_data=data2
Final_data=Final_data%>% 
  add_count(CHROM, POS, sort=T) %>%
  mutate(score_1=ifelse(n >=3, 1, 0)) %>%
  mutate(score_3=ifelse(grepl('VQSR', FILTER), 1, 0)) %>%
  mutate(score_4=ifelse(grepl('ExcessH', FILTER), 1, 0)) %>%
  mutate(score_5=ifelse(grepl('deleterio', SIFT) & grepl('benig', PolyPhen), 1, 
                        ifelse(grepl('tolerate', SIFT) & grepl('damag', PolyPhen), 1, 0))) 

#Prepare a variable using AD_ALT and AD_Ref; which can be used as a predictor
#to classify the false positives

Final_data$AD_alt=as.numeric(as.character(Final_data$AD_alt))
Final_data$AD_ref=as.numeric(Final_data$AD_ref)

#The following needs to run only if the inheritance pattern will be provided
v=Final_data%>%
  filter(!grepl('Homo', Mut_type_x)) %>%
  mutate(prop=round((AD_alt/AD_ref)*100, 2))

j=Final_data%>%
  mutate(prop=100) %>%
  filter(grepl('Homo', Mut_type_x))

Final_data=rbind(v,j)

#-----------------------------------------------Need to run this-----------

#Replace missing observation with 0 under the newly created proportion variable
Final_data$prop[is.na(Final_data$prop)]=0

#Replace missing observation with 0 under the newly created proportion variable
Final_data$prop[is.infinite(Final_data$prop)]=0

Final_data$Repeat_pos=as.numeric(Final_data$n)#Before running it, please run score_1 above
Final_data$DP=as.numeric(Final_data$DP)
Final_data$QD=as.numeric(Final_data$QD)
Final_data$QUAL=as.numeric(Final_data$QUAL)
Final_data$SOR=as.numeric(Final_data$SOR)
Final_data$MQ=as.numeric(Final_data$MQ)
Final_data$FS=as.numeric(Final_data$FS)
Final_data$MQRankSum=as.numeric(Final_data$MQRankSum)
Final_data$ReadPosRankSum=as.numeric(Final_data$ReadPosRankSum)

#Since ReadPosRankSum has missing observation we need to impute 0's for those observations.
#This will give advantage in the model building process
Final_data$Repeat_pos[is.na(Final_data$Repeat_pos)]=0
Final_data$DP[is.na(Final_data$DP)]=0
Final_data$QD[is.na(Final_data$QD)]=0
Final_data$QUAL[is.na(Final_data$QUAL)]=0
Final_data$SOR[is.na(Final_data$SOR)]=0
Final_data$MQ[is.na(Final_data$MQ)]=0
Final_data$FS[is.na(Final_data$FS)]=0
Final_data$MQRankSum[is.na(Final_data$MQRankSum)]=0
Final_data$ReadPosRankSum[is.na(Final_data$ReadPosRankSum)]=0

Final_data$score_3=as.factor(Final_data$score_3)
levels(Final_data$score_3)=c('Absent', 'present')
Final_data$score_4=as.factor(Final_data$score_4)
levels(Final_data$score_4)=c('Absent', 'present')
Final_data$score_5=as.factor(Final_data$score_5)
levels(Final_data$score_5)=c('Absent', 'present')

#Lets create a new variable SOR_FS by combining 2 highly collinear variables
#------------------------------Need to run-------------
Final_data$SOR_FS=Final_data$SOR+Final_data$FS

# Prediction on unknown validation (uncertain validation of mutation) set
#This step will be performed once we will finalize our best fitted model

#Now we can load our trained model in future as follow (which is stored in the server)
#Give a new name to the loaded model
mod_fit=readRDS("model_GBM.rds") #Add the location of the file to import

#Some Diagnostics about our trained model which we loaded here
mod_fit$results
mod_fit$finalModel

#Now provide predictions
#install.packages("caret")
library(caret) #Need to load the caret package
#install.packages("gbm")
library(gbm)
Final_data$pred_prob=round(predict(mod_fit, type="prob",newdata=Final_data),3)
P <- predict(mod_fit, Final_data)

Final_data$Pred_Label=P

#Applying some filters to remove the potential false positives
xx=Final_data %>%
  select(-n) %>%
  filter(Pred_Label=="Real" & FILTER== "[]") %>%
  add_count(SampleID, CHROM, Symbol,sort=T) %>%
  select(-c(prop, Repeat_pos, SOR_FS, SOR, score_1, score_3, score_4, score_5,
            QUAL, MQ, QD, FS, MQRankSum, ReadPosRankSum))

#Exporting the refined results of gene-panel analysis
write.table(xx, file = "Psi_Variant_output.table", #Add the location of the file to export
            sep = "\t", row.names = TRUE, col.names = TRUE)
			
# ------------------------------------------------------------The End------------------------------------------------