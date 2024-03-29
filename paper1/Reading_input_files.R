##Start of R session for HCM/DCM analysis
#Wipe all objects from R memory
#rm(list=ls())

library("dplyr")
library("tidyr")
library("ggplot2")
library("readxl")
library("forcats")

## STARTS SECTION THAT NEEDS CHANGING BEFORE EXECUTION
#Data to load True is the 
AllGenes <- T
TTN <- F
TTN_PAN <- F       
sSNV <- F

#Location where the R scripts are
#Rfolder='/Users/lpesce/Documents - MacBook Pro/Beth/Howtos/R/paper1'
#Location of work folder where the dataframes are
#setwd('/Users/lpesce/Documents - MacBook Pro/Beth/Papers/HCM_DCM')
setwd('/home/lpe159/HCM_DCM/sSNV_Bootstrap') #For server

#Depending upon the GTEx database being used the formatting can be different, see below


# data to exclude 
ALMS1 <- T
MultiAllelic <- T
Indels <- F
gnomAD <- F
NU_FQ <- F
Remove_CBL <- F

#Set the interval of ExAC whose variants will be included
FqMin <- 0.00
FqMax <- 1.00


dir() #List files in the workdirectory

#Time stamp
dateStr <- Sys.Date() %>% format(format="%B_%d_%Y")


##Read phenotype Data
#Original
#MD5 (Cohort Patient Data MASTERmrp_machine_cluster.xlsx) = 50ccaa5c9ca1473a1705f3e43f460f1c
#File corrected on 1/19/2018 
#MD5 (Cohort Patient Data MASTERmrp_machine_cluster.xlsx) = 1e4ac601b744bab7be7d18a98f2705ba
#File corrected on 1/22/2018 
#MD5 (Cohort_Patient_Data_MASTERUPDATE01222018.xlsx) = f6ed4d594e5987cc661770fe323a0d26
#File corrected on 1/22/2018 
#MD5 (Cohort_Patient_Data_MASTERUPDATE01232018.xlsx) = 5e538345dd35a79a53a5295c0c4e7b56
#File connected on August 26th 2019
#MD5 (Cohort_Patient_Data_MASTERUPDATE090418.xlsx) = a1fab4ee06840421ae48b1f535cea0bc

subjPheno <-read_excel("Cohort_Patient_Data_MASTERUPDATE090418.xlsx")
subjPheno[! is.na(subjPheno$MegaSeqID),] -> subjPheno #Eliminate columns with gibberish if any
names(subjPheno)[14] <-"Machine_Recode" # Fix name, there are multiple columns with the same name in the original file
subjPheno$Id <- as.factor(subjPheno$ID) #match names to later loaded phenotype file
subjPheno$ID <- NULL
#Reorder fields as necessary for readability
refCols <- c("Id")
subjPheno <- subjPheno[, c(refCols, setdiff(names(subjPheno), refCols))]
str(subjPheno)
#Check some fields for sanity and fix if necessary, NA errors show missing data
is.numeric(subjPheno$EF)
subjPheno$'age 1st present' <- as.numeric(subjPheno$'age 1st present')
#Rename fields to remove names that might make R break
names(subjPheno) <- names(subjPheno) %>% make.names
#Fix Ethnicity gibbersih
subjPheno$Race_Recode_Cluster2 <- as.factor(subjPheno$Race_Recode_Cluster2)
subjPheno %>% group_by(Race.ethnicity,Race_Recode_Cluster2) %>% summarize(n())
# A tibble: 13 x 3
# Groups:   Race.ethnicity [13]
#   Race.ethnicity               Race_Recode_Cluster2 `n()`
#   <chr>                        <fct>                <int>
# 1 AA                           1                        2
# 2 African American             1                       21
# 3 Caucasian                    0                       61
# 4 CAUCASIAN                    0                       25
# 5 Caucasian (French descent)   0                        1
# 6 Caucasian/Hispanic           0                        1
# 7 Hispanic                     0                        9
# 8 HISPANIC                     0                        1
# 9 Hispanic(more than one race) 0                        1
#10 Indian                       0                        1
#11 MIDDLE EASTERN               0                        1
#12 NOT FOUND                    1                        1
#13 Unknown                      0                        2
subjPheno$Race.ethnicity <- subjPheno$Race.ethnicity %>% fct_collapse(AA = c("AA","African American","NOT FOUND"))
subjPheno$Race.ethnicity <- subjPheno$Race.ethnicity %>% fct_collapse(EA = c("Caucasian","CAUCASIAN","Caucasian (French descent)"))
subjPheno$Race.ethnicity <- subjPheno$Race.ethnicity %>% fct_collapse(OA = c("Indian","MIDDLE EASTERN","Unknown"))
subjPheno$Race.ethnicity <- subjPheno$Race.ethnicity %>% fct_collapse(HA = c("Hispanic","HISPANIC","Caucasian/Hispanic","Hispanic(more than one race)"))
subjPheno %>% group_by(Race.ethnicity,Race_Recode_Cluster2) %>% summarize(n())
# A tibble: 4 x 3
# Groups:   Race.ethnicity [4]
#  Race.ethnicity Race_Recode_Cluster2 `n()`
#  <fct>          <fct>                <int>
#1 AA             1                       24
#2 EA             0                       87
#3 HA             0                       12
#4 OA             0                        4


##Read Gene group data
#Date sent by Megan: 10/13/2017
#MD5 (Pan_Card_LMM_GeneDx_Invitae_groups.xlsx) = 504a85ac91e44ca089f6fd3be8c09d5f
geneGroups <- read_excel("Pan_Card_LMM_GeneDx_Invitae_groups.xlsx", sheet=1, col_names=FALSE)
names(geneGroups) <- c("Gene","Group")
geneGroups$Gene <- as.factor(geneGroups$Gene)
geneGroups$Group <- as.factor(geneGroups$Group)



##Read segment duplicaton nr file
#MD5 (Genomes_MegaSeq6.0_KnownGene_totalSegmentDuplicates.dataframe_results) = #819ec847bf0025aaa2930e073adde358
segDupByGene <- read.table('Genomes_MegaSeq6.0_KnownGene_totalSegmentDuplicates.dataframe_results',colClasses = "character", header=TRUE) #dataframe containing all the data
#Fix names (remove ":") and turn number of duplications into numeric type
names(segDupByGene) <- c("Gene","nSegDup")
segDupByGene$Gene <- substr(segDupByGene$Gene,1,nchar(segDupByGene$Gene)-1) %>% as.factor
segDupByGene$nSegDup <- as.numeric(segDupByGene$nSegDup)

##Select a tab separated variants file  -- File has some subjects already removed
# scp lpe159@mcnallylabwkstn01.fsm.northwestern.edu:/data/VariantAnalysis/HCM_DCM_Comparison/HCM_DCM_Analysis/HCM_DCM_August2017/PanCardio/Genomes_MegaSeq6.0_Pan_Card_LMM_GeneDx_Invitae.dataframe .
# 67ab62a36248aa0bc25acb6847f42e77  /data/VariantAnalysis/HCM_DCM_Comparison/HCM_DCM_Analysis/HCM_DCM_August2017/PanCardio/Genomes_MegaSeq6.0_Pan_Card_LMM_GeneDx_Invitae.dataframe
#scp lpe159@mcnallylabwkstn01.fsm.northwestern.edu:/data/VariantAnalysis/HCM_DCM_Comparison/HCM_DCM_Analysis/HCM_DCM_August2017/TTN/Genomes_MegaSeq6.0_TTN.dataframe .
#MD5 (Genomes_MegaSeq6.0_TTN.dataframe) = e337d65e1410745389f08c18176a4ced
#scp lpe159@mcnallylabwkstn01.fsm.northwestern.edu:/data/VariantAnalysis/HCM_DCM_Comparison/HCM_DCM_Analysis/HCM_DCM_August2017/Cardio_TTN/Genomes_MegaSeq6.0_Pan_Card_LMM_GeneDx_Invitae_TTN.dataframe .
#MD5 (Genomes_MegaSeq6.0_Pan_Card_LMM_GeneDx_Invitae_TTN.dataframe) = 23313f9d7d2777fe096a8ce8362cfec9
#scp lpe159@mcnallylabwkstn01.fsm.northwestern.edu:/data/VariantAnalysis/HCM_DCM_Comparison/HCM_DCM_Analysis/HCM_DCM_October2017/PanCardio_TTN/Genomes_MegaSeq6.0_Pan_Card_LMM_GeneDx_Invitae_TTN.dataframe .
#MD5 (Genomes_MegaSeq6.0_Pan_Card_LMM_GeneDx_Invitae_TTN.dataframe) = a51fd2a17e0d43a1c74592c71d7f3457
#scp lpe159@mcnallylabwkstn01.fsm.northwestern.edu: /data/VariantAnalysis/HCM_DCM_Comparison/HCM_DCM_Analysis/HCM_DCM_October2017/Autosomes/...
#MD5 (Genomes_MegaSeq6.0_RefSeqGeneListAutosomes.dataframe) = f25bca2a8e5a31154ef0e998407feb11
#ALL GENES both locally and server
#MD5 (Genomes_MegaSeq6.0_KnownGene_total.dataframe) = be3f261a4705dd487b951b92c2cdac5d
#sSNVs - Done on server for bootstrap only
#6a0d828259831fb96c2f19b1f305ae81  Genomes_MegaSeq6.0_RefSeqGeneListAutosomes.sSNV.dataframe


if(AllGenes){
  #variantFile.dataframe <- "Genomes_MegaSeq6.0_RefSeqGeneListAutosomes.dataframe" #All Autosomes
  variantFile.dataframe <- "Genomes_MegaSeq6.0_KnownGene_total.dataframe" #All genes
}else if(TTN){
  variantFile.dataframe <- "Genomes_MegaSeq6.0_TTN.dataframe" #TTN variants set
}else if(TTN_PAN){
  variantFile.dataframe <- "Genomes_MegaSeq6.0_Pan_Card_LMM_GeneDx_Invitae_TTN.dataframe" #TTN+ Pan Cardio variants set
}else if(sSNV){ #Done on the server
 variantFile.dataframe <- "Genomes_MegaSeq6.0_RefSeqGeneListAutosomes.sSNV.dataframe"
}else{
  variantFile.dataframe <- "Genomes_MegaSeq6.0_Pan_Card_LMM_GeneDx_Invitae.dataframe" #Pan Cardio gene set
}

input.all <- read.csv(variantFile.dataframe, sep="\t",colClasses = "character") #dataframe containing all the data


##Convert variables to the types 
#Factors 
input.all$Id=as.factor(input.all$Id)
input.all$phen=as.factor(input.all$phen)
input.all$Gene=as.factor(input.all$Gene)
input.all$Chr=as.factor(input.all$Chr)
input.all$Ref=as.character(input.all$Ref)
input.all$Alt=as.character(input.all$Alt)
#Ordered factors
input.all$Eff = factor(input.all$Eff,levels=c("L","M","H"),ordered=TRUE) #Ordered! added low 9/23/2019 LLP, NU

#Conversion of numbers, notice that some of them will produce NA because there is no number for them in the dataframe
input.all$Loc=as.integer(input.all$Loc)
#input.all$gnomAD_AF_NFE_E=as.numeric(input.all$gnomAD_AF_NFE_E)
#input.all$gnomAD_AF_AFR_E=as.numeric(input.all$gnomAD_AF_AFR_E)
#input.all$gnomAD_AF_E=as.numeric(input.all$gnomAD_AF_E)
#input.all$gnomAD_COV_E=as.numeric(input.all$gnomAD_COV_E)
input.all$PP2=as.numeric(input.all$PP2)
#input.all$ExAC=as.numeric(input.all$ExAC)


# We don't convert GERP because it has "Indel" information
#input.all$GERP=as.numeric(input.all$GERP)

#Split GTEX into numerical columns NOTE that the number of GTEx columns changes with sets 
input.all <- input.all %>% extract(GTEx,into=c("GTExLV","GTExAA"),"LV,(.+),AA,(.+)") 
# after version 6p data is structured differently: median, max, min, standard deviation
# and for some reason names have been changed by someone who can't think straight.
#HAA,25.37,150.8,6.55,21.269006,HLV,12.89,61.43,1.911,8.469271
#input.all <- input.all %>% extract(GTEx,into=c("GTExAA","GTExLV"),"HAA,([^,]+),.+HLV,([^,]+),.*") 

#Create actual frequencies for the 900 genomes, HCM & DCM
input.all <- input.all %>%
 separate(NUGENOMES_FREQ, c("num","den"), sep = "\\/") %>%
 mutate(NU_FREQ = as.numeric(num)/as.numeric(den)) %>%
 select(-num,-den)

if("DCM_FREQ" %in% names(input.all)){
input.all <- input.all %>%
 separate(DCM_FREQ, c("num","den"), sep = "\\/") %>%
 mutate(DCM_FREQ = as.numeric(num)/as.numeric(den)) %>%
 select(-num,-den)

#Changed because the cardiomyopathy and whole gene set have different names 
input.all <- input.all %>%
 separate(HCM_FREQ, c("num","den"), sep = "\\/") %>%
 mutate(HCM_FREQ = as.numeric(num)/as.numeric(den)) %>%
 select(-num,-den)

}else{
input.all <- input.all %>%
 separate(DCMFreq, c("num","den"), sep = "\\/") %>%
 mutate(DCM_FREQ = as.numeric(num)/as.numeric(den)) %>%
 select(-num,-den)

input.all <- input.all %>%
 separate(HCMFreq, c("num","den"), sep = "\\/") %>%
 mutate(HCM_FREQ = as.numeric(num)/as.numeric(den)) %>%
 select(-num,-den)
}

# NOTE THAT THE NAME PANCARDIO REMAINED to indicate THE ANALYSIS SET
#Create filtered dataset, note that it assumes that only one type of input gene set can be selected
runTitle <- ""
if(gnomAD == T & NU_FQ == F){
   panCardio <- input.all[ FqMin <= input.all$gnomAD_AF_E  & input.all$gnomAD_AF_E <= FqMax & ! is.na(input.all$gnomAD_AF_E) ,]
   runTitle <- paste('gnomAD::',FqMin ,' <= gnomAD_E <= ', FqMax, sep = "")
}else if(gnomAD == F & NU_FQ == T){
   panCardio <- input.all[ FqMin <= input.all$NU_FREQ  & input.all$NU_FREQ <= FqMax  ,]
   runTitle <- paste('NU_GENOMES_FQ::',FqMin ,' <= NU_GENOMES_FQ <= ', FqMax, sep = "")
}else{
   panCardio <- input.all
   runTitle <- "Freq::No_filter"
}

rm(input.all) # Eliminate the raw input


if (MultiAllelic == F) {
 panCardio$AF <- as.numeric(panCardio$AF) # Multiallelic are coereced to NA because they contain double numbers
 panCardio <-  panCardio[! is.na(panCardio$AF),]
 runTitle <- paste(runTitle, " Multiallelic::removed", sep = "")
}else{
 runTitle <- paste(runTitle, " Multiallelic::included", sep = "")
}

if (ALMS1 == F) {
panCardio <-  panCardio[! panCardio$Gene == "ALMS1",]
 runTitle <- paste(runTitle, " ALMS1::removed", sep = "")
}else{
 runTitle <- paste(runTitle, " ALMS1::included", sep = "")
}


if (Indels == F) {
 panCardio <- panCardio %>% filter ( nchar(panCardio$Ref) == 1 & nchar(panCardio$Alt) ==1)
 #panCardio <-  panCardio[! panCardio$GERP == "INDEL",]
 runTitle <- paste(runTitle, " Indels::removed", sep = "")
}else{
 runTitle <- runTitle
}


if ( Remove_CBL == T){
 panCardio <-  panCardio[! panCardio$Gene == "CBL",]
 runTitle <- paste(runTitle, " CBL removed", sep = "")
}

if(AllGenes){
 runTitle=paste("Genes::All",runTitle) 
}else if (TTN){
 runTitle=paste("Genes::TTN",runTitle) 
}else if (TTN_PAN){
 runTitle=paste("Genes::TTN+PanCardio",runTitle) 
}else if(sSNV){ #Done on the server
 runTitle=paste("Genes::sSNVs",runTitle) 
}else{
 runTitle=paste("Genes::PanCardio",runTitle) 
}

