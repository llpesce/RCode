## Create tables counting variants belonging to different groups (of IDs or variants)
## Functions that produce dataframes to be analyzed by various functions 

## Tables

variantsBySubject <- panCardio %>% filter(Chr!="chrX") %>% group_by (Id) %>% summarise(N=n())
summary(variantsBySubject)

IdCount <-subjPheno[subjPheno$Id %in% panCardio$Id,] %>% group_by(Race_Recode_Cluster2,type.of.CM) %>% summarize ('Tot'=n())


#Build table with counts by group
Count <- merge(panCardio, subjPheno, by = 'Id') %>% group_by(Gene,Race_Recode_Cluster2,type.of.CM,Eff) %>% summarize ('Tot'=n()) %>% spread(Eff, Tot, fill= 0 )

check_df <- function(df){
    print(paste("Nr of variants:",nrow(df)))
    print(paste("Nr of Genes:",df %>% select(Gene) %>% unique %>% nrow))
}

 
#Attach blocks of variants only by gnomAD_AF_E bucket

## Functions 

#Create count dataframe 

create.count.df <- function ( dataframe) 
{
 # Variant dataframe with multiple rows with different types of variants

  #start with all variants in the dataframe
  A <- as.data.frame(table( dataframe  %>% select(Id) ))
  names(A) <- c("Id", "V.all")

  #Add the variants that have a gnomAD value
  M <- as.data.frame(table(dataframe %>% filter( ! is.na(gnomAD_AF_E) ) %>% select(Id) %>% select(Id) ))
  names(M) <- c("Id", "V.gnomAD")
  #Merge
  A <- merge(A,M,by="Id",suffixes="")

  #Filter by GERP & PP2 and Eff
  M <- as.data.frame(table( dataframe %>% filter(PP2 >= .9 | GERP == "INDEL" | as.numeric(GERP) >= 3.0 | Eff == 'H' ) %>% select(Id) ))
  names(M) <- c("Id", "V.PP2_9_GERP3_Eff_H")
  A <- merge(A,M,by="Id",suffixes="")

  M <- as.data.frame(table(dataframe %>% filter(gnomAD_AF_E <= .0025) %>% select(Id) ))
  names(M) <- c("Id", "V.Freq.0025")
  A <- merge(A,M,by="Id",suffixes="")

  M <- as.data.frame(table(dataframe %>% filter(gnomAD_AF_E <= .01) %>% select(Id) ))
  names(M) <- c("Id", "V.Freq.01")
  A <- merge(A,M,by="Id",suffixes="")

  M <- as.data.frame(table(dataframe %>% filter(gnomAD_AF_E <= .1) %>% select(Id) ))
  names(M) <- c("Id", "V.Freq.1")
  A <- merge(A,M,by="Id",suffixes="")

  M <- as.data.frame( table(dataframe %>% filter( gnomAD_AF_E <= .25 ) %>% select(Id) ))
  names(M) <- c("Id", "V.Freq.25")
  A <- merge(A,M,by="Id",suffixes="")

  M <- as.data.frame( table(dataframe %>% filter( gnomAD_AF_E > .1 & gnomAD_AF_E <= .25 ) %>% select(Id) ))
  names(M) <- c("Id", "V.Freq.1_.25")
  A <- merge(A,M,by="Id",suffixes="")

  M <- as.data.frame( table(dataframe %>% filter( gnomAD_AF_E > .25 & gnomAD_AF_E <= .5 ) %>% select(Id) ))
  names(M) <- c("Id", "V.Freq.25_.5")
  A <- merge(A,M,by="Id",suffixes="")

  return(A)
}

create.norm.df <-  function (df1,df2,var)
{
  #Return ratio dataframe with new label
  X <- merge( df1[,c("Id", var)],df2[,c("Id", var)],by="Id") 
  ret_names <- names(X)
  names(X) <- c("Id","first","second")
  X <- X %>% mutate( ratio = first/second) %>% select(Id, ratio)
  names(X) <- c("Id",paste(var,"ratio",sep="."))
  return(X)

}

create.GLM.df <- function (A){
  X.GLM <- merge( A,subjPheno[,c("Id","type.of.CM",'Race_Recode_Cluster2','Machine_Recode')],by="Id",suffixes="")
  #recode race as factor
  X.GLM$'Race_Recode_Cluster2' <- factor(X.GLM$'Race_Recode_Cluster2',levels=c("0","1"))
  X.GLM$Machine_Recode <- factor(X.GLM$Machine_Recode, levels=c("1","2"))
  #Fix names
  names(X.GLM) <- names(X.GLM) %>% make.names
  #Order them so that probability is probability of DCM
  X.GLM$type.of.CM <- factor(X.GLM$type.of.CM, levels=c("HCM","DCM"))
  return(X.GLM)
}



