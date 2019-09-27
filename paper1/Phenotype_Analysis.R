## Phenotype analysis

# Load the data

# Collect and write out phenotype stats: Note that there is no data for a group 
subjPheno %>% group_by(Race_Recode_Cluster2,Machine_Recode,type.of.CM) %>% summarise('<EF>'=mean(EF,na.rm=TRUE),'<LVPWd_BSA>'=mean(LVPWd_BSA,na.rm=TRUE),'<LVIDd_BSA>'=mean(LVIDd_BSA,na.rm=TRUE),'<IVSd_BSA>'=mean(IVSd_BSA,na.rm=TRUE))
#  Race_Recode_Cluste… Machine_Recode type.of.CM `<EF>` `<LVPWd_BSA>` `<LVIDd_BSA>` `<IVSd_BSA>`
#                <dbl>          <dbl> <chr>       <dbl>         <dbl>         <dbl>        <dbl>
#1                   0              1 DCM          28.0         0.527          3.15        0.538
#2                   0              1 HCM          63.5         0.661          2.24        1.04 
#3                   0              2 DCM          28.3         0.524          3.40        0.516
#4                   0              2 HCM          63.4         0.595          2.33        0.886
#5                   1              1 DCM          24.1         0.574          3.43        0.533
#6                   1              1 HCM         NaN         NaN            NaN         NaN    
#7                   1              2 DCM          13.7         0.470          4.20        0.479
#8                   1              2 HCM          68.2         0.679          2.56        0.718


library("ggplot2")
library(ggrepel)
(panCardio %>% filter( Id != '5167') %>% select(Id) %>% unique )$Id -> ID # vector with the 127 IDs in the HCM/DCM paper minus the juvenile with Rheumathoid arthritis 
 length(ID)
#[1] 126


outDir="~/Dropbox/DCM_HCM_Landscape_MRP/PhenotypePlots"
dir.create(outDir)

##PCA analysis 

#Filter the IDs that are in the dataframe and also kept
#QTc is missing from too many cases
PCA.df <- subjPheno %>% filter( Id %in% ID) %>% 
   select('Id','type.of.CM','Race_Recode_Cluster2', 'Race.ethnicity', 'Sex', 'Machine_Recode', #factor variables
   'IVSd_BSA','LVPWd_BSA','LVIDd_BSA','EF') # quantitative variables


PCA.df[complete.cases(PCA.df),] -> PCA.df
nrow(PCA.df)
#[1] 89
#[1] 82 #without 6167
PCA.df$type.of.CM <- as.factor(PCA.df$type.of.CM)
PCA.df$Race_Recode_Cluster2 <- factor(PCA.df$Race_Recode_Cluster2, 
                           levels = c("0", "1"), 
                           labels = c("EA", "AA"))
PCA.df$Sex <- as.factor(PCA.df$Sex)
PCA.df$Machine_Recode <- factor(PCA.df$Machine_Recode, 
                           levels = c("1", "2"), 
                           labels = c("XTEN", "Other"))

summary(PCA.df)
#      Id     type.of.CM Race_Recode_Cluster2 Race.ethnicity Sex    Machine_Recode    IVSd_BSA        LVPWd_BSA        LVIDd_BSA           EF       
# 4001   : 1   DCM:38     EA:70                AA :12         F:36   XTEN :39       Min.   :0.2979   Min.   :0.3379   Min.   :1.375   Min.   : 4.00  
# 4003   : 1   HCM:44     AA:12                EA :58         M:46   Other:43       1st Qu.:0.5000   1st Qu.:0.4940   1st Qu.:2.238   1st Qu.:25.10  
# 4004   : 1                                   HA: 9                               Median :0.6642   Median :0.5590   Median :2.635   Median :55.00  
# 4005   : 1                                   OA: 3                               Mean   :0.7426   Mean   :0.5744   Mean   :2.774   Mean   :47.28  
# 4009   : 1                                                                        3rd Qu.:0.8820   3rd Qu.:0.6478   3rd Qu.:3.142   3rd Qu.:65.00  
# 4012   : 1                                                                        Max.   :2.2330   Max.   :0.8389   Max.   :4.778   Max.   :81.00  
# (Other):76                         

#Merge OA into EA to avoid diluting power or relative small differences in ancestra origin
PCA.df$Race.ethnicity <- PCA.df$Race.ethnicity %>% fct_collapse(EA = c("EA","OA")) 
summary(PCA.df$Race.ethnicity)
#AA EA HA 
#12 61  9 


#Collect the numerical parts for the PCA
X.pca <- PCA.df %>% select(IVSd_BSA,LVPWd_BSA,LVIDd_BSA,EF) # Leave only the numeric data
f2 <- prcomp(as.matrix(X.pca), scale=T, center=T)
summary(f2)
#Importance of components:
#                          PC1    PC2    PC3     PC4
#Standard deviation     1.5341 1.0021 0.6605 0.45411
#Proportion of Variance 0.5883 0.2511 0.1091 0.05155
#Cumulative Proportion  0.5883 0.8394 0.9485 1.00000

plot(f2) # Scree plot
biplot(f2,xlabs=rep("·", nrow(X.pca)), cex=c(5, 1), xlim=c(-.3,+.3))

# To print as png
# dev.copy(png,'Scree.PCA.png')
# dev.off()


##Better plot made with ggplot
#Change pallette
#palette.pca <- c( "#3366CC","#339900") # Colors for DCM and HCM x race Before August 2019
palette.pca <- c( "#000000","#FF0000") # Colors for DCM and HCM x race

pca <- cbind(
  PCA.df %>% select(-IVSd_BSA,-LVPWd_BSA,-LVIDd_BSA,-EF), #Eliminate numerical variables
  predict(f2,X.pca)[,c(1,2)] #Add predicted PC1 and PC2 
 ) 
pca$PC2 <- -pca$PC2 #Reverse the direction of PC1, so that DCM is for larger values

#Plot
g <-  ggplot(pca,aes(x=PC1,y=PC2,color=type.of.CM))+
    geom_point(size=5) +
scale_colour_manual(values=palette.pca)

#Axe legend and plot
g +theme_bw() + theme(legend.position="none")

outPng <- paste(outDir,'/PCA.png',sep='')
ggsave(outPng)

# Make pca plot with Race and sex

palette2.pca <- c( "#000000","#A0A0A0","#FF0000","#FA8072") # Colors for DCM and HCM x race as Race_Recode_Cluster2
palette4.pca <- c( "#CCCCCC","#999999","#000000","#CCCC99","#FFCCCC","#FF6666","#FF0000","#CC6600") # Colors for DCM and HCM x race as Collapsed Race.ethnicity with 4
#palette3.pca <- c( "#CCCCCC","#999999","#000000","#FFCCCC","#FF6666","#FF0000") # Colors for DCM and HCM x race as Collapsed Race.ethnicity with 3
#palette3.pca <- c( "#666666","#000000","#999999","#FFCCCC","#FF6666","#FF0000") # Colors for DCM and HCM x race as Collapsed Race.ethnicity with 3
palette3.pca <- c( "#666666","#000000","#CCCCCC","#FF6600","#FF0000","#FF9900") # Colors for DCM and HCM x race as Collapsed Race.ethnicity with 3 order is : DCM -> AA,EA,HA ; HCM same
palette3_2.pca <- c( "#000000","#000000","#000000","#FF0000","#FF0000","#FF0000") # Colors for DCM and HCM x race as Collapsed Race.ethnicity with 3 order is : DCM -> AA,EA,HA ; HCM same

 
g <- ggplot(pca,aes(x=PC1,y=PC2,color=type.of.CM:Race.ethnicity,shape=Sex)) +
  geom_point(size=5) +
  scale_colour_manual(values=palette3.pca) +
  theme_bw() 

# Add an easy to remove, but helpful for plots legend
g + theme(legend.position="bottom") + theme(legend.title = element_blank())

outPng <- paste(outDir,'/PCA.Sex.Race.png',sep='')
ggsave(outPng)

##Regress first PCA coordinate on nsSNVs

#Build dataframe with PCAs and total variant counts
X <-merge( counts %>% filter(! Id == '5167'), pca,by="Id",suffixes="")

Y='Echo PC1'

# Fit PC1 vs V.all, marginalize over all the rest
model.0 <- lm (PC1 ~ V.all,X)
#Only EU, according to Race_Recode_Cluster2
model.onlyEA.RRC2 <- lm (PC1 ~ V.all,X %>% filter(Race_Recode_Cluster2 =='EA'))
#Only EU, according to self declared Race plus Race_Recode_Cluster when didn't fit EU/AA/HA 
model.onlyEA.RE <- lm (PC1 ~ V.all,X %>% filter( Race.ethnicity=='EA'))
#Everyone, excluding HA
model.noHA.RE <- lm (PC1 ~ V.all,X %>% filter( Race.ethnicity !='HA'))
#Everyone, excluding HA & OA, if OA is present (if it isn't of course the logical will have no effect
model.noHAO.RE <- lm (PC1 ~ V.all,X %>% filter( Race.ethnicity !='HA' & Race.ethnicity !='OA'))


#Sink the fits and reset the fits file
sink(paste(outDir,"/","pheno_fits.txt",sep=''))
summary(model.0)
summary(model.onlyEA.RRC2)
summary(model.onlyEA.RE)
summary(model.noHA.RE)
summary(model.noHAO.RE)
sink()

#Predicted coordinates
pred.0 <-  data.frame(Y.pred = predict(model.0, X),var=X$V.all)


# Plot only with HCM and DCM
ggplot(X,aes(x=V.all,y=PC1,color=type.of.CM))+  
    geom_point(size=6, show.legend=F) + 
    labs(title = paste('Scatterplot of ',Y,'by type of CM'), x = 'V.all', y = 'PC1') +
    geom_line(color='black',data = pred.0,size=2, linetype=6, aes(x=var, y=Y.pred)) +
    scale_colour_manual(values=palette.pca) +theme_bw()+
    theme(legend.position="bottom") + theme(legend.title = element_blank())


outPng <- paste(outDir,'/PC1VsV.all.png',sep='')
ggsave(outPng)


#Plot with Race.ethnicity and then Sex
pred.EA.RRC2 <-  data.frame(Y.pred = predict(model.onlyEA.RRC2, X %>% filter(Race_Recode_Cluster2 =='EA') ),var=(X%>% filter(Race_Recode_Cluster2 =='EA'))$V.all)
pred.EA.RE <-  data.frame(Y.pred = predict(model.onlyEA.RE, X %>% filter(Race.ethnicity=='EA') ),var=(X%>% filter(Race.ethnicity=='EA'))$V.all)
pred.noHA.RE <-  data.frame(Y.pred = predict(model.noHA.RE, X %>% filter(Race.ethnicity!='HA') ),var=(X%>% filter(Race.ethnicity!='HA'))$V.all)
pred.noHAO.RE <-  data.frame(Y.pred = predict(model.noHAO.RE, X %>% filter(Race.ethnicity!='HA' & Race.ethnicity!='OA')  ),var=(X%>% filter(Race.ethnicity!='HA' & Race.ethnicity!='OA') )$V.all)


#Add lines, switch between palette3 and four and remove line depending upon what to do with OA.
# With different colors by ethnicity and condition
ggplot(X,aes(x=V.all,y=PC1,color=type.of.CM:Race.ethnicity))+  
     geom_point(size=6, show.legend=T) + 
     labs(title = paste('Scatterplot of ',Y,'by type of CM'), x = 'V.all', y = 'PC1') +
     geom_line(color='blue',data = pred.0,size=1, linetype=1, aes(x=var, y=Y.pred)) +
     geom_line(color=palette3.pca[2],data = pred.EA.RRC2,size=1, linetype=5, aes(x=var, y=Y.pred)) +
     geom_line(color=palette3.pca[2],data = pred.EA.RE,size=1, linetype=6, aes(x=var, y=Y.pred)) +
     geom_line(color='blue',data = pred.noHA.RE,size=1, linetype=6, aes(x=var, y=Y.pred)) +
     scale_colour_manual(values=palette3_2.pca) +theme_bw()+
     theme(legend.position="bottom") + theme(legend.title = element_blank())
outPng <- paste(outDir,'/PC1VsV.all.various.race.fits.png',sep='')


#Just two colors
ggplot(X,aes(x=V.all,y=PC1,color=type.of.CM:Race.ethnicity))+  
     geom_point(size=6, show.legend=T) + 
     labs(title = paste('Scatterplot of ',Y,'by type of CM'), x = 'V.all', y = 'PC1') +
     geom_line(color='gray',data = pred.0,size=1, linetype=1, aes(x=var, y=Y.pred)) +
     geom_line(color='gray',data = pred.noHA.RE,size=1, linetype=6, aes(x=var, y=Y.pred)) +
     geom_line(color='gray',data = pred.EA.RRC2,size=1, linetype=2, aes(x=var, y=Y.pred)) +
     geom_line(color='gray',data = pred.EA.RE,size=1, linetype=3, aes(x=var, y=Y.pred)) +
     scale_colour_manual(values=palette3_2.pca) +theme_bw()+
     theme(legend.position="bottom") + theme(legend.title = element_blank())
outPng <- paste(outDir,'/PC1VsV.all.racefits.noRace.png',sep='')


ggsave(outPng)

# Race_Recode_cluster2 & Sex
ggplot(X,aes(x=V.all,y=PC1,color=type.of.CM:Race_Recode_Cluster2,shape=Sex)) +  
    geom_point(size=6, show.legend=T) + 
    labs(title = paste('Scatterplot of ',Y,'by type of CM'), x = 'V.all', y = 'PC1') +
    geom_line(data = pred.0,  mapping= aes(x=var,y=Y.pred),color='black',size=2,linetype=6,inherit.aes=FALSE) +
    scale_colour_manual(values=palette2.pca) +theme_bw()+
    theme(legend.position="bottom") + theme(legend.title = element_blank())

outPng <- paste(outDir,'/PC1VsV.all.raceNsex.png',sep='')
ggsave(outPng)

##Checking machine and PC1

Y='Echo PC1'


model.m <- lm (PC1 ~ Machine_Recode+ V.all,X)
model.r <- lm (PC1 ~ Race_Recode_Cluster2 + V.all,X)
model.mr <- lm (PC1 ~ Machine_Recode + Race_Recode_Cluster2 + V.all,X)

#Add to the fits file
sink(paste(outDir,"/","pheno_fits.txt",sep=''), append = TRUE)
 
summary(model.m)
summary(model.r)
summary(model.mr)

#DO nested model ANOVA
anova(model.0, model.m)
anova(model.0, model.r)
anova(model.m, model.mr)
sink()

#Predicted coordinates including machine
pred <-  data.frame(Y.pred = predict(model.m, X),V.all=X$V.all,Machine_Recode=X$Machine_Recode)

# Plot only with HCM and DCM
ggplot(X,aes(x=V.all,y=PC1,color=type.of.CM,shape=Machine_Recode))+  
    geom_point(size=6, show.legend=T) + 
    scale_shape_manual(values=c(15,0))+
    labs(title = paste('Scatterplot of ',Y,'by type of CM'), x = 'V.all', y = 'PC1') +
    geom_line(data=pred,aes(x=V.all,y=Y.pred,linetype=Machine_Recode),color='gray',size=1,inherit.aes=FALSE) +
    geom_line(data=pred.0,aes(x=var,y=Y.pred),linetype=6,color='black',size=2,inherit.aes=FALSE) +
    scale_colour_manual(values=palette.pca) +theme_bw()+
    theme(legend.position="bottom") + theme(legend.title = element_blank())

outPng <- paste(outDir,'/PC1VsV.all.A.Machine.png',sep='')
ggsave(outPng)


# Same as before for technicity, but add machine
#Only EU, according to Race_Recode_Cluster2
model.onlyEA.m.RRC2 <- lm (PC1 ~ Machine_Recode + V.all,X %>% filter(Race_Recode_Cluster2 =='EA'))
#Only EU, according to self declared Race plus Race_Recode_Cluster when didn't fit EU/AA/HA 
model.onlyEA.m.RE <- lm (PC1 ~ Machine_Recode+V.all,X %>% filter( Race.ethnicity=='EA'))
#Everyone, excluding HA
model.noHA.m.RE <- lm (PC1 ~ Machine_Recode+V.all,X %>% filter( Race.ethnicity !='HA'))
#Everyone, excluding HA & OA
model.noHAO.m.RE <- lm (PC1 ~ Machine_Recode+V.all,X %>% filter( Race.ethnicity !='HA' & Race.ethnicity !='OA'))

sink(paste(outDir,"/","pheno_fits.txt",sep=''), append = TRUE)
summary(model.m)
summary(model.onlyEA.m.RRC2)
summary(model.onlyEA.m.RE)
summary(model.noHA.m.RE)
summary(model.noHAO.m.RE)
#Makes no material difference

sink()




#Other stuff tested for various reviewers and for fun

## checking age at echo measurements to see if it has anything to do with the experiments
##Plot of age vs nsSNVs 

X <-left_join(X, subjPheno %>% select(Id, years.old),by="Id",suffixes="")

Y='Age at echo'

model <- lm (years.old ~ V.all,X)
pred <-  data.frame(Y.pred = predict(model, X),var=X$V.all)


ggplot(X,aes(x=V.all,y=years.old,color=type.of.CM:Race_Recode_Cluster2,shape=Sex)) +  
    geom_point(size=6, show.legend=F) + 
    labs(title = paste('Scatterplot of ',Y,'by type of CM'), x = 'V.all', y = 'PC1') +
    geom_line(data = pred,  mapping= aes(x=var,y=Y.pred),color='black',size=2,linetype=6,inherit.aes=FALSE) +
    scale_colour_manual(values=palette2.pca) +theme_bw()+ theme(legend.position="none")

outPng <- paste(outDir,'/AgeAtEchoVsV.all.raceNsex.png',sep='')
ggsave(outPng)


#Predict age from PC1, to see if PC1 is just telling us about age

model <- lm (years.old ~ PC1,X)
pred <-  data.frame(Y.pred = predict(model, X),var=X$PC1)

ggplot(X,aes(x=PC1,y=years.old,color=type.of.CM:Race_Recode_Cluster2,shape=Sex)) +  
    geom_point(size=6, show.legend=F) + 
    labs(title = paste('Scatterplot of ',Y,'by type of CM'), x = 'PC1', y = 'years.old') +
    geom_line(data = pred,  mapping= aes(x=var,y=Y.pred),color='black',size=2,linetype=6,inherit.aes=FALSE) +
    scale_colour_manual(values=palette2.pca) +theme_bw()+ theme(legend.position="none")

outPng <- paste(outDir,'/AgeAtEchoVsPC1.raceNsex.png',sep='')
ggsave(outPng)





# Regress variables as a function of each other
#x.name<-'EF'
#x.name<-'largest.LVIDd'
#x.name <- 'Largest.ivsd.value.cm.'
x.name <- 'Largest.LVPWd.cm.'
#y.name<-'LVIDd_BSA'
#y.name <- 'IVSd_BSA'
y.name <- 'LVPWd_BSA'
#y.name<-'bsa.m2'

xy.df <- subjPheno[,c('Id','type.of.CM','Sex',x.name,y.name)] 
xy.df$type.of.CM <-as.factor(xy.df$type.of.CM )
xy.df <- xy.df[complete.cases(xy.df),]

xy.df <- xy.df %>% filter (Id %in% ID) #Only selected IDs are kept

cor(xy.df[[y.name]],xy.df[[x.name]])
model <- lm(data=xy.df,xy.df[[y.name]] ~ xy.df[[x.name]])
summary(model)

xy.df %>% group_by(type.of.CM,Sex) %>% summarise_at(c(x.name,y.name),funs(mean,sd))

#basic scatterplot layer
p <- ggplot(xy.df, aes(x=xy.df[,x.name],y=xy.df[,y.name],col=type.of.CM,shape=Sex))+
geom_point(size=5) +
labs(title = paste('Scatterplot'),x = x.name, y = y.name)+
geom_smooth(method=lm, size=2) +
scale_colour_manual(values=palette) 
#Add type.of.CM to plot
model1 <- lm(data=xy.df,xy.df[[y.name]] ~ type.of.CM + xy.df[[x.name]])
summary(model1)

p

outPng <- paste(outDir,'/X_', x.name,'_Y_',y.name,'.png',sep='')
ggsave(outPng)


# To Add points for LVIDd & EF vs bsa
#p + geom_text_repel(data=xy.df, aes(x=xy.df[,x.name],y=xy.df[,y.name],label=ifelse( xy.df[,y.name] < 1.53 | #xy.df[,y.name] > 2.3, as.character(Id),'')))
# add points for LVIDd vs. EF
#p1 <- p + geom_text_repel(data=xy.df, aes(x=xy.df[,x.name],y=xy.df[,y.name],label=ifelse( xy.df$Id=='5173' |  xy.df$Id =='5115' | xy.df$Id=='4433' | xy.df$Id =='5161' | xy.df$Id=='5171' | xy.df$Id=='5238' | xy.df$Id=='5168'| xy.df$Id=='4013'| xy.df$Id=='5142'| xy.df$Id=='5123' | xy.df$Id=='4456'| xy.df$Id=='4008', as.character(Id),'')))


#Density plots over coordinates
hist(pca$PC1, seq(-4,4,.4))
d <- density(pca$PC1, bw=0.4)
d$y <- d$y*10/max(d$y)
lines(d,lwd=5)
d.HCM <- density(PC1.HCM, bw=0.4)
d.HCM$y <- d.HCM$y*10/max(d.HCM$y)
lines( d.HCM,col="#339900",lwd=5 )
d.DCM <- density(PC1.DCM, bw=0.4)
d.DCM$y <- d.DCM$y*8.8/max(d.DCM$y)
lines( d.DCM,col="#3366CC",lwd=5 )

dev.copy(png,'PCA1.hist.png')
dev.off()


