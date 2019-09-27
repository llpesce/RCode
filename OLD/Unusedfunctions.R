

plotSaveGLM <- function(X.GLM,x, name,GLMdir)
{ 
#Function to fit a GLM and plot it for the variable chosen with titles and saving file
    DCM <- X.GLM %>%  filter(type.of.CM == 'DCM') %>% select (Id, x,Race_Recode_Cluster2)
    DCM <- cbind(DCM, c(1))
    names(DCM) <- c("Id","NumVar","Race", "PofDCM")
    #Create matrix with number of subjects for number of variants and use it to shift points
    cnt <- table(DCM$NumVar)
    for (j in names(cnt)){
       shift=0.0
       dx=0.01
       for (Id in DCM[DCM$NumVar==j,]$Id){
          DCM[DCM$Id == Id,]$PofDCM <- DCM[DCM$Id == Id,]$PofDCM+shift
          shift <- shift+dx
       }
    }
    #Points with the averages by race by phenotype
    DCM.0 <- data.frame(mean=mean(DCM[DCM$Race == 0,]$NumVar),loc=.96)
    DCM.1 <- data.frame(mean=mean(DCM[DCM$Race == 1,]$NumVar),loc=.96)

    HCM <- X.GLM %>%  filter(type.of.CM == 'HCM') %>% select (Id,x,Race_Recode_Cluster2)
    HCM <- cbind(HCM, c(0))
    names(HCM) <- c("Id","NumVar","Race", "PofDCM")
    cnt <- table(HCM$NumVar)
    for (j in names(cnt)){
       shift=0.0
       dx=0.01
       for (Id in HCM[HCM$NumVar==j,]$Id){
          HCM[HCM$Id == Id,]$PofDCM <- HCM[HCM$Id == Id,]$PofDCM+shift
          shift <- shift+dx
       }
    }
    #Points with the averages by race by phenotype

    HCM.0 <- data.frame(mean=mean(HCM[HCM$Race == 0,]$NumVar),loc=-.04)
    HCM.1 <- data.frame(mean=mean(HCM[HCM$Race == 1,]$NumVar),loc=-.04)
    
    title<-paste("GLM "," X::", x)
    subtitle<-runTitle

    #Fit GLM and create multilayered plot
    print (paste('for interval',x)) 
    model.GLM <- glm( type.of.CM ~ X.GLM[[x]],family=binomial,data=X.GLM)
    print(summary(model.GLM))
    intercept <- as.numeric( coef( model.GLM )[1] )
    beta <- as.numeric( coef( model.GLM )[2] )
    rangex <- range(X.GLM[[x]])
    range <- seq(max(rangex[1]-10,0),rangex[2],by=1)
    logodds <- intercept + beta * range
    prob <- exp(logodds) / (1+exp(logodds))
    PLOT <- data.frame(range,prob)
    ggplot(PLOT, aes(range)) +                    # basic graphical object
    geom_line(aes(y=prob, colour="Fit")) +  # first layer
    labs(title = title, subtitle = subtitle, x = '# of Variants', y = 'Probability of DCM') + # Third
    geom_point(data=HCM,mapping = aes(x = NumVar, y = PofDCM, colour="HCM",shape=HCM$Race))+
    geom_point(data=DCM,mapping = aes(x = NumVar, y = PofDCM, colour="DCM",shape=DCM$Race))+
    geom_point(data=DCM.0,mapping = aes(x = mean, y = loc, colour="DCM", shape="0", size=2))+
    geom_point(data=DCM.1,mapping = aes(x = mean, y = loc, colour="DCM", shape="1", size=2))+
    geom_point(data=HCM.0,mapping = aes(x = mean, y = loc, colour="HCM", shape="0", size=2))+
    geom_point(data=HCM.1,mapping = aes(x = mean, y = loc, colour="HCM", shape="1", size=2))+
    scale_colour_manual(name="Legend",values=c(Fit="black",HCM="red",DCM="blue",shape=DCM$Race))

    outPng <- paste(GLMdir,"/", x, name,'.png',sep='')

    ggsave(outPng)


}

autosomes.counts <- create.count.df(autosomes.all)
panCardio.counts <- create.count.df(panCardio.all)

imax <- as.integer(ncol(autosomes.counts))
for (i in 2:imax){
  var <- names(autosomes.counts)[i]
  print(var)
  A <- create.norm.df(panCardio.counts,autosomes.counts,var)
  X.GLM <- create.GLM.df(A)
  x <- names(X.GLM)[2]
  model.GLM <- glm( type.of.CM ~ X.GLM[[x]],family=binomial,data=X.GLM)
  print(summary(model.GLM))
}


#build averages
X.GLM <- create.GLM.df(autosomes.counts)
for (i in 2:imax){
    x <- names(X.GLM)[i]
    print(x)
    XTEN  <- X.GLM %>%  filter(Machine_Recode == 1) %>% select (x)
    OTHER <- X.GLM %>%  filter(Machine_Recode == 2) %>% select (x)
    test <- wilcox.test(XTEN[,1],OTHER[,1], conf.int = TRUE, paired=FALSE, exact=TRUE)
    print(test)
    test <- t.test(XTEN[,1],OTHER[,1], conf.int = TRUE, paired=FALSE, exact=TRUE)
    print(test)

}

