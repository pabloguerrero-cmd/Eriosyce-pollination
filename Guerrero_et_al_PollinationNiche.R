##Title: Floral mimicry induces sympatric speciation in a South American cactus 
##Authors: Pablo C. Guerrero, Claudia A. Antinao, Mary T. K. Arroyo, Deren Eaton, Beatriz Vergara–Meriño, Heidy Villalobos, Gastón O. Carvallo

## Check packages
packages = c("vegan", "ggplot2",
             "png", "grid","lme4","nlstimedist", 
             "nlstools","abind","plyr","ggridges","emmeans", "dplyr","boot")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

##setup directory
setwd("~/data")

##MORPHOSPACE##
##MULTIDIMENSIONAL SCALING (NMDS) 
datos_nmds<-read.table(file="1morphospace.csv", header=T, sep=";", check.names=T, dec=".", row.names=1) 
ord<-read.table(file="1morphospace.csv", header=T, sep=";", check.names=T, dec=".") 
ordenamiento<-ord[,c(1,2)]
attach(ordenamiento)

nmds_traits<-datos_nmds[,c(3:13)] #model01
model01<-metaMDS(nmds_traits, distance="euclidean", binary=F, k=2, zerodist="add", trymax=1000, autotransform=T,plot=F,na.rm=T) 
model01$stress
stressplot(model01) # Produces a Shepards diagram

#Multidimensional analysis
dist_traits<-vegdist(nmds_traits, method="euclidean", binary=F)
adonis_traits<-adonis(dist_traits ~ Taxa, ordenamiento,method="euclidean", binary=F, diag=F, permutations=999)
adonis_traits

# NMDS Plot
scores<-scores(model01, display="sites", shrinf=F)
scores<-cbind(scores,ordenamiento)

temp=c(2,2,2,2)
barplot(temp,col=c("#189E76BF","#DA5D00BF","#736DB3BF","#E7AB00BF"), names.arg=c("chilensis", "mutabilis", "litoralis","albidiflora"))

plot01<-ggplot(scores,aes(x = NMDS1, y = NMDS2), color=Taxa,shape=Taxa) +
  theme_gray(base_size=30)

plot02<-plot01 +stat_ellipse(geom = "polygon", aes(group = Taxa, color = Taxa, fill = Taxa), 
                     alpha = 0.3,show.legend =F)+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),
        axis.text=element_text(size=22), 
        panel.border=element_rect(colour = "black", fill=NA, size=3),
        legend.position="right", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.title = element_blank(),
        legend.text=element_text(size=18))+
  geom_point(aes(fill=Taxa),colour="black",shape=21, size=4)+
  guides(col=FALSE)+
  scale_colour_manual(values = c("#E7AB00BF","#189E76BF", "#736DB3BF","#DA5D00BF"),
                      labels = element_blank()) +
  scale_fill_manual(values = c("#E7AB00BF","#189E76BF","#736DB3BF","#DA5D00BF"),
                    labels=c("chilensis-albidiflora", "chilensis", "litoralis","mutabilis"))
plot02

##PHENOLOGY##
# This code compute nonlinear models describing pollination niche. Niche similarity is based on the estimated variables (days de apertura, flowering span,flower50%, r, t and c)
# Citation: Steer et al (2019) https://doi.org/10.1111/2041-210X.13293

phenology<-read.table(file="2fenologia.csv", header=T, sep=";", check.names=T, dec=".")

# r: maximum growth rate (1/flowering time). Flowering time represents the whole period of phenology.
# c: concentration of cases, number of cases when r becomes maximum. 
# t: measure of the weight of the process duration [time].

# Phenology_albidiflora
alb<-phenology[c(1:23),c(5,14)]
dt_alb<-tdData(alb,x="dfs_corr", y="alb")
model_alb<-timedist(dt_alb,x = "dfs_corr", y = "propMax", 
                    r = 0.011, c = 0.5, t = 42)
attach(model_alb)
confint2(model_alb) #confidence intervals
tdRSS(model_alb) # pseudo-R
tdPdfPlot(model_alb, S = 1, xVals = seq(0, 160, 0.01)) #Probability of flowering of all plants
tdCdfPlot(model_alb, S = 1, xVals = seq(0, 160, 0.01))

# Phenology_chilensis
chi<-phenology[c(1:23),c(5,8)]
chi
dt_chi<-tdData(chi,x="dfs_corr", y="chi_Molles")
model_chi<-timedist(dt_chi,x = "dfs_corr", y = "propMax", 
                    r = 0.0119, c = 0.5, t = 42)
confint2(model_chi) #confidence intervals
tdRSS(model_chi) # pseudo-R
tdPdfPlot(model_chi, S = 1, xVals = seq(0, 160, 0.01))
tdCdfPlot(model_chi, S = 1, xVals = seq(0, 160, 0.01))

# Phenology_litoralis (includes Pichidangui and Chivato, named "Pichidangui")
litPich<-phenology[c(1:23),c(5,15)]
dt_litPich<-tdData(litPich,x="dfs_corr", y="lit")
model_litPich<-timedist(dt_litPich,x = "dfs_corr", y = "propMax", 
                        r = 0.007, c = 0.5, t = 63)
model_litPich
confint2(model_litPich) #confidence intervals
tdRSS(model_litPich) # pseudo-R
tdPdfPlot(model_litPich, S = 1, xVals = seq(0, 160, 0.01))
tdCdfPlot(model_litPich, S = 1, xVals = seq(0, 160, 0.01))

#Phenology_mutabilis
mut<-phenology[c(1:23),c(5,12)]
dt_mut<-tdData(mut,x="dfs_corr", y="mut_Qchiv")
model_mut<-timedist(dt_mut,x = "dfs_corr", y = "propMax", 
                    r = 0.0243, c = 0.0155, t = 20.5)
confint2(model_mut,n = c(0.05,0.95)) #confidence intervals
tdRSS(model_mut) # pseudo-R
tdPdfPlot(model_mut, S = 1, xVals = seq(0, 160, 0.01))
tdCdfPlot(model_mut, S = 1, xVals = seq(0, 160, 0.01))

# All models polot
tdPdfPlot(model_alb, model_chi,model_litPich,model_mut, 
          S = c(0.9, 0.9, 0.9,0.9),
          xVals = seq(-20, 180, 0.1))
tdCdfPlot(model_alb,model_chi,model_litPich,model_mut, 
          S = c(0.9, 0.9, 0.9,0.9),xVals = seq(0, 180, 0.1))#,yVals=seq(0,1,1))

## Edited Plot_phenology (95% confidence intervals)
##parameters
dt_alb$pred<-predict(model_alb) 
dt_mut$pred<-predict(model_mut)
dt_chi$pred<-predict(model_chi)
dt_litPich$pred<-predict(model_litPich)
# estandard error 
se_alb = summary(model_alb)$sigma
se_mut = summary(model_mut)$sigma
se_chi = summary(model_chi)$sigma
se_lit = summary(model_litPich)$sigma
# CI95%
ci_alb = outer(dt_alb$pred, c(outer(se_alb, c(-1,1), '*'))*1.96, '+')
ci_mut = outer(dt_mut$pred, c(outer(se_mut, c(-1,1), '*'))*1.96, '+')
ci_chi = outer(dt_chi$pred, c(outer(se_chi, c(-1,1), '*'))*1.96, '+')
ci_lit = outer(dt_litPich$pred, c(outer(se_lit, c(-1,1), '*'))*1.96, '+')

# x and f are the values of dfs_corr and propMax from "dt_spX"
data_fgy<-read.table(file="3fenologia_plot.csv", header=T, sep=";", check.names=T, dec=".")

plot.new()

plot04<-ggplot(data_fgy, aes(x=x, y=f, color=taxa, shape=taxa)) +
  ylab("Proportion of \n flowering plants") + 
  xlab("Days") +
  xlim(-10,165)+ ylim(-0.1,0.6)+
  geom_vline(xintercept = c(10,41,72,102,133,163))+
  geom_hline(yintercept=0)+
  annotate("text", x=c(14,45,76,106,137), y=0.56, 
           label= c("Jul","Aug","Sep","Oct","Nov"),size=9,hjust=0)+
  theme_gray(base_size=30)

#with confidence intervals
plot05<-plot04 + geom_line(aes(x=x,y=pred,colour=taxa),size=1.5) +
  geom_ribbon(aes(ymax=upr.conf, ymin=lwr.conf, fill=taxa), 
              alpha = 0.3,show.legend =F )+
  theme(axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22,margin=unit(c(0,4,4,0),"mm")),
        axis.text=element_text(size=22), 
        panel.border=element_rect(colour = "black", fill=NA, size=3),
        legend.position="right", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.title = element_blank(),
        legend.text=element_text(size=18))+
  geom_point(aes(fill=taxa),colour="black",shape=21, size=4)+
  guides(col=FALSE)+
  scale_colour_manual(values = c("#E7AB00BF","#189E76BF", "#736DB3BF","#DA5D00BF"),
                      labels = element_blank()) +
  scale_fill_manual(values = c("#E7AB00BF","#189E76BF","#736DB3BF","#DA5D00BF"),
                    labels=c("chilensis-albidiflora", "chilensis", "litoralis","mutabilis"))
plot05

##
##POLLINATORS##
##

#data_including_functionalgroups
SLPat<-read.table(file="4pollinatorsGF.csv", header=T, sep=";", check.names=T, dec=".") # para los datos incluir el row.names (no para el factor de ordenamiento)

modelo_polinizadoresGF<-glmer(visitas~taxa+gf+(1|ID.1), family=poisson(), data=SLPat,
                  control=glmerControl(optimizer="bobyqa"),nAGQ=3 ) #To avoid a warning of nonconvergence, we specify a different optimizer with the argument control=glmerControl(optimizer="bobyqa"). Although the model will produce nearly identical results without the new argument, we prefer to use models without such warnings.https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
summary(modelo_polinizadoresGF)
drop1(modelo_polinizadoresGF, test="Chi")

# Model_pollinatordGF_Plot
myData<- aggregate(SLPat$tv,by = list(taxa = SLPat$taxa, gf = SLPat$gf),
                   FUN = function(x) c(mean = mean(x), sd = sd(x),
                                       n = length(x))) #genera la estad?stica descriptiva que va al gr?fico
myData <- do.call(data.frame, myData)
myData$se <- myData$x.sd / sqrt(myData$x.n) # error estandard
colnames(myData) <- c("taxa", "gf", "mean", "sd", "n", "se")

dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = mean + se,
              ymin = mean - se)

plot06 <- ggplot(data = myData, aes(x = gf, y = mean, fill = taxa))+
  geom_hline(yintercept = c(0,1,2))+labs(y=expression(Visits%.%hour^{-1}),size=18)+
  theme_gray(base_size=30)
plot07<-plot06 + geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.4,lwd=1) +
  theme(axis.text.x=element_text(size=15), 
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=20,margin=unit(c(0,2,0,0),"mm")),
        panel.border=element_rect(colour = "black", fill=NA, size=3.5),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.title = element_blank(),
        legend.text=element_text(size=22),
        plot.background = element_rect(fill = "gray93"))+
  scale_fill_manual(values = c("#E7AB00BF","#189E76BF","#736DB3BF","#DA5D00BF"),
                    labels=c("chilensis-albidiflora", "chilensis", "litoralis","mutabilis"))+
  scale_x_discrete(limits=c("small","large","patagona"), labels=c("Small bees","Large bees","Patagona gigas"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 2.5))

plot07  

#Stacked pollinators plot
df_stack_pol<-data.frame(
  pol=rep(c("Lip (s)","Ant (s)","Humm",	"Dli1 (s)",	"Chi1 (s)",	"Tri (l)",	"Chi2 (s)",
            "Other bees"),4),
  taxa=rep(c("alb","chi","lit","mut"),each=8),
  visits=c(0.163346614,	0.075697211,	0,	0.003984064,	0.001992032,	0.001992032,	0,	0.009960159,
           0.129482072,	0.039840637,	0,	0	,0.023904382,	0.007968127,	0,	0.009960159,
           0,	0,	0.103585657,	0.015936255,	0,	0,	0,	0.007968127,
           0.059760956,	0.185258964,	0,	0.079681275,	0.025896414,	0.015936255,	0.015936255,	0.021912351))

plot08 <- ggplot(df_stack_pol, aes(x = reorder(pol,visits), y = visits))+
  geom_col(aes(fill =taxa), width = 0.9)+
  ylab("Proportion of visits")+
  theme_gray(base_size=30)
plot09<-plot08 + theme(axis.title.x = element_text(size=24,margin=unit(c(4,7,7,0),"mm")),
               axis.title.y = element_blank(),
               axis.text=element_text(size=20), 
               panel.border=element_rect(colour = "black", fill=NA, size=3),
               legend.position="right", panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),legend.title = element_blank(),
               legend.text=element_text(size=18))+
  scale_fill_manual(values = c("#E7AB00BF","#189E76BF","#736DB3BF","#DA5D00BF"),
                    labels=c("chilensis-albidiflora", "chilensis", "litoralis","mutabilis"))+
  scale_x_discrete(labels=c("Chi2","Tri*","Chi1","Other bees","Dli1","Humm","Ant","Lip"))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 0.45))+
  coord_flip()

plot09


##
###MULTIDIMENSIONAL NICHE##
###

# Modified codes from Geange et al (2010) Methods in Ecology and Evolution https://doi.org/10.1111/j.2041-210X.2010.00070.x
##A unified analysis of niche overlap incorporating data of different types
# A Program to calculate niche overlaps over multiple niche axes per individual.

##First of all, you need to load functions used in this script
# Required package abind
# Read in the niche overlap functions:
source("Geange_et_al_2010_niche_functions.txt")

###   
# -----------------------------------------------------------
#          Analysis B: POLLINATORS
# -----------------------------------------------------------
### 

# Analysis for calculating niche overlap on a single niche axis,
# in this instance a categorical variable.
B.df <- read.table("5polinizadores_NMD.txt",T)
#list.files()
# Ensure the first two column names are "id" and "species".
colnames(B.df)[1] <- "id"
colnames(B.df)[2] <- "species"
# Ensure that the first 2 cols are factors.
B.df$id      <- as.factor(B.df$id)
B.df$species <- as.factor(B.df$species)
# Store some vectors of names:
spnames   <- sort(unique(as.character(B.df$species)))
no.spp    <- length(spnames)
varnames <- colnames(B.df)[-(1:2)]    
no.vars  <- length(varnames) 

# Make a vector of variable types to match the variable names:
vartypes <- c("cat")

# Check they are correctly labelled:
cbind(varnames,vartypes)
avail.list <- vector("list",no.vars) ###
names(avail.list) <- varnames ###

# Set up R objects to store results
# alpha.list
Balpha.list <- vector("list",no.vars)
names(Balpha.list) <- varnames

for (vv in 1:no.vars) if (vartypes[vv]=="rsel")
{
  choices <- unique(B.df[,vv+2])
  no.ch   <- length(choices)
  Balpha.list[[vv]] <- matrix(NA,no.spp,no.ch)
  dimnames(Balpha.list[[vv]]) <- list(spnames,choices)
}

# Set up an array of niche overlaps.
# The object no.array is an array of niche overlaps.
# It is a 3-D array, with rows and columns being species 
# (a square symmetric matrix for pairwise niche overlaps), 
# and the layers are the dimensions for the multivariate 
# niche overlap measure (one dimension per variable).
# Rows and columns are species, layers are variables.
Bno.array  <- array(1,c(no.spp,no.spp,no.vars))
dimnames(Bno.array) <- list(spnames,spnames,varnames)

# Run through each variable in turn, identify its type,
# calculate the appropriate NO matrix and store it in
# the right layer of the no.array.
for (vv in 1:no.vars)
{
  y <- B.df[,colnames(B.df)==varnames[vv]]
  if (vartypes[vv] == "bin")
    Bno.array[,,vv] <- no.bin.fn(B.df$species,y)
  if (vartypes[vv] == "cat")
    Bno.array[,,vv] <- no.cat.fn(B.df$species,y)
  if (vartypes[vv] == "count")
    Bno.array[,,vv] <- no.count.fn(B.df$species,y)
  if (vartypes[vv] == "cts")
    Bno.array[,,vv] <- no.cts.fn(B.df$species,y)
  if (vartypes[vv] == "meas")
    Bno.array[,,vv] <- no.cts.fn(B.df$species,log(y))
  if (vartypes[vv] == "pcent")
    Bno.array[,,vv] <- no.cts.fn(B.df$species,
                                 log(y/(100-y)))
  if (vartypes[vv] == "propn")
    Bno.array[,,vv] <- no.cts.fn(B.df$species,
                                 log(y/(1-y)))
  if (vartypes[vv] == "rsel")
  {
    
    # Do Manly's alpha calculations, store.
    avail.vect <- avail.list[[vv]]
    alpha.mat <- alpha.fn(B.df$species,y,avail.vect)
    alpha.list[[vv]] <- alpha.mat         
    
    # Do niche overlaps, as proportions in categories:
    Bno.array[,,vv] <- no.rsel.cat.fn(alpha.mat)
  }
}

# Analysis B  -  Permutation testing.
# -----------------------------------
# Permutation of the species niches.

replic<-999
Bpseudo.no.array  <- array(1,c(no.spp,no.spp,no.vars,replic)) #remplace replic por el n?mero
dimnames(Bpseudo.no.array) <- list(spnames,spnames,varnames,NULL)

# Set a temporary data frame, which will change each time
# through the cycle by having its species column permuted.
temp.df <- B.df

# For each replication, permute the species labels, run the
# niche overlap calculations, and store the results in the
# pseudo NO array.
for (rr in 1:replic)
{
  # Permute the species labels in the temporary dataframe:
  temp.df$species <- sample(temp.df$species)
  for (vv in 1:no.vars)
  {
    
    # Read out the column from this variable:
    y <- temp.df[,colnames(B.df)==varnames[vv]]
    
    # Run through the variable types, do appropriate analyses:
    if (vartypes[vv] == "bin")
      Bpseudo.no.array[,,vv,rr] <- no.bin.fn(temp.df$species,y)
    if (vartypes[vv] == "cat")
      Bpseudo.no.array[,,vv,rr] <- no.cat.fn(temp.df$species,y)
    if (vartypes[vv] == "count")
      Bpseudo.no.array[,,vv,rr] <- no.count.fn(temp.df$species,y)
    if (vartypes[vv] == "cts")
      Bpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,y)
    if (vartypes[vv] == "meas")
      Bpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,log(y))
    if (vartypes[vv] == "pcent")
      Bpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                             log(y/(100-y)))
    if (vartypes[vv] == "propn")
      Bpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                             log(y/(1-y)))
    if (vartypes[vv] == "rsel")
    {
      
      # Do Manly's alpha calculations, store.
      avail.vect <- avail.list[[vv]]
      alpha.mat  <- alpha.fn(temp.df$species,y,avail.vect)
      
      # Do niche overlaps, as proportions in categories:
      Bpseudo.no.array[,,vv,rr] <- no.rsel.cat.fn(alpha.mat)
    }
  }
  print(paste("Rep",rr,"done"))
}

Bpseudo.no.array
apply(Bpseudo.no.array,c(1,2),function(x) sd(as.vector(x))) #SD estimation from the array

###     
# -----------------------------------------------------------
#          Analysis C: TRAITS (MORPHOSPACE)
# -----------------------------------------------------------
### 

# Analysis for calculating niche overlap on a single niche axis,
# in this instance a measurement variable.

# The input data set is a .txt file.
# It needs to have its first two columns
# labelled "id" and "species". Subsequent columns are
# individual-based variables, one per column.

# Input the individual data file:

# ??????????????????????????????????????????????????????????
C.df <- read.table("6firstNMDSaxis.txt",T)
# ??????????????????????????????????????????????????????????

# Ensure the first two column names are "id" and "species".
colnames(C.df)[1] <- "id"
colnames(C.df)[2] <- "species"
# Ensure that the first 2 cols are factors.
C.df$id      <- as.factor(C.df$id)
C.df$species <- as.factor(C.df$species)
# Store some vectors of names:
spnames   <- sort(unique(as.character(C.df$species)))
no.spp    <- length(spnames)
varnames <- colnames(C.df)[-(1:2)]    
no.vars  <- length(varnames) 

# Make a vector of variable types to match the variable names:
vartypes <- c("cts")
# Check they are correctly labelled:
cbind(varnames,vartypes)

# Set up a list of objects which are NULL if this is not
# a resource selection variable, and with the availability
# vector if it is resource selection.
avail.list <- vector("list",no.vars)
names(avail.list) <- varnames

# Set up R objects to store results
# alpha.list
# The object alpha.list has one component per variable.
# The components are NULL for ordinary variables.
Calpha.list <- vector("list",no.vars)
names(Calpha.list) <- varnames
for (vv in 1:no.vars) if (vartypes[vv]=="rsel")
{
  choices <- unique(C.df[,vv+2])
  no.ch   <- length(choices)
  Calpha.list[[vv]] <- matrix(NA,no.spp,no.ch)
  dimnames(Calpha.list[[vv]]) <- list(spnames,choices)
}

# no.array
# Set up an array of niche overlaps.
# The object no.array is an array of niche overlaps.
# It is a 3-D array, with rows and columns being species 
# (a square symmetric matrix for pairwise niche overlaps), 
# and the layers are the dimensions for the multivariate 
# niche overlap measure (one dimension per variable).
# Rows and columns are species, layers are variables.
Cno.array  <- array(1,c(no.spp,no.spp,no.vars))
dimnames(Cno.array) <- list(spnames,spnames,varnames)

# Run through each variable in turn, identify its type,
# calculate the appropriate NO matrix and store it in
# the right layer of the no.array.
for (vv in 1:no.vars)
{
  y <- C.df[,colnames(C.df)==varnames[vv]]
  if (vartypes[vv] == "bin")
    Cno.array[,,vv] <- no.bin.fn(C.df$species,y)
  if (vartypes[vv] == "cat")
    Cno.array[,,vv] <- no.cat.fn(C.df$species,y)
  if (vartypes[vv] == "count")
    Cno.array[,,vv] <- no.count.fn(C.df$species,y)
  if (vartypes[vv] == "cts")
    Cno.array[,,vv] <- no.cts.fn(C.df$species,y)
  if (vartypes[vv] == "meas")
    Cno.array[,,vv] <- no.cts.fn(C.df$species,log(y))
  if (vartypes[vv] == "pcent")
    Cno.array[,,vv] <- no.cts.fn(C.df$species,
                                 log(y/(100-y)))
  if (vartypes[vv] == "propn")
    Cno.array[,,vv] <- no.cts.fn(C.df$species,
                                 log(y/(1-y)))
  if (vartypes[vv] == "rsel")
  {
    
    # Do Manly's alpha calculations, store.
    Cavail.vect <- avail.list[[vv]]
    Calpha.mat <- alpha.fn(C.df$species,y,Cavail.vect)
    alpha.list[[vv]] <- Calpha.mat         
    
    # Do niche overlaps, as proportions in categories:
    Cno.array[,,vv] <- no.rsel.cat.fn(Calpha.mat)
  }
}

# Analysis C  -  Permutation testing.
# -----------------------------------
# Permutation of the species labels would give data 
# satisfying the null model of complete niche overlap, 
# i.e. that none of the variables 
# serves to differentiate species into different niches.

# Hence for each replication, permute the species labels
# and run through all the calculations above.
# Stor NOs in an array with one extra dimension, one
# layer for each replication.
# Then the null distributions are all stored.
# Can use the original availability data, but need a new 
# alpha list each time.

# Set up array to store pseudo niche overlaps:
Cpseudo.no.array  <- array(1,c(no.spp,no.spp,no.vars,replic))
dimnames(Cpseudo.no.array) <- list(spnames,spnames,varnames,NULL)

# Set a temporary data frame, which will change each time
# through the cycle by having its species column permuted.
temp.df <- C.df

# For each replication, permute the species labels, run the
# niche overlap calculations, and store the results in the
# pseudo NO array.
for (rr in 1:replic)
{
  
  # Permute the species labels in the temporary dataframe:
  temp.df$species <- sample(temp.df$species)
  for (vv in 1:no.vars)
  {
    
    # Read out the column from this variable:
    y <- temp.df[,colnames(C.df)==varnames[vv]]
    
    # Run through the variable types, do appropriate analyses:
    if (vartypes[vv] == "bin")
      Cpseudo.no.array[,,vv,rr] <- no.bin.fn(temp.df$species,y)
    if (vartypes[vv] == "cat")
      Cpseudo.no.array[,,vv,rr] <- no.cat.fn(temp.df$species,y)
    if (vartypes[vv] == "count")
      Cpseudo.no.array[,,vv,rr] <- no.count.fn(temp.df$species,y)
    if (vartypes[vv] == "cts")
      Cpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,y)
    if (vartypes[vv] == "meas")
      Cpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,log(y))
    if (vartypes[vv] == "pcent")
      Cpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                             log(y/(100-y)))
    if (vartypes[vv] == "propn")
      Cpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                             log(y/(1-y)))
    if (vartypes[vv] == "rsel")
    {
      
      # Do Manly's alpha calculations, store.
      Cavail.vect <- avail.list[[vv]]
      Calpha.mat  <- alpha.fn(temp.df$species,y,Cavail.vect)
      
      # Do niche overlaps, as proportions in categories:
      Cpseudo.no.array[,,vv,rr] <- no.rsel.cat.fn(Calpha.mat)
    }
  }
  print(paste("Rep",rr,"done"))
}

Cno.array
apply(Cpseudo.no.array,c(1,2),function(x) sd(as.vector(x))) # SD estimation from array

##-----------------------------------------------------------
#          Analysis A: PHENOLOGY
# -----------------------------------------------------------
### 

# Analysis for calculating niche overlap on a single niche axis,
# in this instance a continuous variable.

# The input data set is a .txt file.
# It needs to have its first two columns
# labelled "id" and "species". Subsequent columns are
# individual-based variables, one per column.

# Input the individual data file:

# ?????????????????????????????????????????????????????????? 
# read in the individual data file
A.df <- read.table("7phenologia_MDN.txt",T)
# ??????????????????????????????????????????????????????????

# Ensure the first two column names are "id" and "species".
colnames(A.df)[1] <- "id"
colnames(A.df)[2] <- "species"
# Ensure that the first 2 cols are factors.
A.df$id      <- as.factor(A.df$id)
A.df$species <- as.factor(A.df$species)
# Store some vectors of names:
spnames   <- sort(unique(as.character(A.df$species)))
no.spp    <- length(spnames)
varnames <- colnames(A.df)[-(1:2)]    
no.vars  <- length(varnames) 

# Make a vector of variable types to match the variable names:
vartypes <- c("cts")
# Check they are correctly labelled:
cbind(varnames,vartypes) #das mean "day after solstice"

# Set up a list of objects which are NULL if this is not
# a resource selection variable, and with the availability
# vector if it is resource selection.
avail.list <- vector("list",no.vars)
names(avail.list) <- varnames

# Set up R objects to store results
# ---------------------------------
# alpha.list
# The object alpha.list has one component per variable.
# The components are NULL for ordinary variables.
Aalpha.list <- vector("list",no.vars)
names(Aalpha.list) <- varnames

for (vv in 1:no.vars) if (vartypes[vv]=="rsel")
{
  choices <- unique(A.df[,vv+2])
  no.ch   <- length(choices)
  Aalpha.list[[vv]] <- matrix(NA,no.spp,no.ch)
  dimnames(Aalpha.list[[vv]]) <- list(spnames,choices)
}

# no.array
# Set up an array of niche overlaps.
# The object no.array is an array of niche overlaps.
# It is a 3-D array, with rows and columns being species 
# (a square symmetric matrix for pairwise niche overlaps), 
# and the layers are the dimensions for the multivariate 
# niche overlap measure (one dimension per variable).
# Rows and columns are species, layers are variables.
Ano.array  <- array(1,c(no.spp,no.spp,no.vars))
dimnames(Ano.array) <- list(spnames,spnames,varnames)

# Run through each variable in turn, identify its type,
# calculate the appropriate NO matrix and store it in
# the right layer of the no.array.
for (vv in 1:no.vars)
{
  y <- A.df[,colnames(A.df)==varnames[vv]]
  if (vartypes[vv] == "bin")
    Ano.array[,,vv] <- no.bin.fn(A.df$species,y)
  if (vartypes[vv] == "cat")
    Ano.array[,,vv] <- no.cat.fn(A.df$species,y)
  if (vartypes[vv] == "count")
    Ano.array[,,vv] <- no.count.fn(A.df$species,y)
  if (vartypes[vv] == "cts")
    Ano.array[,,vv] <- no.cts.fn(A.df$species,y)
  if (vartypes[vv] == "meas")
    Ano.array[,,vv] <- no.cts.fn(A.df$species,log(y))
  if (vartypes[vv] == "pcent")
    Ano.array[,,vv] <- no.cts.fn(A.df$species,
                                 log(y/(100-y)))
  if (vartypes[vv] == "propn")
    Ano.array[,,vv] <- no.cts.fn(A.df$species,
                                 log(y/(1-y)))
  if (vartypes[vv] == "rsel")
  {
    
    # Do Manly's alpha calculations, store.
    avail.vect <- avail.list[[vv]]
    alpha.mat <- alpha.fn(A.df$species,y,avail.vect)
    alpha.list[[vv]] <- alpha.mat         
    
    # Do niche overlaps, as proportions in categories:
    Ano.array[,,vv] <- no.rsel.cat.fn(alpha.mat)
  }
}

# Analysis A  -  Permutation testing.
# -----------------------------------
# Permutation of the species niches.

# Set up array to store pseudo niche overlaps:
Apseudo.no.array  <- array(1,c(no.spp,no.spp,no.vars,replic))
dimnames(Apseudo.no.array) <- list(spnames,spnames,varnames,NULL)

# Set a temporary data frame, which will change each time
# through the cycle by having its species column permuted.
temp.df <- A.df

# For each replication, permute the species labels, run the
# niche overlap calculations, and store the results in the
# pseudo NO array
for (rr in 1:replic)
{
  
  # Permute the species labels in the temporary dataframe:
  temp.df$species <- sample(temp.df$species)
  for (vv in 1:no.vars)
  {
    
    # Read out the column from this variable:
    y <- temp.df[,colnames(A.df)==varnames[vv]]
    
    # Run through the variable types, do appropriate analyses:
    if (vartypes[vv] == "bin")
      Apseudo.no.array[,,vv,rr] <- no.bin.fn(temp.df$species,y)
    if (vartypes[vv] == "cat")
      Apseudo.no.array[,,vv,rr] <- no.cat.fn(temp.df$species,y)
    if (vartypes[vv] == "count")
      Apseudo.no.array[,,vv,rr] <- no.count.fn(temp.df$species,y)
    if (vartypes[vv] == "cts")
      Apseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,y)
    if (vartypes[vv] == "meas")
      Apseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,log(y))
    if (vartypes[vv] == "pcent")
      Apseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                             log(y/(100-y)))
    if (vartypes[vv] == "propn")
      Apseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                             log(y/(1-y)))
    if (vartypes[vv] == "rsel")
    {
      
      # Do Manly's alpha calculations, store.
      avail.vect <- avail.list[[vv]]
      alpha.mat  <- alpha.fn(temp.df$species,y,avail.vect)
      
      # Do niche overlaps, as proportions in categories:
      Apseudo.no.array[,,vv,rr] <- no.rsel.cat.fn(alpha.mat)
    }
  }
  print(paste("Rep",rr,"done"))
}

Ano.array
apply(Apseudo.no.array,c(1,2),function(x) sd(as.vector(x))) #SD estimation from the array

###   
# -----------------------------------------------------------
# Combine the three axes into a single dataframe and 
#  calculate mean niche overlap across axes
# en nuestro caso RASGOS, POLINIZADORES Y FENOLOG?A
# ----------------------------------------------------------- 
### 

# Combine the three axes (habitat, position and pfar) into a single dataframe   
#A. Phenology==> type="cts"
#B. Pollinators==> type="cat"
#C. Traits==> type="cts"
# Follow the same order used in individual axis characterization
no.all.mat <- abind(Bno.array,Cno.array,Ano.array,along=3) 

# calculate mean niche overlap across axes  
mean_overlap <- apply(no.all.mat,1:2,mean)
# calculate the associated standard deviation
sd_overlap <- apply(no.all.mat,1:2,sd)

# Combine the three pseudo datasets (habitat, position and pfar) into a single dataframe   
pseudo.no.all.mat <- abind(Bpseudo.no.array,Cpseudo.no.array,Apseudo.no.array,along=3) 

# For each replicate, calculate mean niche overlap across axes 
pseudo.mean_overlap <- apply(pseudo.no.all.mat,c(1:2,4),mean)
pseudo.mean_overlap

###
#--------------------------------------------------------
# Null model analysis determining if the niches
# of two species in niche space differ
#--------------------------------------------------------
###
# Calculate p values for each pair of species 
# separately for each variable.

# calculate number of axes
no.axes <- dim(no.all.mat)[3]

# assign name to each axis
axis.names <- c("Pollinators","Traits", "Phenology") 

# calculate seperate p-values for each axis
sep.pvals     <- array(1,c(no.spp,no.spp,no.axes))
dimnames(sep.pvals) <- list(spnames,spnames,axis.names)

for (spa in 1:(no.spp-1)) for (spb in (spa+1):no.spp)
  for (vv in 1: no.axes)   
  {
    pseudo.nos <- pseudo.no.all.mat[spa,spb,vv,]
    data.no    <- no.all.mat[spa,spb,vv]
    sep.pvals[spa,spb,vv] <- mean(pseudo.nos<data.no) 
    length(pseudo.nos[data.no<pseudo.nos])
    sep.pvals[spb,spa,vv] <- sep.pvals[spa,spb,vv] 
  }

# Also calculate a p-value for overall NO measure averaged across axes
overall.pvals <- matrix(1,no.spp,no.spp)
overall.pvals
dimnames(overall.pvals) <- list(spnames,spnames)

for (spa in 1:(no.spp-1)) for (spb in (spa+1):no.spp)
{
  temp.mat  <- pseudo.no.all.mat[spa,spb,,]
  pseudo.nos <- apply(temp.mat,2,mean)
  data.no    <- mean_overlap[spa,spb]
  overall.pvals[spa,spb] <- mean(pseudo.nos<data.no) 
  length(pseudo.nos[data.no<pseudo.nos])
  overall.pvals[spb,spa] <- overall.pvals[spa,spb] 
}

#--------------------------------------------------------
# Null model analysis determining if the distribution of
# species across niche space are more differentiated
# or more clustered than expected
#--------------------------------------------------------

# First, reformat the observed data to derive a matrix of niche overlaps
# with one row per species combination, and one column for each niche dimension

VV <- no.axes  # Number of axes
RR <- replic   # Number of replications.

no.mat <- matrix(NA,(no.spp*(no.spp-1)/2),VV)
for (vv in 1:VV)
  no.mat[,vv] <- as.vector(as.dist(no.all.mat[,,vv]))

# Next, reformat the pseudo data to derive a matrix of niche overlaps
# with one row per species, and one column for each niche dimension,
# with one extra dimension, one layer for each replication

pseudo.mat <- 	array(NA,c((no.spp*(no.spp-1)/2),VV,replic))
for (vv in 1:VV) for (rr in 1:RR)
  pseudo.mat[,vv,rr] <- as.vector(as.dist(pseudo.no.all.mat[,,vv,rr]))

# --------------------------------------------------------
# For each niche dimension, calculate mean and variance over the species
# pairs, and hence the test statistic ch = coefficient of heterogeneity.
# Note: Need to use variance formula based on n, not n-1.

KK <- ncol(no.mat)      # Number of niche dimensions
SS <- nrow(no.mat)      # Number of species pairs
RR <- replic            # Number of replications.

data.ch <- rep(NA,KK)
pseudo.ch <- matrix(NA,RR,KK)

for (kk in 1:KK)
{
  # Calculate data test statistic:
  x <- mean(no.mat[,kk])
  v <- var(no.mat[,kk])*(SS-1)/SS # Adjust for denom n, not n-1
  data.ch[kk] <- v/x/(1-x)
  
  # Calculate test stats for all pseudo-data:
  for (rr in 1:RR)
  {
    x <- mean(pseudo.mat[,kk,rr])
    v <- var(pseudo.mat[,kk,rr])*(SS-1)/SS
    pseudo.ch[rr,kk] <- v/x/(1-x)
  }
}

# For each niche dimension, see if data more differentiated than random.
p.dims.diff <- rep(NA,KK)
for (kk in 1:KK)
  p.dims.diff[kk] <- mean(data.ch[kk] > pseudo.ch[,kk])
names(p.dims.diff) <- paste("diff.dim",sort(axis.names))

# For each niche dimension, see if data more clustered than random.
p.dims.clus <- rep(NA,KK)
for (kk in 1:KK)
  p.dims.clus[kk] <- mean(data.ch[kk] < pseudo.ch[,kk])
names(p.dims.clus) <- paste("clus.dim",sort(axis.names))

# --------------------------------------------------------
# For average niche overlap, calculate mean and variance over the species
# pairs, and hence the test statistic ch = coefficient of heterogeneity.
# Note: Need to use variance formula based on n, not n-1.

overall.data.ch   <- mean(data.ch)
overall.pseudo.ch <- apply(pseudo.ch,1,mean)

# Test if this community is more differentiated than random:
p.all.diff <- mean(overall.data.ch > overall.pseudo.ch)

# Test if this community is more clustered than random:
p.all.clus <- mean(overall.data.ch < overall.pseudo.ch)

###
#--------------------------------------------------------
# Save all results of the analysis:
#--------------------------------------------------------
###

NOb.results <- list(
  info = list(variables = cbind(axis.names,c("cat","cts","cts")), #"resl","cat","meas"
              perm.reps = replic),
  NOestimates = no.all.mat,
  separate.pvalues = sep.pvals,
  separate.cluster.pvalues = p.dims.clus,
  separate.differentiated.pvalues = p.dims.diff,
  ests.overall = mean_overlap,
  ests.overall.sd = sd_overlap,
  overall.pvalues = overall.pvals,
  overall.cluster.pvalues = p.all.clus,
  overall.differentiated.pvalues = p.all.diff)

# To inspect results later, type in
#    names(NOb.results)
# to decide what to look at. Then type (e.g.)
#    NOb.results$NOestimates
# to see that component of the list.

# Save the NOb.results with a more informative names for 
# your own data set.
Eriosyce.results<-NOb.results
Eriosyce.results

###
###REPRODUCTIVE SUCCESSS##
###

RS<-read.table("8Eriosyce_RS.csv", header=T, sep=";",dec=".")

# FRUIT_SET_MODEL
RS$Species<-factor(RS$Species) 
is.factor(RS$Species)
str(RS)

model_FS<-glm(fruit.set~Species, data=RS,family=binomial())
drop1(model_FS,test="Chisq")

#PLOT_FRUIT_set
mu <- ddply(RS, "Species", summarise, grp.mean=mean(fruit.set))
head(mu)

plot.new()
plot11<-ggplot(RS, aes(x=fruit.set,y=Species,fill=Species)) + 
  geom_density_ridges(alpha=0.8,scale=1,jittered_points=T,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1) + 
  xlab('Fruit set (proportion)')+
  scale_fill_manual(values = c("#E7AB00BF","#189E76BF","#736DB3BF","#DA5D00BF"))+
  geom_vline(xintercept = c(0,1))+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),
        axis.text=element_text(size=22),
        panel.border=element_rect(colour = "black", fill=NA, size=5),
        legend.position="none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.title = element_blank())+
  theme_gray(base_size=30)
plot11

# Number_of_seeds_MODEL
model_SN<-glm(Number.of.seeds~Species, 
              data=RS,family=poisson())
drop1(model_SN,test="Chisq")
comparacion = emmeans(model_SN, ~ Species)
pairs(comparacion)

plot12<-ggplot(RS, aes(x=Number.of.seeds,y=Species,fill=Species)) + 
  geom_density_ridges(alpha=0.8,scale=1,jittered_points=T,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1) + 
  xlab('Number of seeds')+
  scale_fill_manual(values = c("#E7AB00BF","#189E76BF","#736DB3BF","#DA5D00BF"))+
  geom_vline(xintercept = c(0,1))+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),
        axis.text=element_text(size=22),
        panel.border=element_rect(colour = "black", fill=NA, size=5),
        legend.position="none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.title = element_blank())+
  theme_gray(base_size=30)
plot12

#GERMINATION
germ<-read.table("9Eriosyce_germination.csv", header=T, sep=";",dec=".")
model_germ<-glm(germination~species, data=germ, 
                family=quasibinomial)
drop1(model_germ,test="Chisq")

comparacion = emmeans(model_germ, ~ species)
pairs(comparacion)

# Germination_plot
detach(package:plyr)
df.summary <- germ %>%
  group_by(species) %>%
  summarise(
    sd = sd(germination, na.rm = TRUE),
    germination = mean(germination),
  )
df.summary

ggplot(df.summary, aes(species, germination)) +
  geom_pointrange(
    aes(ymin = germination-sd, ymax = germination+sd,color=species))+
  geom_point(aes(color=species),
             position=position_dodge(0.7),size=4.5)+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1.1))+
  scale_color_manual(values = c("#E7AB00BF","#189E76BF","#736DB3BF","#DA5D00BF", "black"))+
  ylab("Germination (proportion)")+
  theme_gray(base_size=30)+theme(legend.position="none")

##
###POLLEN LIMITATION##
##

pol_lim<-read.table("11Eriosyce_PL.csv", header=T, sep=";",dec=".")
pl_alb<-pol_lim[c(17:27),c(7:8)]
pl_chi<-pol_lim[c(1:16),c(7:8)]
pl_lit<-pol_lim[c(48:67),c(7:8)]
pl_mut<-pol_lim[c(28:47),c(7:8)]

#PL, based on Larson & Barrett (2000)
fc_seeds <- function(data,i){
  d <- data[i,]
  Po_s<-sum(na.omit(d[,1]))/length(na.omit(which(d[,1]>0)))
  Pc <-sum(na.omit(d[,2]))/length(na.omit(which(d[,2]>0)))        
  return(1-(Po_s/Pc))
} #PL in seeds

#Indice de LArson y Barret (2000)
# PL = 1 - (Po/Pc)
#Po fruit set/seed number in the open flowers
#Pc fruit set7seed number in the closed
fc_seeds(pl_alb) #

set.seed(626)
#library(boot)
bootcorr<-boot(pl_alb, fc_seeds, R=999) # 
boot.ci(bootcorr,conf=0.95,type="norm") 

# Test the distribution (PL)
t.test(bootcorr$t[,1],alternative="greater", mu=0, paired=F,conf.level=0.95)
#Descriptive statistics
pol_lim2<-pol_lim[,c(9:11)]
df.pl<-pol_lim2 %>%
  group_by(taxa,treatemt) %>%
  summarise(
    sd = sd(seeds, na.rm = T),
    N_seeds = mean(seeds)
  )
df.pl

ggplot(aes(y = seeds, x = treatemt,fill=treatemt), data = pol_lim2) + geom_boxplot()+
  facet_grid(.~ taxa)
