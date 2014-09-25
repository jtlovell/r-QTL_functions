#QTL workflow:
##################################
#1 Get the environment set up:
rm(list=ls())
pkg <- c("qtl","snow","rlecuyer","plyr","ggplot2","reshape","reshape2","qvalue","lme4","nlme","car","doBy","scales")
invisible(lapply(pkg, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))
foo <- function() { if (!is.null(seed <- getOption("myseed"))); set.seed(seed)}
options(myseed = 42)
sessionInfo()

##################################
# Required data input
# cross: the cross object, with genotype probabilities calculated, phenotype normalized and genotyping errors solved
cross<-read.cross("csvs",dir="", genotypes=c("a","b"),
                  genfile="....csv", phefile= "....csv",
                  na.strings=c("NA","-", "#VALUE!","#NUM!"))
cross<-calc.genoprob(...)

# physical positions of markers on reference... may need to be culled by the markers in the cross if some are excluded
physpos.marker<-read.csv(...)
physpos.marker<-physpos.marker[physpos.marker$Marker %in% markernames(cross),]
markernames(cross)[!markernames(cross) %in% physpos.marker$Marker] #remove markers without positions
cross.markers<-as.data.frame(pull.map(cross,as.table=T))
cross.markers$Marker<-rownames(cross.markers)
marker.info<-merge(physpos.marker,cross.markers, by="Marker")[,-4]

#physical position of all genes
cmgenpos<-read.csv("....csv",header=T, na.strings="NA") 
physpos.gene<-cmgenpos[,c(1,4:5)]

#gene expression data from experiment
geneexp<-read.csv("....csv",header=T) 
geneexp<-merge(idnames,geneexp, by="line")

#eqtls if you have them...
eqtls.wet<-read.csv("....csv", header=T)
eqtls.dry<-read.csv("....csv", header=T)

#marker ids
idnames<-read.csv("....csv",header=T)

#covariate data
cyto<-read.csv("....csv", header=T)
cyto.covar<-as.data.frame(...)

colnames(idnames)[2]<-"line"

#if needed, to quantile normalize phenotype data
quantnorm<-function(x) {  n=sum(!is.na(x),na.rm=T);  x=rank(x)/(n+1);  x=qnorm(x);  x[is.infinite(x)]=NA;  x}
cross.qn<-transformPheno(cross, pheno.col="FT.June", transf=quantnorm)

allphes<-phenames(cross)[phenames(cross)!="id"] 

#read in the cross
cross<-read.cross("csvs",dir="", genotypes=c("a","b"),
                  genfile="TKrils_map55.csv", phefile= "TKrils_full_phe_final_nogo.csv",
                  na.strings=c("NA","-", "#VALUE!","#NUM!"))
cross<-convert2riself(cross)
cross<-calc.genoprob(cross, step=1, stepwidth="max",error.prob=0.01, map.function="kosambi")


##################################
#2 Run scanone perms- figure out which qtl phenotypes to retain. 
cross.qn<-transformPheno(cross, pheno.col=2:nphe(cross), transf=quantnorm)
# run permutations for a single phenotype- since it is qned, all other phenotypes should follow
s1.perms<-scanone(cross, pheno.col=allphes, method="hk",n.perm=1000, n.cluster=22, addcovar=NULL)
s1.perms.qn<-scanone(cross.qn, pheno.col=allphes, method="hk",n.perm=1000, n.cluster=22, addcovar=NULL)

s1.perms.cov<-scanone(cross, pheno.col=allphes, method="hk",n.perm=1000, n.cluster=22, addcovar=cyto.covar)
s1.perms.qn.cov<-scanone(cross.qn, pheno.col=allphes, method="hk",n.perm=1000, n.cluster=22, addcovar=cyto.covar)
save(s1.perms, s1.perms.cov,s1.perms.qn, s1.perms.qn.cov, file="...")

##################################
#3 Determine the traits that have QTL and what model gives them to you...
s1<-scanone(cross, pheno.col=allphes, method="hk", addcovar=NULL)
s1.cov<-scanone(cross, pheno.col=allphes, method="hk",addcovar=cyto.covar)
s1.qn<-scanone(cross.qn, pheno.col=allphes, method="hk", addcovar=NULL)
s1.qn.cov<-scanone(cross.qn, pheno.col=allphes, method="hk",addcovar=cyto.covar)
# this function needs 2-8 inputs... the results from above
qtl.test<- s1.qtlcheck(s1=s1, s1.perms=s1.perms,
                       s1.cov=s1.cov, s1.perms.cov=s1.perms.cov,
                       s1.qn=s1.qn, s1.perms.qn=s1.perms.qn,
                       s1.qn.cov=s1.qn.cov, s1.perms.qn.cov=s1.perms.qn.cov,
                       alpha=0.1)
#check to see which phenotypes have QTL
qtlphes<-as.character(qtl.test$phenotype[qtl.test$s1.raw.qtl==TRUE | 
                                           qtl.test$s1.qn.qtl==TRUE |
                                           qtl.test$s1.qn.qtl==TRUE |
                                           qtl.test$s1.qn.qtl==TRUE ])

##################################
#4 Run the permutations w/ and w/o covariates
#this can be batched normally, but there is a bug in R/qtl that does not permit stratified permutations when batching
#so, run through a loop- slower, but gets the job done
perms.out<-list()
perms.out.cov<-list()
for (i in qtlphes){
  print(i)
  perms.out[[i]]<-scantwo(cross, pheno.col=i, method="hk", n.perm=..., n.cluster=.., addcovar=NULL)
  print("covar")
  perms.out.cov[[i]]<-scantwo(cross, pheno.col=i, method="hk",n.perm=..., n.cluster=.., addcovar=cyto.covar, perm.strata=cyto$cyto.num)
}
save(perms.out, perms.out.cov, file="...")

##################################
#5 Process scantwo and generate penalties
pens.out<-data.frame()
for (i in qtlphes){pens<-calc.penalties(perms.out[[i]], alpha=.1); pens.out<-rbind(pens.out,pens)}
pens.out.cov<-data.frame()
for (i in qtlphes){pens<-calc.penalties(perms.out.cov[[i]], alpha=.1); pens.out.cov<-rbind(pens.out.cov,pens)}
rownames(pens.out)<-qtlphes; colnames(pens.out)<-c("main","heavy","light")
rownames(pens.out.cov)<-qtlphes; colnames(pens.out.cov)<-c("main","heavy","light")

##################################
#6 Run stepwise model selection
#get the function, "stepwise.summary", which conducts stepwise model selection then outputs the statistics
if(!exists("stepwiseqtl.summary", mode="function")) source("stepwiseqtl.summary_function.R")
model.out<-stepwiseqtl.summary(cross=cross, phes=qtlphes, max.qtls=6, 
                               pens=pens.out, covar=NULL, printout=T, ci.method="bayes",prob=.9)
nocovar.modelstats<-read.csv("...",header=T)
write.csv(nocovar.modelstats, file="...")
model.out.cyto<-stepwiseqtl.summary(cross=cross, phes=qtlphes, max.qtls=6, 
                                    pens=pens.out.cov, covar=as.data.frame(cyto$cyto.num), printout=T, ci.method="bayes",prob=.9)
cytocovar.modelstats<-read.csv("...",header=T)
write.csv(cytocovar.modelstats, file="...")

save(model.out,nocovar.modelstats,model.out.cyto,cytocovar.modelstats,
     file="...")

##################################
#7 extract information from models, choose which models are best
# then extract statistics from these models and output
#four inputs are the two models  and two dataframes of statistics from summary.stepwise
#output can be "stats"- for a dataframe of best statistics, "models"- for a list of best models, or "types"- for a dataframe with pLOD scores and categories for each phenotype 
if(!exists("choose.best.model", mode="function")) source("choose.best.model.function.R")
bestmodels<-choose.best.model(model.out=model.out,  model.out.cyto=model.out.cyto,
                              nocovar.modelstats=nocovar.modelstats, cytocovar.modelstats=cytocovar.modelstats,   output="models")
bestmodel.types<-choose.best.model(model.out=model.out,  model.out.cyto=model.out.cyto,
                                   nocovar.modelstats=nocovar.modelstats, cytocovar.modelstats=cytocovar.modelstats, output="dataframe")
beststats<-choose.best.model(model.out=model.out,  model.out.cyto=model.out.cyto,
                             nocovar.modelstats=nocovar.modelstats, cytocovar.modelstats=cytocovar.modelstats, output="stats")
write.csv(beststats, file="...", row.names=F)
beststats<-read.csv("...",header=T)
qtl.ids<-paste(beststats$phenotype,"_",as.numeric(beststats$chromosome),"_",round(as.numeric(beststats$position)),sep="")
beststats$qtl.id<-qtl.ids

##################################
#8  extract all lps into a dataframe
#one input is necessary, a list of models, typically this comes from choose.best.model, but can also come from summary.stepwiseqtl
if(!exists("extract.lps", mode="function")) source("extract.lps.function.R")
lps.df<-extract.lps(bestmodels=bestmodels, beststats=beststats,qtlphes=qtlphes)
##################################
#9 make plots of positions and lod profiles
#this function makes a nice matrix plot of confidence intervals and standardized lod scores.
#it requires two inputs, the stats and models, typically from above, but from any summary.stepwise function would be fine
if(!exists("plot.lp.ci", mode="function")) source("plot.lp.ci.function.R")
plot.lp.ci(beststats=beststats,bestmodels=bestmodels)

##################################
#10 extract all genes within these intervals
#associate genes with markers
genes.by.markers.out<-gene.by.marker(physpos.gene=physpos.gene,marker.info=marker.info)
#based on eQTLs
if(!exists("cis.eQTL.select", mode="function")) source("cis.eQTL.select.function.R")
candidates.byeqtl<-cis.eQTL.select(cross=cross,
                                   beststats=beststats,
                                   eqtls.wet=eqtls.wet,
                                   eqtls.dry=eqtls.dry,
                                   method="marker.bound",
                                   genes.by.markers.out=genes.by.markers.out)
#other method not needing eQTLs are underway. 


##################################
#11 run covariate scan on all candidate genes
if(!exists("covariate.scan", mode="function")) source("covariate.scan.function.R")
covscan.out<-covariate.scan(beststats=beststats,
                            cross=cross,
                            bestmodel.types=bestmodel.types,
                            candidates.byeqtl=candidates.byeqtl,
                            qtl.ids=beststats$qtl.id[-c(5:6)])

