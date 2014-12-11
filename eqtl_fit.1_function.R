#fit1.2 function

# this is a function to determine the best model for model selection 
#data.model can be either "normal" or "binary
#perms in needs to be a list of all permutations with the names "nocovar", addcovar, intcovar, "_" model type
fit.1<-function(cross=cross, pheno.in,covar.in, data.model="normal", nocovar.perms, addcovar.perms, intcovar.perms){
  model.out<-list()
  #run 3 separate models:1)no covariate, 2)additive covariate, 3)additive and interactive covariate
  s1.no<-scanone(cross, pheno.col=pheno.in, method="hk", addcovar=NULL, intcovar=NULL, penalties.in, model=data.model)
  s1.add<-scanone(cross, pheno.col=pheno.in, method="hk", addcovar=cov, intcovar=NULL, model=data.model)
  s1.int<-scanone(cross, pheno.col=pheno.in, method="hk", addcovar=covar.in, intcovar=covar.in, model=data.model)
  #get permutations for each model type
  perm.no<-subset(perms.in[[paste("nocovar",data.model, sep="_"]]
  perm.add<-perms.in[[paste("addcovar",data.model, sep="_"]]
  perm.int<-perms.in[[paste("intcovar",data.model, sep="_"]]
  #return pvalues for each model
  p.no<-min(summary(s1.no, perms=perm.no, pvalues=T)$pval)
  p.add<-min(summary(s1.add, perms=perm.add, pvalues=T)$pval)
  p.int<-min(summary(s1.int, perms=perm.int, pvalues=T)$pval)
  #which is the minimum
  minp<-which(c(p.no,p.add,p.int)==min(c(p.no,p.add,p.int)))
  #which is the best model?
  if(length(minp)==1){
    bestmod<-c("no","add","int")[minp]
  }else{
    l.no<-max(summary(s1.no, perms=perm.no, pvalues=T)$lod)
    l.add<-max(summary(s1.add, perms=perm.add, pvalues=T)$lod)
    l.int<-max(summary(s1.int, perms=perm.int, pvalues=T)$lod)
    maxl<-which(c(l.no,l.add,l.int)==max(c(l.no,l.add,l.int)))
    bestmod<-c("no","add","int")[maxl]
    minp<-maxl
  }
  #generate the best model depending on what was found above
  if(bestmod=="no"){
    chr.out<-max(s1.no)$chr;  pos.out<-max(s1.no)$pos
    mod.out<-makeqtl(cross, chr=chr.out,pos=pos.out, what="prob")
    form.out<-"y ~ q1";  addcovar.out<-NULL
  }else{
    if(bestmod=="add"){
      chr.out<-max(s1.add)$chr;  pos.out<-max(s1.add)$pos
      mod.out<-makeqtl(cross, chr=chr.out,pos=pos.out, what="prob")
      form.out<-"y ~ trt + q1";   addcovar.out<-covar.in
    }else{
      chr.out<-max(s1.int)$chr;  pos.out<-max(s1.int)$pos
      mod.out<-makeqtl(cross, chr=chr.out,pos=pos.out, what="prob")
      form.out<-"y ~ trt + q1 + trt*q1";  covar.out<-covar.in
    }
  }
  #fit the best model
  fit.out<-fitqtl(cross, qtl=mod.out, formula=form.out, pheno.col=pheno.in, covar=covar.out, method="hk", dropone=T, model=data.model)
  ref<- refineqtl(cross, qtl=mod.out, formula=form.out, pheno.col=pheno.in, covar=covar.out, method="hk", model=data.model)
  #return the 0.05 threshold corrected lod scores (here, "qtlplod")
  qtlplod.out<-max(summary(s1.no, perms=perm.no, pvalues=T)$lod)-summary(perm.no)[1]
  lod.out<-fit.out$result.full[1,4]
  plod.out<-calc.plod(lod.int, countqtlterms("y ~ q1", ignore.covar=T), penalties=penalties.in["int",])
  #generate convidence interval
  ciout<-bayesint(ref,qtl.index=1, expandtomarkers=F, prob=0.95)
  chr.ci<-ref$chr
  lowmarker<-rownames(ciout)[1]
  highmarker<-rownames(ciout)[3]
  lowposition<-ciout[1,2]
  highposition<-ciout[3,2]
  cis<-rbind(cis,cbind(chr.ci,lowmarker,highmarker,lowposition,highposition))
  #output model attributes
  model.out[["model"]]<-mod.out
  model.out[["formula"]]<-formula(fit.out)
  model.out[["plod"]]<-plod.out
  model.out[["fit"]]<-summary(fit.add)
  model.out[["cis"]]<-cis
  model.out[["ref"]]<-ref
  names(model.out)<-c(paste("model", pheno.in, sep="."),
                      paste("plod", pheno.in, sep="."),
                      paste("formula", pheno.in, sep="."),
                      paste("fit", pheno.in, sep="."),
                      paste("cis", pheno.in, sep="."),
                      paste("ref", pheno.in, sep="."))
  return(model.out)
}