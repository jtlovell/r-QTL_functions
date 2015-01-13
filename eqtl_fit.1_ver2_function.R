# this is a function to determine the best model for model selection 

fit.1<-function(cross=cross, pheno.in,covar.in,penalties.in){
  model.out<-list()
  s1.no<-scanone(cross, pheno.col=pheno.in, method="hk", addcovar=NULL, intcovar=NULL)
  chr.no<-max(s1.no)$chr
  pos.no<-max(s1.no)$pos
  mod.no<-makeqtl(cross, chr=chr.no,pos=pos.no, what="prob")
  fit.no<-fitqtl(cross, qtl=mod.no, formula="y ~ q1", pheno.col=pheno.in, covar=NULL, method="hk", dropone=T)
  lod.no<-fit.no$result.full[1,4]
  plod.no<-calc.plod(lod.no, countqtlterms("y ~ q1", ignore.covar=F), penalties=penalties.in["no",])
  
  s1.add<-scanone(cross, pheno.col=pheno.in, method="hk", addcovar=covar.in, intcovar=NULL)
  chr.add<-max(s1.add)$chr
  pos.add<-max(s1.add)$pos
  mod.add<-makeqtl(cross, chr=chr.add,pos=pos.add, what="prob")
  fit.add<-fitqtl(cross, qtl=mod.add, formula="y ~ trt + q1", pheno.col=pheno.in, covar=covar.in, method="hk", dropone=T)
  lod.add<-fit.add$result.full[1,4]
  plod.add<-calc.plod(lod.add, countqtlterms("y ~ trt + q1", ignore.covar=F), penalties=penalties.in["add",])
  
  s1.int<-scanone(cross, pheno.col=pheno.in, method="hk", addcovar=covar.in, intcovar=covar.in)
  chr.int<-max(s1.int)$chr
  pos.int<-max(s1.int)$pos
  mod.int<-makeqtl(cross, chr=chr.int,pos=pos.int, what="prob")
  fit.int<-fitqtl(cross, qtl=mod.int, formula="y ~ trt + q1 + trt*q1", pheno.col=pheno.in, covar=covar.in, method="hk", dropone=T)
  lod.int<-fit.int$result.full[1,4]
  plod.int<-calc.plod(lod.int, countqtlterms("y ~ trt + q1 + trt*q1", ignore.covar=F), penalties=penalties.in["int",])
  plods.out<-c(plod.no,plod.add,plod.int)  
  cis<-data.frame()
  if(plod.no==max(plods.out)){
    model.out[["model"]]<-mod.no
    model.out[["s1maxlod"]]<-max(s1.no)
    model.out[["plod"]]<-plod.no
    model.out[["formula"]]<-"y ~ q1"
    ref<- refineqtl(cross, qtl=mod.no, formula="y ~ q1", pheno.col=pheno.in, covar=NULL, method="hk", verbose=F)
    ciout<-bayesint(ref,qtl.index=1, expandtomarkers=F, prob=0.95)
    chr.ci<-ref$chr
    lowmarker<-rownames(ciout)[1]
    highmarker<-rownames(ciout)[3]
    lowposition<-ciout[1,2]
    highposition<-ciout[3,2]
    cis<-rbind(cis,cbind(chr.ci,lowmarker,highmarker,lowposition,highposition))
    model.out[["fit"]]<-summary(fit.no)
    model.out[["cis"]]<-cis
    model.out[["ref"]]<-ref
  }else{
    if(plod.add==max(plods.out)){
      model.out[["model"]]<-mod.add
      model.out[["s1maxlod"]]<-max(s1.add)
      model.out[["plod"]]<-plod.add
      model.out[["formula"]]<-"y ~ trt + q1"
      ref<- refineqtl(cross, qtl=mod.add, formula="y ~ trt + q1", pheno.col=pheno.in, covar=covar.in, method="hk", verbose=F)
      ciout<-bayesint(ref,qtl.index=1, expandtomarkers=F, prob=0.95)
      chr.ci<-ref$chr
      lowmarker<-rownames(ciout)[1]
      highmarker<-rownames(ciout)[3]
      lowposition<-ciout[1,2]
      highposition<-ciout[3,2]
      cis<-rbind(cis,cbind(chr.ci,lowmarker,highmarker,lowposition,highposition))
      model.out[["fit"]]<-summary(fit.add)
      model.out[["cis"]]<-cis
      model.out[["ref"]]<-ref
    }else{
      if(plod.int==max(plods.out)){
        model.out[["model"]]<-mod.int
        model.out[["s1maxlod"]]<-max(s1.int)
        model.out[["plod"]]<-plod.int
        model.out[["formula"]]<-"y ~ trt + q1 + trt*q1"
        ref<- refineqtl(cross, qtl=mod.int, formula="y ~ trt + q1 + trt*q1", pheno.col=pheno.in, covar=covar.in, method="hk", verbose=F)
        ciout<-bayesint(ref,qtl.index=1, expandtomarkers=F, prob=0.95)
        chr.ci<-ref$chr
        lowmarker<-rownames(ciout)[1]
        highmarker<-rownames(ciout)[3]
        lowposition<-ciout[1,2]
        highposition<-ciout[3,2]
        cis<-rbind(cis,cbind(chr.ci,lowmarker,highmarker,lowposition,highposition))
        model.out[["fit"]]<-summary(fit.int)
        model.out[["cis"]]<-cis
        model.out[["ref"]]<-ref
      }
    }
  }
  names(model.out)<-c(paste("model", pheno.in, sep="."),
                      paste("s1maxlod", pheno.in, sep="."),
                      paste("plod", pheno.in, sep="."),
                      paste("formula", pheno.in, sep="."),
                      paste("fit", pheno.in, sep="."),
                      paste("cis", pheno.in, sep="."),
                      paste("ref", pheno.in, sep="."))
  return(model.out)
}


