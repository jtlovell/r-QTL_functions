#this function returns the best model from a set of formulae
#used to select a complex model that includes epistasis and interactions with the environment
add.1<-function(cross, qtl.in=model.in, formula.in=formulas, phenos.in, 
                covar=cov,formula.types=formula.types, penalties.in=penalties.in, data.model="normal"){
  output<-list()
  bestmodels<-list()
  formula.out<-list()
  model.fit<-list()
  cis<-list()
  plods<-vector()
  for (i in 1:length(formula.in)){
    print(i)
    scan <- addqtl(cross, qtl=qtl.in, formula=formula.in[[i]], method="hk", covar=covar, pheno.col=phenos.in, model=data.model)
    mod <- addtoqtl(cross, qtl.in, max(scan)$chr, max(scan)$pos)
    ref<-refineqtl(cross, mod, formula=formula.in[[i]], pheno.col=phenos.in, covar=covar, method="hk", model=data.model)
    fit <- fitqtl(cross, qtl=ref, formula=formula.in[[i]], pheno.col=phenos.in, covar=covar, method="hk", model=data.model)
    
    cis<-data.frame()
    for (j in 1:countqtlterms(formula.in[[i]], ignore.covar=T)){
      ciout<-bayesint(ref,qtl.index=j, expandtomarkers=F, prob=0.95)
      chr.ci<-ref$chr[j]
      lowmarker<-rownames(ciout)[1]
      highmarker<-rownames(ciout)[3]
      lowposition<-ciout[1,2]
      highposition<-ciout[3,2]
      cis<-rbind(cis,cbind(chr.ci,lowmarker,highmarker,lowposition,highposition))
    }
    lod<-fit$result.full[1,4]
    bestmodels[[i]]<-mod
    formula.out[[i]]<-formula.in[[i]]
    if(formula.types[i]=="no"){
      plods[[i]]<-calc.plod(lod, countqtlterms(formula.in[[i]], ignore.covar=T), penalties=penalties.in["no",])
    }else{
      if(formula.types[i]=="add"){
        plods[[i]]<-calc.plod(lod, countqtlterms(formula.in[[i]], ignore.covar=T), penalties=penalties.in["add",])
      }else{
        plods[[i]]<-calc.plod(lod, countqtlterms(formula.in[[i]], ignore.covar=T), penalties=penalties.in["int",])
      }
    }
  } 
  
  output[["model"]]<-bestmodels[[which(plods==max(plods))[1]]]
  output[["plod"]]<-plods[which(plods==max(plods))[1]]
  output[["formula"]]<-formula.out[which(plods==max(plods))[1]]
  
  output[["fit"]]<-summary(fit)
  output[["cis"]]<-cis
  output[["ref"]]<-ref
  
  names(output)<-c(paste("model", phenos.in, sep="."),
                   paste("plod", phenos.in, sep="."),
                   paste("formula", phenos.in, sep="."),
                   paste("fit", phenos.in, sep="."),
                   paste("cis", phenos.in, sep="."),
                   paste("ref", phenos.in, sep="."))
  return(output)
}