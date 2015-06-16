stepwiseqtl.summary.simple.4way<-function(step, phename=1, cross=fake.4way, covar=NULL){
  #fit the qtl model
  stepout<-step
  fit<-fitqtl(cross,
              pheno.col=phename,
              qtl=stepout,
              formula=formula(stepout),
              get.ests=F,dropone=T,covar=covar,
              method="hk")
  #determine the number of terms/qtl in the model
  nqtls<-nqtl(stepout)
  nterms<-sum(countqtlterms(formula(stepout), ignore.covar=F)[c(1,4)])
  ncovar<-length(covar)
  #extract information if the qtl model has multiple qtls/covariates
  if(nterms==1){
    ciout<-lodint(stepout,qtl.index=1, expandtomarkers=F, drop=1.5)
    lowmarker<-rownames(ciout)[1]
    highmarker<-rownames(ciout)[3]
    lowposition<-ciout[1,2]
    highposition<-ciout[3,2]
    cis<-data.frame(lowmarker,highmarker,lowposition,highposition)
    chrs<-stepout$chr
    pos<-stepout$pos
    fit.out<-data.frame(summary(fit)$result.full)[1,c(1,2,4,5,3,6,7)]
    fit.out[1,5]<-NA
    colnames(fit.out)<-c("df","Type.III.SS","LOD","X.var","F.value","Pvalue.Chi2.","Pvalue.F.")
    stats<-data.frame(rep(phename,nterms),chrs,pos,fit.out,cis)
  }else{
    if(nterms>1){
      cis<-data.frame()
      #calculate confidence intervals for each qtl
      for (j in 1:nqtls){
        attr(stepout, "lodprofile")
        ciout<-lodint(stepout,qtl.index=j, expandtomarkers=F, drop=1.5)
        lowmarker<-rownames(ciout)[1]
        highmarker<-rownames(ciout)[3]
        lowposition<-ciout[1,2]
        highposition<-ciout[3,2]
        cis<-rbind(cis,cbind(lowmarker,highmarker,lowposition,highposition))
      }
      chrs<-stepout$chr
      pos<-stepout$pos
      #if covariates, add in NAs into the dataframe
      if(ncovar>0) {
        covar.ci<-(rep(NA,4))
        cis<-rbind(covar.ci, cis)
        chrs<-c(names(covar), stepout$chr)
        pos<- c(names(covar), stepout$pos)
      }
      #if epistasis, add NAs into the dataframe
      if(as.numeric(countqtlterms(formula(stepout))[4])>0) {
        epi.ci<-data.frame()
        for (j in 1:countqtlterms(formula(stepout))[4]){
          epi.ci.out<-(rep(NA,4))
          cis<-rbind(cis,epi.ci.out)
        }
        chrs<-c(chrs,rep(NA,countqtlterms(formula(stepout))[4]))
        pos<- c(pos,rep(NA,countqtlterms(formula(stepout))[4]))
      }
      stats<-data.frame(rep(phename,nterms),chrs,pos,fit$result.drop[1:nterms,],cis)
    }
  }
  colnames(stats)[1:2]<-c("phe","chr")
  return(stats)
}
