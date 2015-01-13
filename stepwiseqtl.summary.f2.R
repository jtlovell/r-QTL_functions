#run stepwise qtl for your phenotypes... this is a wrapper that fits QTL and takes stats out. 
stepwiseqtl.summary.f2<-function(cross, phe, model.in, covar=NULL){
  #set up the environment
  stepout<-model.in
  if(nqtl(stepout)==0){
    return(NULL)
  }else{
    #plot the lod profile to the window
    #fit the qtl model
    nqtls<-nqtl(stepout)
    nterms<-sum(countqtlterms(formula(stepout), ignore.covar=F)[c(1,4)])
    ncovar<-length(covar)
    fit<-fitqtl(cross,
                pheno.col=phe,
                qtl=stepout,
                formula=formula(stepout),
                get.ests=T,dropone=T,covar=covar,
                method="hk")
    #part 2 output estimates
    ests.all<-summary(fit)$ests
    rows.dom<-grep("d",rownames(ests.all))
    rows.add<-grep("a",rownames(ests.all))
    if(nterms==1){
      dom.ests_out<-ests.all[rows.dom,1:3]; names(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
      add.ests_out<-ests.all[rows.add,1:3]; names(add.ests_out)<-c("est.add","SE.add","t.add")
      ests.out<-data.frame(t(data.frame(c(dom.ests_out,add.ests_out))))
    }else{
      if(nterms>nqtls){
        rows.epi<-grep(":",rownames(ests.all))
        rows.dom<-rows.dom[!rows.dom %in% rows.epi]
        rows.add<-rows.add[!rows.add %in% rows.epi]
        dom.ests_out<-data.frame(ests.all[rows.dom,1:3]); colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
        add.ests_out<-data.frame(ests.all[rows.add,1:3]); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
        ests.out<-cbind(dom.ests_out,add.ests_out)
        #add rows of nas to data frame for epis
        n.epi<-nterms-nqtls
        for (i in 1:n.epi) ests.out<-rbind(ests.out,NA)     
      }else{
        dom.ests_out<-data.frame(ests.all[rows.dom,1:3]); colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
        add.ests_out<-data.frame(ests.all[rows.add,1:3]); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
        ests.out<-cbind(dom.ests_out,add.ests_out)
      }
    }
    #part 3 get confidence intervals
    cis<-data.frame()
    #calculate confidence intervals for each qtl
    for (j in 1:nqtls){
      ciout<-lodint(stepout,qtl.index=j, expandtomarkers=F, drop=1.5)
      lowmarker<-rownames(ciout)[1]
      highmarker<-rownames(ciout)[3]
      lowposition<-ciout[1,2]
      highposition<-ciout[3,2]
      cis<-rbind(cis,cbind(lowmarker,highmarker,lowposition,highposition))
    }
    #add nas if epistasis
    if(nterms>nqtls){
      for (i in 1:n.epi) cis<-rbind(cis,NA)     
    }
    #part 4 full model stats out
    chrs<-data.frame(stepout$chr); colnames(chrs)<-"chr"
    if(nterms>nqtls){
      for (i in 1:n.epi) chrs<-rbind(chrs,NA)     
    }
    poss<-data.frame(stepout$pos); colnames(poss)<-"pos"
    if(nterms>nqtls){
      for (i in 1:n.epi) poss<-rbind(poss,NA)     
    }
    if(nqtls==1){
      dropone<-data.frame(t(data.frame(fit$result.full[1,c(1,2,4,5,7)])))
      rownames(dropone)<-paste(chrs,"@",round(poss,1), sep="")
    }else{
      dropone<-data.frame(fit$result.drop[,-c(5,6)])
      colnames(dropone)[2]<-"SS"
    }
    qtlids<-data.frame(rownames(dropone)); colnames(qtlids)<-"qtlid"
    stats<-cbind(phe,qtlids,chrs,poss,dropone,ests.out,cis)
    return(stats) 
  }
}
