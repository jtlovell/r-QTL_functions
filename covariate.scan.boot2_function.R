#covariate scan function with bootstrapping to determine significance. Extremely computationally intensive.
#if possible just run for 1 QTL with very few genes... ~3minutes/gene, so if you have 10 genes in the QTL region, 30 minutes to run
#>1000 bootstraps needed to achieve convergence. 
#most of this function is preparing the covariate matrix for each type of QTL
#right now, this can incorporate all types of qtl, but these need to be specified
#this is not optimized for other environments- proline.dry is treated specially 
covariate.scan.boot2<-function(
  beststats,
  bestmodels,
  cross,
  cyto,
  idnames,
  bestmodel.types,
  candidates.byeqtl,
  qtl.ids,
  nboot,
  interactor.df){
  
  cross.sub<-subset(cross,ind=idnames$id)
  cyto.covar<-data.frame(cyto$cyto.num[cyto$id %in% idnames$id]); colnames(cyto.covar)<-"cyto"
  all.stats.out<-data.frame()
  all.out<-data.frame()
  for (i in qtl.ids){
    qtl.id.in<-i
    cat("\n","running covariate scans for genes in", i, sep=" ")
    pheno.in<-as.character(beststats$phenotype[beststats$qtl.id==qtl.id.in])
    chr.in<-as.numeric(beststats$chromosome[beststats$qtl.id==qtl.id.in])
    pos.in<-as.numeric(beststats$position[beststats$qtl.id==qtl.id.in])
    mar.in<-find.marker(cross.sub,chr=chr.in,pos=pos.in)
    trt.in<-as.character(beststats$phenotype[beststats$qtl.id==qtl.id.in])
    trt.in<-ifelse(strsplit(trt.in,"[.]")[[1]][2]=="dry","dry","wet")
    trt.in[is.na(trt.in)]<-"wet"
    cov.chr<-summary(bestmodels[[pheno.in]])$chr
    cov.pos<-summary(bestmodels[[pheno.in]])$pos
    model.info<-as.character(bestmodel.types$best.model.type[rownames(bestmodel.types)==pheno.in])
    
    #get the name of the marker under the qtl and remove it
    if(length(cov.chr)>1){
      cov.markers<-find.marker(cross.sub,chr=cov.chr,pos=cov.pos)
      cov.markers<-cov.markers[!cov.markers==mar.in]
      cov.markers.topull<-paste(cov.markers,"BB",sep=":")
      genoprobs<-pull.genoprob(cross.sub,omit.first.prob=T)
      genoprobs.ascov<-as.data.frame(genoprobs[,cov.markers.topull])
      if (model.info=="add.covar"){
        genoprobs.ascov<-cbind(genoprobs.ascov,cyto.covar)
      }
    }
    if(model.info=="add.covar" & length(cov.chr)==1){
      genoprobs.ascov<-cyto.covar
    }
    if (model.info=="no.covar"& length(cov.chr)==1){
      genoprobs.ascov<-NULL
    }    
    #epistasis into model
    if (i %in% interactor.df$QTL.1){
      interactor.qtl.in<-interactor.df$QTL.2[interactor.df$QTL.1== qtl.id.in]
      epi.qtl.stats<-beststats[beststats$qtl.id==interactor.qtl.in,]
      marker<-find.marker(cross.sub,chr=epi.qtl.stats$chromosome, pos=epi.qtl.stats$position)
      genotype.intcov<-pull.genoprob(cross.sub, chr=epi.qtl.stats$chromosome, omit.first.prob=T)
      genotype.intcov<-as.data.frame(genotype.intcov[,paste(marker,":BB", sep="")]); colnames(genotype.intcov)<-marker
    }else{genotype.intcov<-NULL}
    s1.full<-scanone(cross.sub,pheno.col=pheno.in,method="hk",addcovar=genoprobs.ascov, intcovar=genotype.intcov)
    s1.est<-as.data.frame(s1.full)
    s1.est<-s1.est[rownames(s1.est)==mar.in,"lod"]
    if(qtl.id.in %in% focal.qtls){
      genes.to.run<-as.character(top.cov.out.epi$gene.in[top.cov.out.epi$qtl.id.in==qtl.id.in])
    }else{
      if(qtl.id.in=="proline.dry_2_74"){
        genes.to.run<-proline.dry.candidates.wcyto
      }else{
        genes.to.run<-as.character(candidates.byeqtl$all.eqtl.genes[candidates.byeqtl$qtl==qtl.id.in])
      }
    }
    genes.to.run<-genes.to.run[!is.na(genes.to.run)]
    colnames.in<-c("qtl.id.in","sample","gene.in","pheno.in","chr.in","pos.in","mar.in","trt.in","model.info","s1.est","s1.cov.est")
    gene.out<-data.frame()
    stats.out<-data.frame()
    for (j in 1:length(genes.to.run)){
      nsamples=nboot
      gene.in<-genes.to.run[j]
      geneexp.culled<-geneexp[,c("treatment","line",gene.in)]
      geneexp.culled<-geneexp.culled[geneexp.culled$treatment==trt.in,]
      names(idnames)[2]<-"line"
      exp.covars<-merge(idnames,geneexp.culled, by="line", all=T)
      exp.cov<-data.frame(exp.covars[,4])
      if(model.info=="no.covar" & length(cov.chr)==1){
        all.cov<-exp.cov
      }else{
        all.cov<-cbind(genoprobs.ascov,exp.cov)
      }
      
      s1.cov<-scanone(cross.sub,pheno.col=pheno.in,method="hk",addcovar=all.cov, intcovar=genotype.intcov)
      s1.cov.est<-as.numeric(s1.cov[rownames(s1.cov)==mar.in,"lod"])
      bounds<-as.numeric(beststats[beststats$qtl.id==qtl.id.in,][c("lowCIpos","hiCIpos")])
      boot.out<-data.frame()
      hist.loc<-c(pos.in,seq(from=20, to=chrlen(cross.sub)[chr.in], by=20))
      cat("\n","running",nsamples, "bootstraps for", gene.in, sep=" ")
      
      for (k in 1:nsamples){
        geneexp.culled.samp<-geneexp.culled
        geneexp.culled.samp[,gene.in]<-sample(geneexp.culled.samp[,gene.in], replace=T)
        names(idnames)[2]<-"line"
        exp.covars.samp<-merge(idnames,geneexp.culled.samp, by="line", all=T)
        exp.cov.samp<-data.frame(exp.covars.samp[,4])
        if(model.info=="no.covar"& length(cov.chr)==1){
          all.cov<-exp.cov.samp
        }else{
          all.cov<-cbind(genoprobs.ascov,exp.cov.samp)
        }
        s1.cov.samp<-scanone(cross.sub,pheno.col=pheno.in,method="hk",addcovar=all.cov, intcovar=genotype.intcov)
        s1.cov.est.samp<-as.data.frame(s1.cov.samp)
        s1.cov.est.samp<-as.numeric(s1.cov.est.samp[rownames(s1.cov.est.samp)==mar.in ,"lod"])
        output.samp<-data.frame(qtl.id.in,paste("samp",k,sep=""),gene.in,pheno.in,chr.in,pos.in,mar.in,trt.in,model.info,s1.est,s1.cov.est.samp)
        colnames(output.samp)<-colnames.in
        boot.out<-rbind(boot.out,output.samp)
        samp.diff<-s1.full-s1.cov.samp
      }
      output.cov<-data.frame(qtl.id.in,"main.covscan",gene.in,pheno.in,chr.in,pos.in,mar.in,trt.in,model.info,s1.est,s1.cov.est)
      colnames(output.cov)<-colnames.in
      gene.out<-rbind(gene.out,output.cov,boot.out)
      pval.out<-length(boot.out$s1.cov.est[boot.out$s1.cov.est<s1.cov.est])/nboot
      effect.out<-(s1.est-s1.cov.est)/s1.est
      cat("... P-value= ", pval.out, "; effect= ", effect.out)
      stats<-data.frame(qtl.id.in,"main.covscan",gene.in,pheno.in,chr.in,pos.in,mar.in,trt.in,model.info,pval.out,effect.out)
      stats.out<-rbind(stats.out,stats)
    }
    all.stats.out<-rbind(all.stats.out,stats.out)
  }
  return(all.stats.out)
}

