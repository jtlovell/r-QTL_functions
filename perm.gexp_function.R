# perm gene expression function

perm.gexp<-function(
  geneexp.culled=geneexp.culled,
  idnames=idnames,
  model.info=model.info,
  cov.chr=cov.chr,
  genoprobs.ascov=genoprobs.ascov,
  cross=cross,
  gene.in=gene.in){
  
  boot.out<-data.frame()
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
  return(boot.out)
}

