qtl2df<-function(cross,chrs, poss, phe){
  mar<-find.pseudomarker(cross, chr=chrs, pos=poss)
  gp<-pull.genoprob(fake.4way, chr=chrs)
  gp<-gp[,grep(mar,colnames(gp))]
  
  als<-as.character(sapply(colnames(gp),function(x) strsplit(x,":")[[1]][2]))
  gp<-data.frame(gp)
  gp$call<-as.character(apply(gp,1,function(x)  ifelse(max(x)>.95, als[which(x==max(x))], NA)))
  calls<-as.character(gp$call)
  phe1<-pull.pheno(cross, pheno.col=phe)
  id<-cross$phe$id
  df<-data.frame(id,calls,phe1)
  colnames(df)[1:3]<-c("id",paste("chr",chrs, ".pos",poss, sep=""),phe)
  return(df)
}

qtl2tuk<-function(cross, chrs, poss, phe, covar=NULL, wt=NULL, alpha=0.05){
  n.qtls<-length(chrs)
  out<-list()
  for(i in 1:n.qtls){
    chr=chrs[i]; pos=poss[i]
    out[[i]]<-qtl2df(cross=cross, chrs=chr, poss=pos, phe=phe)
  }
  df<-Reduce(function(...) merge(..., by=c("id",phe),all=T), out)
  if(!is.null(covar)){
    df<-data.frame(df,covar)
  }
  qtls<-colnames(df)[-which(colnames(df) %in% c("id",phe))]
  form<-as.formula(paste(phe,"~",paste(qtls,collapse="+")))
  av<-lm(form,data=df)
  cld.out<-list()
  for(i in qtls){
    if(is.null(wt)){
      lsm<-lsmeans(av,specs=i)
    }else{
      lsm<-lsmeans(av,specs=i, weights=wt)
    }
    if(length(unique(df[,i]))<3){
      cld.out[[i]]<-lsm
    }else{
      cld.out[[i]]<-cld(lsm, alpha=alpha)
    }
  }
  cld.out
}

#implementation on a single trait
data(fake.4way)
fake.4way <- sim.geno(fake.4way, step=1)
fake.4way <- calc.genoprob(fake.4way, step=1)
id<-1:nind(fake.4way)
fake.4way$pheno<-cbind(fake.4way$pheno, id)
test<-qtl2tuk(cross=fake.4way, n.qtls=2, chrs=c(2,7), poss=c(9,42), phe="phenotype")

#looping though multiple traits
#stepout should be a dataframe with all epistatic elements removed. 
out<-list()
for (i in unique(stepout$phenotype)){
  dat<-stepout[stepout$phenotype==i,]
  chrs=dat$chr
  poss=dat$pos
  out[[i]]<-qtl2tuk(cross=cross, n.qtls=length(chrs), chrs=chrs, poss=poss, phe=i)
}

mname1<-find.marker(fake.4way, 2,9)
out<-eft(fake.4way, pheno.col="phenotype", mname1=mname1)
>>>>>>> upstream/master
