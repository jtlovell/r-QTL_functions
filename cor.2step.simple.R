cor.2step_simple<-function(test){
  for(j in seqin){
    dat<-data.matrix(all.snps[,(j-1):j])
    c<-cor(dat, use="pairwise.complete.obs")[1,2]
    if(c==1){
      na1<-length(dat[,1][is.na(dat[,1])])
      na2<-length(dat[,2][is.na(dat[,2])])
      if(na1==na2){
        index<-1
      }else{
        index<-which(c(na1,na2)==min(c(na1,na2)))[1] 
      }
      out<-colnames(dat)[index]
    }else{
      out<-colnames(dat)
    }
    out.temp<-c(out.temp,out)
    print(j)
  }  
}