# cross: the R/qtl cross object
# stats: output from stepwise.summary.simple. if not provided, the following five vectors are neede
# phes=NULL,chrs=NULL, peak=NULL, right=NULL, left=NULL: vectors of phenotypes, chromosomes, the QTL positions (peak) and confidence interval bounds
# cols="black": line and point color, indexed by the line in stats, or position in vectors. Overridden by colbychr.
# pointsize=1: size (cex) of points, if FALSE, points are not plotted
# pointshape=1: shape (pch) of points
# linetype=1: lty of lines
# linethickness=1: lwd of lines
# mark.epi=FALSE: place a symbol where epistatic effects reside... not implemented yet
# plotQTLdensity=TRUE
# adj.ylabsize=TRUE: alter the size of the ylabel text based on the number of phenotypes?
# colbychr=TRUE: color points, lines and boxes by chromosome
# palette=terrain.colors: the palette to use for coloring
# outline=FALSE: put a box around the plot?
# background=FALSE: put a very transparent rectangle over each chromosome
# title="QTL position plot": title of the plot


plotMultiQTL<-function(cross, stats=NULL, phes=NULL,chrs=NULL, peak=NULL, right=NULL, left=NULL, cols=NULL, 
                       chr.subset=NULL, ylabelcex=NULL, rugsize=NULL,
                       pointsize=NULL, pointshape=19,linetype=1,linethickness=1,
                       plotQTLdensity=TRUE, binwidth=1, mark.epi=FALSE, 
                       colbychr=TRUE, palette=rainbow, 
                       outline=FALSE, background=TRUE, title="QTL position plot"){
  #set up env.
  if(!is.null(chr.subset)){
    cross<-subset(cross, chr=chr.subset)
    stats<-stats[stats$chromosome==chr.subset,]
  }
  pal<-palette(nchr(cross))
  if(is.null(stats)){
    stats<-data.frame(phes, chrs, peak, right, left)
    colnames(stats)<-c("phenotype","chromosome","position","lowCIpos","hiCIpos")
  }else{
    stats1<-stats[complete.cases(stats),c("phenotype","chromosome","position","lowCIpos","hiCIpos")]
    if(!is.null(cols)){
      cols<-cols[complete.cases(stats1)]
    }
    stats<-stats1
  }
  for(i in c("chromosome","position","lowCIpos","hiCIpos")){
    stats[,i]<-as.numeric(as.character(stats[,i]))
  }
  
  nqtl<-dim(stats)[1]
  nphes<-length(unique(stats$phenotype))
  if(is.null(rugsize)){
    rugsize<-(1/(nphes^2))+.01
  }
  if(is.null(pointsize)){
    pointsize<-(1/(1+(nphes*.01)))
  }
  if(is.null(ylabelcex)){
    ylab.adj<-(1/(1+(nphes*.015)))
  }
  a<-(sqrt(2*max(sapply(as.character(unique(stats$phenotype)),nchar)))*(ylab.adj^2))+4
  par(mar=c(5,a,4,2)+.1)
  map<-pull.map(cross, as.table=T)
  #determine plotting positions
  totlen<-sum(chrlen(cross))
  gapsize=totlen/40
  gaps<-c(0,cumsum(rep(gapsize,nchr(cross)-1)))
  chrslens<-c(0, cumsum(chrlen(cross)[-length(chrlen(cross))]))
  chrn<-chrnames(cross)
  corr<-gaps+chrslens
  corr<-data.frame(chrn, chrslens,corr, stringsAsFactors=F)
  
  #convert positions to plotting positions
  dat<-stats
  map2<-map
  for (i in corr$chrn){
    for(j in c("position", "lowCIpos", "hiCIpos")){
      dat[dat$chromosome==i,j]<-dat[dat$chromosome==i,j]+corr$corr[corr$chrn==i]
    }
    map2[map2$chr==i,"pos"]<-map2[map2$chr==i,"pos"]+corr$corr[corr$chrn==i]
  }
  
  
  if(plotQTLdensity){
    plot(0,0, ylim=c(0,(nphes+(.1*nphes)+1)), xlim=c(0,max(map2$pos)), type="n", bty="n", yaxt="n",xaxt="n",ylab="", xlab="Chromosome")
  }else{
    plot(0,0, ylim=c(0,nphes), xlim=c(0,max(map2$pos)), type="n", bty="n", yaxt="n",xaxt="n",ylab="", xlab="Chromosome")
  }
  if(outline) box()
  tloc<-ddply(map2,"chr",summarize,mean=mean(pos))$mean
  for(i in 1:nchr(cross)){
    axis(side=1, at=tloc[i], labels=chrnames(cross)[i])
  }
  phes<-as.character(unique(stats$phenotype))
  
  if(length(cols)==1) {
    cols<-rep(cols, nphes)
  }

  if(length(pointsize)==1) {
    pointsize<-rep(pointsize, nphes)
  }
    
  for(i in 1:nphes){
    p<-phes[i]
    tp<-dat[dat$phenotype==p,]
    if(is.null(ylabelcex)){
      axis(side=2, at=i, labels=p, las=2, cex.axis=ylab.adj)
    }else{
      axis(side=2, at=i, labels=p, las=2, cex.axis=ylabelcex)
    }
    
    for(j in 1:nrow(tp)){
      if(colbychr & is.null(cols)){
        segments(tp$lowCIpos[j],i,tp$hiCIpos[j],i,col=pal[which(tp$chr[j]==chrnames(cross))], lty=linetype, lwd=linethickness)
      }else{
        segments(tp$lowCIpos[j],i,tp$hiCIpos[j],i,col=cols[i], lty=linetype, lwd=linethickness)
      }

    }
    if(pointshape){
      if(colbychr & is.null(cols)){
        points(tp$position,rep(i,length(tp$position)), col=pal[tp$chr],cex=pointsize[i], pch=pointshape)
      }
      else{
        points(tp$position,rep(i,length(tp$position)), col=cols[i],cex=pointsize[i], pch=pointshape)
      }
    }
  }
  rug(map2$pos, ticksize=rugsize)
  if(background){
    if(colbychr){
      min.pos<-tapply(map2$pos, map2$chr, min)
      max.pos<-tapply(map2$pos, map2$chr, max)
      pal2<-palette(nchr(cross), alpha=.05)
      rect(min.pos,rep(0,nchr(cross)),max.pos, rep(nphes), col=pal2, border = pal2)
    }else{
      min.pos<-tapply(map2$pos, map2$chr, min)
      max.pos<-tapply(map2$pos, map2$chr, max)
      rect(min.pos,rep(0,nchr(cross)),max.pos, rep(nphes), col=rgb(0,0,0,.05), border = rgb(0,0,0,.05))
    }
  }
  title(title)
  if(plotQTLdensity){
    denpos<-density(dat$position,bw=binwidth)$x
    den<-density(dat$position,bw=binwidth)$y
    den<-den/max(den)
    den<-den*(.1*nphes)+1
    den<-den+nphes
    lines(denpos,den)
  }
}
  