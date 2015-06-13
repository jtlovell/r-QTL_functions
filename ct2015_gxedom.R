#load datasets
#final.stats
stats<-read.csv("ct2015_finalstats.csv")
stats.cull<-stats[,c("id","best.mode","fig2_colorcat","trt_pvalue","cis_sig","trans_sig","trt.cis_sig","trt.trans_sig","cis_propvar","trt_propvar","trans_propvar")]

#dominance pvalues
load("ct2015_dominancecatuntransformed.RData")
res.parents<-dominance
for(i in colnames(res.parents)[-11]) res.parents[,i]<-as.numeric(unlist(res.parents[,i]))
res.parents$id<-rownames(res.parents)
modes<-c("codomominant","dominant","recessive","overdominant","log.additive")
res.parents$best.mode<-apply(res.parents[,6:10],1,function(x) modes[which(x==min(x))][1])
res.parents$min.p<-apply(res.parents[,1:5],1,function(x) x[which(x==min(x))][1])
res.parents$best.mode[res.parents$min.p>0.05]<-"ambiguous"
dom<-res.parents[,c("id",colnames(res.parents)[grep("pvalue",colnames(res.parents))])]

#g+e+gxe pvalues
gxe<-read.csv("ct2015_gxestats.csv")
gxe.cull<-gxe[,c("id","trt.pvalue","geno.pvalue","genotrt.pvalue","main.sig")]
colnames(gxe.cull)[5]<-"gxe.sig"

#merge
dat<-merge(stats.cull, res.parents, by="id")
dat<-merge(dat,gxe.cull, by="id")

#xy plots of pvalue ranks
with(dat, plot(rank(-log10(Dominant.pvalue)),rank(-log10(trt_pvalue)), pch="."))

#barplots with proportions and raw counts
pdf("test.pdf")
par(mfrow=c(2,2))
tab1<-table(dat[,c("fig2_colorcat","best.mode")])
tx<-barplot(tab1[,c(1,2,4,3,6,5)], col=c("grey","red","blue","green","orange"), xaxt="n")
text(cex=1, x=x+.1, y=-500, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("counts of genes (trt in model)")
legend("topright",c("conserved","dry only","wet only","both opp.","both same"), fill=c("grey","red","blue","green","orange"))

tab1n<-apply(tab1,2,function(x) x/sum(x))
x<-barplot(tab1n[,c(1,2,4,3,6,5)], col=c("grey","red","blue","green","orange"), xaxt="n")
text(cex=1, x=x+.1, y=-.05, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("proportions of fig2 cats")

tab1<-table(dat[,c("gxe.sig","best.mode")])
x<-barplot(tab1[,c(1,2,4,3,6,5)], col=c("grey","red","blue","green","orange"), xaxt="n")
text(cex=1, x=x+.1, y=-500, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("counts of genes (by GxE sig)")
legend("topright",c("","E","E+G","G","GxE"), fill=c("grey","red","blue","green","orange"))

tab1n<-apply(tab1,2,function(x) x/sum(x))
x<-barplot(tab1n[,c(1,2,4,3,6,5)], col=c("grey","red","blue","green","orange"), xaxt="n")
text(cex=1, x=x+.1, y=-.05, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("proportions  (by GxE sig)")

tab1<-table(dat[,c("cis_sig","best.mode")])
x<-barplot(tab1[,c(1,2,4,3,6,5)], col=c("grey","red","blue","green","orange"), xaxt="n")
text(cex=1, x=x+.1, y=-500, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("counts of genes (by cis sig)")
legend("topright",c("not sig","cis P<0.05"), fill=c("grey","red"))

tab1n<-apply(tab1,2,function(x) x/sum(x))
x<-barplot(tab1n[,c(1,2,4,3,6,5)], col=c("grey","red","blue","green","orange"), xaxt="n")
text(cex=1, x=x+.1, y=-.05, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("proportions (by cis sig)")

tab1<-table(dat[,c("trans_sig","best.mode")])
x<-barplot(tab1[,c(1,2,4,3,6,5)], col=c("grey","red","blue","green","orange"), xaxt="n")
text(cex=1, x=x+.1, y=-500, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("counts of genes (by trans sig)")
legend("topright",c("not sig","trans p<0.05"), fill=c("grey","red"))

tab1n<-apply(tab1,2,function(x) x/sum(x))
x<-barplot(tab1n[,c(1,2,4,3,6,5)], col=c("grey","red","blue","green","orange"), xaxt="n")
text(cex=1, x=x+.1, y=-.05, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("proportions (by trans sig)")


dev.off()



dom.genes<-as.character(dat$id[dat$best.mode=="dominant"])
rec.genes<-as.character(dat$id[dat$best.mode=="recessive"])
c<-data.frame(t(counts(dds.all3)))
d<-data.frame(colData(dds.all3))
cd<-c[,dom.genes]
cr<-c[,rec.genes]
rd<-cbind(d,cd)
rr<-cbind(d,cr)
test<-  ddply(rd, c("Treatment","id"), numcolwise(mean))
test1<-test[,-which(colnames(test) %in% c("VWC_July","time.sampled"))]
test2<-data.frame(test1[,1:2],apply(test1[,dom.genes],2,scale))

test<-  ddply(rr, c("Treatment","id"), numcolwise(mean))
test1<-test[,-which(colnames(test) %in% c("VWC_July","time.sampled"))]
test2<-data.frame(test1[,1:2],apply(test1[,rec.genes],2,scale))

test2$id<-as.character(test2$id)
test2$Treatment<-as.character(test2$Treatment)


dat1<-test2[,]
pdf("ct2015_lineplots_recscaled.pdf")
par(mfrow=c(3,2))
plot(dat1[1:3,4], ylim=c(min(dat1[,-c(1:2)]),max(dat1[,-c(1:2)])), xlim=c(.8,3.2), type="n", bty="n", axes=F, ylab="scaled exp", xlab="genotype")
axis(2, at=c(-2,-1,0,1,2)); axis(1,at=c(1,2,3),labels=c("Fil","F1","Hal"))
count<-1
for(i in 4:dim(dat1)[2]){
  p<-dat1[,i]
  if(min(p)==p[1]){
    lines(p[c(2,1,3)], col=rgb(1,0,0,alpha=.2));lines(p[c(5,4,6)], col=rgb(0,0,1,alpha=.2))
    count<-count+1
  }
}
title(paste("F1 dry is lowest .. n=", count))

plot(dat1[1:3,4], ylim=c(-2,2), xlim=c(.8,3.2), type="n", bty="n", axes=F, ylab="scaled exp", xlab="genotype")
axis(2, at=c(-2,-1,0,1,2)); axis(1,at=c(1,2,3),labels=c("Fil","F1","Hal"))
count<-1
for(i in 4:dim(dat1)[2]){
  p<-dat1[,i]
  if(min(p)==p[4]){
    lines(p[c(2,1,3)], col=rgb(1,0,0,alpha=.2));lines(p[c(5,4,6)], col=rgb(0,0,1,alpha=.2))
    count<-count+1
  }
}
title(paste("F1 wet is lowest .. n=", count))

plot(dat1[1:3,4], ylim=c(-2,2), xlim=c(.8,3.2), type="n", bty="n", axes=F, ylab="scaled exp", xlab="genotype")
axis(2, at=c(-2,-1,0,1,2)); axis(1,at=c(1,2,3),labels=c("Fil","F1","Hal"))
count<-1
for(i in 4:dim(dat1)[2]){
  p<-dat1[,i]
  if(min(p)==p[2]){
    lines(p[c(2,1,3)], col=rgb(1,0,0,alpha=.2));lines(p[c(5,4,6)], col=rgb(0,0,1,alpha=.2))
    count<-count+1
  }
}
title(paste("Fil dry is lowest .. n=", count))

plot(dat1[1:3,4], ylim=c(-2,2), xlim=c(.8,3.2), type="n", bty="n", axes=F, ylab="scaled exp", xlab="genotype")
axis(2, at=c(-2,-1,0,1,2)); axis(1,at=c(1,2,3),labels=c("Fil","F1","Hal"))
count<-1
for(i in 4:dim(dat1)[2]){
  p<-dat1[,i]
  if(min(p)==p[5]){
    lines(p[c(2,1,3)], col=rgb(1,0,0,alpha=.2));lines(p[c(5,4,6)], col=rgb(0,0,1,alpha=.2))
    count<-count+1
  }
}
title(paste("Fil wet is lowest .. n=", count))

plot(dat1[1:3,4], ylim=c(-2,2), xlim=c(.8,3.2), type="n", bty="n", axes=F, ylab="scaled exp", xlab="genotype")
axis(2, at=c(-2,-1,0,1,2)); axis(1,at=c(1,2,3),labels=c("Fil","F1","Hal"))
count<-1
for(i in 4:dim(dat1)[2]){
  p<-dat1[,i]
  if(min(p)==p[3]){
    lines(p[c(2,1,3)], col=rgb(1,0,0,alpha=.2));lines(p[c(5,4,6)], col=rgb(0,0,1,alpha=.2))
    count<-count+1
  }
}
title(paste("Hal dry is lowest .. n=", count))


plot(dat1[1:3,4], ylim=c(-2,2), xlim=c(.8,3.2), type="n", bty="n", axes=F, ylab="scaled exp", xlab="genotype")
axis(2, at=c(-2,-1,0,1,2)); axis(1,at=c(1,2,3),labels=c("Fil","F1","Hal"))
count<-1
for(i in 4:dim(dat1)[2]){
  p<-dat1[,i]
  if(min(p)==p[6]){
    lines(p[c(2,1,3)], col=rgb(1,0,0,alpha=.2));lines(p[c(5,4,6)], col=rgb(0,0,1,alpha=.2))
    count<-count+1
  }
}
title(paste("Hal wet is lowest .. n=", count))
dev.off()


dat2<-test2[,c(1,2,which(type=="haldryup"))]
dat1<-melt(dat1,id=c("Treatment","id"))
dat1$value<-as.numeric(dat1$value)
dat1$id<-factor(dat1$id, levels=c("HAL2","F1","FIL2"))
dat2<-melt(dat2,id=c("Treatment","id"))
dat2$id<-factor(dat2$id, levels=c("HAL2","F1","FIL2"))
ggplot(dat1,aes(x=id, y=value, col=Treatment, group=interaction(Treatment,variable)))+
  geom_line()
ggplot(dat1,aes(x=id, y=value, col=Treatment, group=interaction(Treatment,variable)))+
  geom_line()

# plot with qualitative labels on y-axis
out<-res.parents[,c("id","Dominant.pvalue","Recessive.pvalue","Overdominant.pvalue","best.mode")]
stats.cull<-stats[,c("best.mode","fig2_colorcat")]
fig2a.cull<-

# plot with quantitative y-axis
iplotMScanone(out, times=times)

table(dom2[,c("best.mode","best.mode.trtIgnored")])
pdf("ct2015_countsbyinher.pdf")


tab1<-table(dom2[,c("fig2_colorcat","best.mode.trtIgnored")])
x<-barplot(tab1[,c(1,2,4,3,6,5)], col=c("grey","red","blue","green","orange"), xaxt="n")
text(cex=1, x=x+.1, y=-500, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("counts of genes (trt ignored)")

tab1n<-apply(tab1,2,function(x) x/sum(x))
x<-barplot(tab1n[,c(1,2,4,3,6,5)], col=c("grey","red","blue","green","orange"), xaxt="n")
text(cex=1, x=x+.1, y=-.05, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("proportions by fig2 cat (trt ignored)")

tab1<-table(dom2[,c("main.sig","best.mode.trtIgnored")])
x<-barplot(tab1[,c(1,2,4,3,6,5)], col=c("pink","lightgreen","pink","blue","lightgreen","grey","green","blue","blue","grey"), xaxt="n")
text(cex=1, x=x+.1, y=-500, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("counts of genes (trt ignored)")

tab1n<-apply(tab1,2,function(x) x/sum(x))
x<-barplot(tab1n[,c(1,2,4,3,6,5)], col=c("pink","lightgreen","pink","blue","lightgreen","grey","green","blue","blue","grey"), xaxt="n")
text(cex=1, x=x+.1, y=-.05, colnames(tab1[,c(1,2,4,3,6,5)]), xpd=TRUE, srt=30, pos=2)
title("proportions by fig2 cat (trt ignored)")

dev.off()
