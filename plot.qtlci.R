plot.qtlci<-function(  cross,
                       phenames,
                       chr,
                       pos.left,
                       pos.center,
                       pos.right,
                       cat,
                       legend.title){
  #get other custom functions in
  roundUp <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
  }
  

  toplot<-data.frame(phenames,chr,pos.left,pos.center,pos.right,cat)
  
  #make ordered numerical categorization for each phenotype
  toplot<-toplot[with(toplot, order(cat, phenames)), ]
  numnames<-as.numeric(toplot$phenames)
  toplot$phenames <- factor(toplot$phenames, 
                            levels = unique(toplot$phenames)[order(unique(toplot$phenames))])
  toplot$cat <- factor(toplot$cat, 
                       levels = unique(toplot$cat)[order(unique(toplot$cat))])
  phe1<-unique(toplot$phenames)
  num1<-length(unique(toplot$phenames)):1
  numnames<-sapply(toplot$phenames,function(x) num1[phe1==x])
  toplot$numnames<-numnames

  #make dataset for marker ticks
  marker.info<-data.frame(pull.map(cross, as.table=T))
  seg.ht<-(max(toplot$numnames)-1)*.05
  seg.start<-1-(max(toplot$numnames)*.1)
  seg.end<-seg.start-seg.ht
  marker.info$seg.end<-seg.end
  marker.info$seg.start<-seg.start
  
  #make dataset for axes info
  chr.ys<-length(num1)
  chr.len<-rep(as.numeric(chrlen(cross)),length(num1))
  chr<-rep(as.numeric(names(chrlen(cross))),length(num1))
  chr.beg<-rep(0,length(chr))
  chr.ys<-rep(num1,,each=nchr(cross))
  chr.info<-data.frame(chr.len,chr,chr.beg, chr.ys)
  
  #make the plot
  ggplot(marker.info)+
    geom_segment(aes(x=pos, xend = pos, y=seg.start, yend =seg.end))+ #lines for each confidence interval
    geom_segment(aes(x=chr.beg, xend = chr.len, y=chr.ys, yend =chr.ys),
                 alpha=.1, data=chr.info)+ #faint lines across each phenotype
    geom_point(aes(x=pos.center,y=numnames, col=cat), 
               data=toplot)+ #points for each QTL estimate
    geom_segment(aes(x=pos.left, xend=pos.right, yend=numnames,y=numnames,col=cat), 
                 data=toplot)+ #lines for each marker
    facet_grid(.~chr, scale="free_x",space="free_x")+ #split by chromosome

    theme_bw() +
    theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,panel.border = element_blank()
      ,strip.background = element_blank()
    ) +

    scale_color_discrete(name = legend.title)+
    ggtitle("chromosome")+
    theme(axis.line = element_line(color = 'black'))+
    #change axes
    scale_y_continuous("phenotypes",labels=phe1, breaks=num1)+
    scale_x_continuous("mapping position (cM)", 
                       breaks=seq(from=0,to=roundUp(max(marker.info$pos)),by=roundUp(max(marker.info$pos))/4), 
                       labels=c(0,roundUp(max(marker.info$pos))/4,roundUp(max(marker.info$pos))/2,"",""))
  
}
