# Author. JT Lovell
# Version. 2.1
# The plot.qtlci function takes QTL models and statistics and makes plots
# These plots are split by phenotype and chromosome to illustrate colocalization 

# The function takes four main datasets
# cross- the R/qtl cross object 
# toplot- a dataframe of statistics. This needs to include
#       phenotype names,
#       chromosomes,
#       point estimate (cM),
#       upper and lower confidence interval bounds (cM)
# models (if you want to make lod profile plots)-
#       this needs to be a list of QTL models that have associated LodProfiles
#       most easily generated from stepwiseQTL w/ keep.lodprofile=T

# In the QTL workflow associated with this function, these plots can directly utilize output from:
#       summary.stepwiseqt (part VI)-- this is the toplot file
#       stepwiseqtl(part III)-- this is the list of models for models

# The plottype argument asks if you want to plot:
#       lines: confidence interval lines around the point estimate
#       lodprofile: lodprofiles split across chromosomes and phenotypes
plot.qtlci<-function(  cross,
                       toplot,
                       phenames,
                       chr,
                       pos.left,
                       pos.center,
                       pos.right,
                       models=NULL,
                       plottype="lines",
                       standardize=TRUE){
  #get other custom functions in
  roundUp <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
  }
  library("scales")
  number_ticks <- function(n) {function(limits) pretty(limits, n)}
  
  toplot<-toplot[,c(phenames,chr,pos.left,pos.center,pos.right)]
  toplot$numnames<-as.numeric(toplot[,phenames])
  colnames(toplot)<-c("phe","chr","pos.left","pos","pos.right","numnames")
  phe1<-unique(toplot$phe)
  num1<-1:length(unique(toplot$phe))
  
  toplot$pos<-as.numeric(as.character(toplot$pos))
  toplot$pos.left<-as.numeric(as.character(toplot$pos.left))
  toplot$pos.right<-as.numeric(as.character(toplot$pos.right))
  
  #make dataset for marker ticks
  marker.info<-data.frame(pull.map(cross, as.table=T))[,1:2]
  colnames(marker.info)<-c("chr","pos")
  seg.end<-1-1/(max(toplot$numnames)*.25)
  seg.start<-seg.end-4/(max(toplot$numnames))
  marker.info$seg.end<-seg.end
  marker.info$seg.start<-seg.start
  marker.info$chr<-as.numeric(marker.info$chr)
  
  #make dataset for axes info
  chr.ys<-length(num1)
  if(class(cross)[1]=="4way"){
    chr.len<-rep(as.numeric(chrlen(cross)[1,]),length(num1))
  }else{
    chr.len<-rep(as.numeric(chrlen(cross)),length(num1))
  }
  chrs<-rep(as.numeric(chrnames(cross)),length(num1))
  chr.beg<-rep(0,length(chrs))
  chr.ys<-rep(num1,each=nchr(cross))
  chr.info<-data.frame(chr.len,chrs,chr.beg, chr.ys)
  colnames(chr.info)[2]<-"chr"
  
  #make the plot
  if(plottype=="lines"){
    ggplot(marker.info)+
      geom_rug(aes(x=pos),alpha=.2)+ #lines for each confidence interval
      geom_segment(aes(x=chr.beg, xend = chr.len, y=chr.ys, yend =chr.ys),
                   alpha=.05, data=chr.info)+ #faint lines across each phenotype
      geom_segment(aes(x=pos.left, xend=pos.right, yend=numnames,y=numnames), 
                   data=toplot, size=2, ,alpha=.5)+
      geom_point(aes(x=pos,y=numnames, color=phe),size=3, 
                 data=toplot)+ #points for each QTL estimate
      facet_grid(.~chr, scale="free_x",space="free_x")+ #split by chromosome
    
      theme(
        panel.background = element_rect(fill = "ghostwhite")
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,strip.background = element_blank()
        ,axis.text.x = element_blank()
        ,axis.ticks = element_blank()
        ,axis.text.y = element_text(colour="black")
      ) +
      theme(legend.position="none")+
      ggtitle("chromosome")+
      #change axes
      scale_y_continuous("phenotypes",labels=phe1, breaks=num1)+
      scale_x_continuous("")
  }else{
    if(plottype=="lodprofile"){
      #extract lod profiles
      lpdf<-data.frame()
      for.model<-intersect(as.character(unique(toplot$phe)), namesforlist)
      for(i in for.model){
        toext<-models[[i]]
        all.lps<-attr(toext, "lodprofile")
        std.lps<-sapply(all.lps,function(x) x$lod/max(x$lod))
        std.lps<-as.numeric(unlist(std.lps))
        lpdf.out <- ldply(all.lps, data.frame)
        colnames(lpdf.out)<-c("qtlname","chr","pos","lod.profile")
        lpdf.out$standardized.lod.profile<-std.lps
        lpdf.out$phenotype<-i
        lpdf<-rbind(lpdf,lpdf.out)
      }
      lpdf$chr <- factor(lpdf$chr, 
                                levels = 1:nchr(cross))
      if(standardize){
        marker.info$phenotype<-unique(lpdf$phenotype)[length(unique(lpdf$phenotype))]
        marker.info$phenotype<-unique(lpdf$phenotype)[length(unique(lpdf$phenotype))]
        marker.info$chr<-as.numeric(as.character(marker.info$chr))
        ggplot(lpdf)+
          geom_line(aes(x=pos,y=standardized.lod.profile, group=qtlname, color=phenotype))+
          geom_rug(data=marker.info,aes(x=pos),ticksize=.1,alpha=.2)+ #lines for each confidence interval
          facet_grid(phenotype~chr, scale="free",space="free_x")+
          theme(
            panel.background = element_rect(fill = "ghostwhite")
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,panel.border = element_blank()
            ,strip.background = element_blank()
            ,axis.text.x = element_blank()
            ,axis.ticks.x = element_blank()
            ,axis.text.y = element_text(colour="black")
          ) +
          theme(strip.text.y = element_text(angle = 0,hjust=0))+
          theme(legend.position="none")+
          ggtitle("chromosome")+
          #theme(axis.line = element_line(color = 'black'))+
          scale_y_continuous(breaks=number_ticks(2))+
          scale_x_continuous("")
      }else{
        marker.info$phenotype<-unique(lpdf$phenotype)[length(unique(lpdf$phenotype))]
        marker.info$chr<-as.numeric(as.character(marker.info$chr))
        lpdf$chr<-as.numeric(as.character(lpdf$chr))
        ggplot(lpdf)+
          geom_line(aes(x=pos,y=lod.profile, group=qtlname, color=phenotype))+
          geom_rug(data=marker.info,aes(x=pos),ticksize=.1,alpha=.2)+ #lines for each confidence interval
          facet_grid(phenotype~chr, scale="free",space="free_x")+
          theme(
            panel.background = element_rect(fill = "ghostwhite")
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,panel.border = element_blank()
            ,strip.background = element_blank()
            ,axis.text.x = element_blank()
            ,axis.ticks.x = element_blank()
            ,axis.text.y = element_text(colour="black")
          ) +
          theme(strip.text.y = element_text(angle = 0,hjust=0))+
          theme(legend.position="none")+
          ggtitle("chromosome")+
          #theme(axis.line = element_line(color = 'black'))+
          scale_y_continuous(breaks=number_ticks(2))+
          scale_x_continuous("")
      }
    }else{cat("do not know how to make plottype",plottype)}
  } 
}



