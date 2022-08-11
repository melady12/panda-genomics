library(plyr)
library(gtools)
library(vioplot)
library(viridis)
library(scales)

d1<-read.table("/Users/melady/Desktop/PANDA/Projects/03.panda.assembly/04.AB-compartment/JJ.whole.100Kb.AB.compartment.txt",stringsAsFactors=F)
d2<-read.table("/Users/melady/Desktop/PANDA/Projects/03.panda.assembly/04.AB-compartment/JJ.genome.100kb.TSS.count.chr.bed",stringsAsFactors=F)
names(d1)[4]<-'pc1'
names(d2)[4]<-'tss'
dd<-merge(d1,d2,by=c('V1','V2','V3'))
dd$type<-ifelse(dd$tss>6,">6",dd$tss)

lev<-mixedsort(unique(dd$type))

#col<-viridis(length(lev))
#col<-col_numeric(c('#a9252d','#2969a6'),range(1:length(lev)))(1:length(lev))
col<-col_numeric(c('#5585B6','#E24F98'),range(1:length(lev)))(1:length(lev))
n<-c()
med<-c()
mea<-c()

pdf(paste("PC1.and.genedensity.pdf",sep=""),w=4,h=4)
par(xaxs='i',yaxs='i',lend=2,mgp=c(2,0.6,0),tcl=-0.3,mar=c(5,4,5,2),bty='l')
plot(1:10,type='n',axes=FALSE,xlab="Number of PCGs in each 100Kb bin",ylab="PC1 values",xlim=c(0.3,length(lev)+0.5),ylim=c(-40,40))

for (j in seq(lev)){
	idat<-dd[dd$type==lev[j],]
	vioplot(idat$pc1,at=j,col=col[j],border=col[j],add=TRUE)
	n<-c(n,length(idat$pc1))
	med<-c(med,median(idat$pc1))
	mea<-c(mea,mean(idat$pc1))
}
smed<-sprintf("%.2f",med)
smea<-sprintf("%.2f",mea)
maxy<-max(dd$pc1,na.rm=TRUE)*1.7
unit<-diff(range(idat$pc1,na.rm=TRUE))/30
text(1:length(lev),y=maxy,labels=n,cex=0.7,xpd=NA,col=col)
text(1:length(lev),y=maxy-3*unit,labels=smed,cex=0.7,xpd=NA,col=col)
text(1:length(lev),y=maxy-6*unit,labels=smea,cex=0.7,xpd=NA,col=col)
text(0.8,y=maxy,labels=expression(italic(N)~.bins==NA),pos=2,cex=0.8,xpd=NA,offset=1)
text(0.8,y=maxy-3*unit,labels="Median=",pos=2,cex=0.8,xpd=NA,offset=1)
text(0.8,y=maxy-6*unit,labels="Mean=",pos=2,cex=0.8,xpd=NA,offset=1)
axis(1,at=1:length(lev),labels=lev,tcl=-0.25)
#text(x=1:length(type),y=par('usr')[3],labels=type,cex=0.7,srt=45,xpd=NA,adj=c(1,1))
#text(x=median(1:length(type)),y=maxy*1.2,labels="Whole genome",cex=1.2,xpd=NA,font=2)
axis(2,las=2)
dev.off()



