library(plyr)
library(gtools)
library(vioplot)
library(viridis)
library(scales)



d1<-read.table("map_cov.release_GPv1.xls",stringsAsFactors=F,header=TRUE,sep="\t")
d2<-read.table("map_cov.release_Qingling.xls",stringsAsFactors=F,header=TRUE,sep="\t")
d3<-read.table("map_cov.release_Sichuan.xls",stringsAsFactors=F,header=TRUE,sep="\t")
d4<-read.table("map_cov.release_ASMv2.xls",stringsAsFactors=F,header=TRUE,sep="\t")
d5<-read.table("map_cov.release_AME.xls",stringsAsFactors=F,header=TRUE,sep="\t")

AAA<-function(x,y){
	nn<-names(x)[2:length(names(x))]
	names(x)[2:length(names(x))]<-paste(nn,y,sep="")
	return(x)
}
d1<-AAA(d1,".A")
d2<-AAA(d2,".B")
d3<-AAA(d3,".C")
d4<-AAA(d4,".D")
d5<-AAA(d5,".E")

dat<-merge(d1,d2,by='Sample')
dat<-merge(dat,d3,by='Sample')
dat<-merge(dat,d4,by='Sample')
dat<-merge(dat,d5,by='Sample')

deal<-function(x){as.numeric(gsub("%","",x))}
dat[,2:dim(dat)[2]]<-apply(dat[,2:dim(dat)[2]],2,deal)
dat<-dat[order(dat$unique.and.mismatch..8.rate.over.all.A),]
#dat<-dat[order(dat$Sample),]
dat$xx<-1:length(dat$Sample)


lev<-c('.A','.B','.C','.D','.E')
col<-adjustcolor(col_numeric(c('#a9252d','#2969a6'),range(1:length(lev)))(1:length(lev)),alpha.f=0.8)
n<-c()
med<-c()
mea<-c()

pdf(paste("Mapping.ratio.vioplot.pdf",sep=""),w=3,h=4)
par(xaxs='i',yaxs='i',lend=2,mgp=c(2,0.6,0),tcl=-0.3,mar=c(5,4,5,2),bty='l')
plot(1:10,type='n',axes=FALSE,xlab='',ylab="Mapping ratio (%)",xlim=c(0.3,length(lev)+0.5),ylim=c(92.5,99))

for (j in seq(lev)){
	pat<-paste("mapping_ratio",lev[j],sep="")
	idat<-dat[,pat]
	vec<-as.numeric(as.vector(idat))
	vioplot(vec,at=j,col=col[j],border=col[j],add=TRUE)
	n<-c(n,length(vec))
	med<-c(med,median(vec))
	mea<-c(mea,mean(vec))
}
smed<-sprintf("%.2f",med)
smea<-sprintf("%.2f",mea)
#maxy<-max(dd$pc1,na.rm=TRUE)*1.7
maxy<-101
#unit<-diff(range(idat$pc1,na.rm=TRUE))/30
unit<-0.2
text(1:length(lev),y=maxy,labels=n,cex=0.7,xpd=NA,col=col)
text(1:length(lev),y=maxy-3*unit,labels=smed,cex=0.7,xpd=NA,col=col)
text(1:length(lev),y=maxy-6*unit,labels=smea,cex=0.7,xpd=NA,col=col)
text(0.8,y=maxy,labels=expression(italic(N)~.==NA),pos=2,cex=0.8,xpd=NA,offset=1)
text(0.8,y=maxy-3*unit,labels="Median=",pos=2,cex=0.8,xpd=NA,offset=1)
text(0.8,y=maxy-6*unit,labels="Mean=",pos=2,cex=0.8,xpd=NA,offset=1)
axis(1,at=1:length(lev),labels=NA,tcl=-0.25)

text(x=1:length(lev),y=par('usr')[3]*0.999,labels=c('GPv1','Qingling','Sichuan','ASM200744v2','AilMel_1.0'),cex=0.9,srt=45,xpd=NA,adj=c(1,1),col=col)

#text(x=1:length(type),y=par('usr')[3],labels=type,cex=0.7,srt=45,xpd=NA,adj=c(1,1))
#text(x=median(1:length(type)),y=maxy*1.2,labels="Whole genome",cex=1.2,xpd=NA,font=2)
axis(2,las=2)
dev.off()


n<-c()
med<-c()
mea<-c()

pdf(paste("mismatch_ratio.vioplot.pdf",sep=""),w=3,h=4)
par(xaxs='i',yaxs='i',lend=2,mgp=c(2,0.6,0),tcl=-0.3,mar=c(5,4,5,2),bty='l')
plot(1:10,type='n',axes=FALSE,xlab='',ylab="Mismatch ratio (%)",xlim=c(0.3,length(lev)+0.5),ylim=c(0.5,2))

for (j in seq(lev)){
    pat<-paste("mismatch_ratio",lev[j],sep="")
    idat<-dat[,pat]
    vec<-as.numeric(as.vector(idat))
    vioplot(vec,at=j,col=col[j],border=col[j],add=TRUE)
    n<-c(n,length(vec))
    med<-c(med,median(vec))
    mea<-c(mea,mean(vec))
}
smed<-sprintf("%.2f",med)
smea<-sprintf("%.2f",mea)
#maxy<-max(dd$pc1,na.rm=TRUE)*1.7
maxy<-2.5
#unit<-diff(range(idat$pc1,na.rm=TRUE))/30
unit<-0.05
text(1:length(lev),y=maxy,labels=n,cex=0.7,xpd=NA,col=col)
text(1:length(lev),y=maxy-3*unit,labels=smed,cex=0.7,xpd=NA,col=col)
text(1:length(lev),y=maxy-6*unit,labels=smea,cex=0.7,xpd=NA,col=col)
text(0.8,y=maxy,labels=expression(italic(N)~.==NA),pos=2,cex=0.8,xpd=NA,offset=1)
text(0.8,y=maxy-3*unit,labels="Median=",pos=2,cex=0.8,xpd=NA,offset=1)
text(0.8,y=maxy-6*unit,labels="Mean=",pos=2,cex=0.8,xpd=NA,offset=1)
axis(1,at=1:length(lev),labels=NA,tcl=-0.25)

text(x=1:length(lev),y=par('usr')[3]*0.95,labels=c('GPv1','Qingling','Sichuan','ASM200744v2','AilMel_1.0'),cex=0.9,srt=45,xpd=NA,adj=c(1,1),col=col)

#text(x=1:length(type),y=par('usr')[3],labels=type,cex=0.7,srt=45,xpd=NA,adj=c(1,1))
#text(x=median(1:length(type)),y=maxy*1.2,labels="Whole genome",cex=1.2,xpd=NA,font=2)
axis(2,las=2)
dev.off()







