######
library(viridis)
options(scipen=200)
d1<-data.frame(reso=c(1000,2000,5000,10000,20000,40000,100000,500000,1000000),stringsAsFactors=F)
d1$cols<-adjustcolor(viridis(length(d1$reso)),alpha.f=0.8)
#d2<-data.frame(tp=c('','.cis','.trans'),type=c(1,2,3),leg=c('Cis+Trans','Cis','Trans'),stringsAsFactors=F)
d2<-data.frame(tp=c('.cis'),type=c(1),leg=c('Cis'),stringsAsFactors=F)


dat<-data.frame(stringsAsFactors=F)
for (i in d1$reso){
	for (j in d2$tp){
		fil<-paste("JJ",j,".",i,".bed.stat.reso.interval.xls",sep="")
		idat<-read.table(fil,stringsAsFactors=F,sep="\t",header=TRUE)
		idat$reso<-i
		idat$tp<-j
		dat<-rbind(dat,idat)
}
}

dat<-merge(dat,d1,by='reso')
dat<-merge(dat,d2,by='tp')

pdf("resolution.pdf",w=4,h=4)
par(xaxs='i',yaxs='i',mgp=c(2,0.5,0),lend=2,tcl=-0.25)
plot(1:10,type='n',axes=F,xlim=c(1,6),ylim=c(-0.01,1.01),xlab=expression(Number~~of~~Hi-C~~contacts~(10^{axis-x})),ylab="Cumulative percentage\nof bins (>1000 contacts)",main="Resolution evaluation by only Cis contacts",cex.main=0.9)
#grid(lwd=0.5,col='gray')
for (i in d1$reso){
    for (j in d2$tp){
        idat<-dat[dat$tp==j&dat$reso==i,]
		lines(seq(0,6,by=0.1)[1:60],idat$cpct,lwd=1,lty=idat$type,col=idat$cols)
}
}
axis(1)
axis(2,las=2)
box(bty='o')
abline(v=3,col='black',lty=2,lwd=0.5)

legend("topleft",legend=d1$reso,col=d1$cols,lwd=4,nc=1,bty='n',xpd=NA,cex=0.7)
dev.off()
