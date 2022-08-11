dat<-read.table("add.panda.txt",stringsAsFactors=F,header=TRUE)
dat<-dat[order(dat$ord),]

c1<-1200000
c2<-1200000


cols<-rep(c('#3CB371','#EF9D55','#D8272A','#DF5693','black'),times=c(5,2,3,5,1))
pdf("genome.compare.alt.pdf",w=7,h=4)
layout(matrix(2:1,nr=2),w=c(6,0),h=c(2.5,2.5))
par(mar=c(8,5,0.05,5),xaxs='i',yaxs='i',tcl=-0.25)
plot(1:10,axes=F,type="n",xlim=c(1,length(dat$Specy))+c(-1,1),ylim=c(0,c1),xlab="",ylab="")

wid<-0.35
for (i in seq(dat$Specy)){
	if(dat$Contig[i]<=c1){
			rect(i-wid,0,i+wid,dat$Contig[i],lwd=1,col=cols[i],border='black')
			text(i,dat$Contig[i],pos=3,labels=sprintf("%.2f",dat$Contig[i]/1000000),xpd=NA,cex=0.7,col=cols[i])
	}else{
			rect(i-wid,0,i+wid,c1,lwd=1,col=cols[i],border='black')
	
	}
}
axis(2,las=2,at=pretty(par('usr')[3:4],n=4),labels=pretty(par('usr')[3:4],n=4)/1000000)
axis(1,at=1:length(dat$Specy),labels=NA)
box(bty='l')
xut<-diff(par('usr')[1:2])/100
yut<-diff(par('usr')[3:4])/50
#segments(x0=0.5-xut,y0=c1-yut,x1=0.5+xut,y1=c1+yut,xpd=NA,lwd=1)
#text(0.51,c1-yut,srt=45,labels='-',adj=1,xpd=NA)
text(x=1:length(dat$Specy)+0,y=-yut*10,labels=dat$Specy,srt=45,adj=1,xpd=NA,cex=0.8,col=cols)


#st<--3000000
st<--2700000
ut<-700000
ce<-0.75
text(1:length(dat$Specy),y=st,labels=sprintf("%.2f",dat$genome/1000000000),cex=ce,xpd=NA,col=cols)
text(1:length(dat$Specy),y=st-1*ut,labels=dat$chr,cex=ce,xpd=NA,col=cols)
#text(1:length(dat$Specy),y=st-2*ut,labels=dat$version,cex=0.5,xpd=NA,col=cols,srt=18)
#text(1:length(dat$Specy),y=st-3*ut,labels=dat$Protein_coding,cex=ce,xpd=NA,col=cols)
#text(1:length(dat$Specy),y=st-4*ut,labels=dat$lncRNA,cex=ce,xpd=NA,col=cols)
#text(1:length(dat$Specy),y=st-5*ut,labels=dat$miRNA,cex=ce,xpd=NA,col=cols)	

text(1,y=st,labels="Genome size\n(Gb)       ",pos=2,cex=ce,xpd=NA,offset=1)
text(1,y=st-1*ut,labels=expression(italic(N.)~Chr),pos=2,cex=ce,xpd=NA,offset=1)
#text(1,y=st-2*ut,labels="Version",pos=2,cex=ce,xpd=NA,offset=1)
#text(1,y=st-3*ut,labels=expression(italic(N.)~Protein~~coding),pos=2,cex=ce,xpd=NA,offset=1)
#text(1,y=st-4*ut,labels=expression(italic(N.)~lncRNA),pos=2,cex=ce,xpd=NA,offset=1)
#text(1,y=st-5*ut,labels=expression(italic(N.)~miRNA),pos=2,cex=ce,xpd=NA,offset=1)


#####
par(mar=c(0.05,5,2,5),xaxs='i',yaxs='i',tcl=-0.25)
plot(1:10,axes=F,type="n",xlim=c(1,length(dat$Specy))+c(-1,1),ylim=c(c2,max(dat$Contig)),ylab="")

for (i in seq(dat$Specy)){
    if(dat$Contig[i]>=c2){
            rect(i-wid,c2,i+wid,dat$Contig[i],lwd=1,col=cols[i],border='black',xpd=NA)
			text(i,dat$Contig[i],pos=3,labels=sprintf("%.2f",dat$Contig[i]/1000000),xpd=NA,cex=0.7,col=cols[i])
    }
}
axis(2,las=2,at=pretty(par('usr')[3:4],n=5),labels=pretty(par('usr')[3:4],n=5)/1000000)
box(bty='l')
xut<-diff(par('usr')[1:2])/100
yut<-diff(par('usr')[3:4])/50
#segments(x0=0.5-xut,y0=c2-yut,x1=0.5+xut,y1=c2+yut,xpd=NA,lwd=1)
#text(0.51,c2,srt=45,labels='-',adj=1,xpd=NA)
rect(-0.1,c2-1*yut,length(dat$Specy)+1,c2+yut,border=NA,col='white',xpd=NA)
text(0.01,c2-yut,srt=25,labels='-',adj=1,xpd=NA)
text(0.01,c2+yut*1.1,srt=25,labels='-',adj=1,xpd=NA)


par(fig = c(0,1,0,1), new=TRUE,mar=c(8,5,2,5),xaxs='i',yaxs='i')
plot(1:10,axes=F,type="n",xlim=c(1,length(dat$Specy))+c(-1,1),ylim=c(0,max(dat$scaffold)),xlab="",ylab='')
points(1:length(dat$Specy),y=dat$scaffold,pch=19,col='black',xpd=NA,cex=0.6)

lines(spline(x=1:length(dat$Specy),y=dat$scaffold,n=1000),col=adjustcolor('black',alpha.f=0.5),xpd=NA)
axis(4,las=2,at=pretty(par('usr')[3:4],n=5),labels=pretty(par('usr')[3:4],n=5)/1000000)
mtext(4,text='Scaffold N50 (Mb)',line=2.5)
mtext(2,text='Contig N50 (Mb)',line=2.5)
dev.off()


