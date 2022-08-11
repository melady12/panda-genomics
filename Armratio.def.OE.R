library(plyr)
library(reshape2)
args<-commandArgs(TRUE)
inp<-args[1]
reso<-as.numeric(args[2])
oup<-args[3]
chr<-args[4]

dat<-read.table(gzfile(inp),stringsAsFactors=F,sep="\t")
dat<-as.matrix(dat)
nc<-dim(dat)[1]
nn<-as.integer(nc/10)
##########def_outlier
filt<-function(x){
        rx<-range(boxplot(x,plot=FALSE)$stats)
        return(x[(x>=rx[1])&(x<=rx[2])])
}


#########def_function
deal<-function(x){
        nr<-dim(x)[1]
        rownames(x)<-1:nr
        colnames(x)<-1:nr
        mdat<-melt(x)
        mdat$dis<-abs(mdat$Var1-mdat$Var2)
        odat<-ddply(mdat,.(dis),summarize,mea=mean(filt(value)),med=median(filt(value)))
        return(odat)
}


pdat<-data.frame(stringsAsFactors=F)
for (i in seq(nn+1,nc-nn,by=1)){
        #imat<-log(dat[i:(i+nn),i:(i+nn)]+1,10)
        imat<-dat[i:(i+nn),i:(i+nn)]
        tdat<-deal(imat)
        idat<-data.frame(pos=i*reso,mea=mean(filt(tdat$mea)),med=median(filt(tdat$med)),sd=sd(filt(tdat$mea)),stringsAsFactors=F)
        pdat<-rbind(pdat,idat)
}

########Armratio
ii<-pdat[order(pdat$mea),c('mea','pos')]
mdis<-diff(range(ii$pos[1:4]))
if(mdis<=3000000){
	arm<-mean(ii$pos[1:4])/(nc*reso)
	arm<-sprintf("%.3f",arm)
}else{
	arm<-"NA"
}
pdf(oup,w=4,h=4)
par(lend=2,tcl=-0.25,mgp=c(2,0.5,0))
plot(1:10,type='n',axes=F,xlim=range(pdat$pos),ylim=range(c(pdat$mea-pdat$sd,rev(pdat$mea+pdat$sd))),main=paste("Centromere ratio of chromosome ",chr,": ",arm,sep=""),xlab="Length of chromosome (Mb)",ylab="Hi-C contacts\n(Observed/Expected)",cex.main=0.9)
grid(lwd=0.5)
polygon(x=c(pdat$pos,rev(pdat$pos)),y=c(pdat$mea-pdat$sd,rev(pdat$mea+pdat$sd)),col=adjustcolor("gray",alpha.f=0.5),border=NA,xpd=NA)
lines(pdat$pos,pdat$mea,lwd=1)
axis(2,las=2,gap.axis=-1)
axis(1,at=pretty(par('usr')[1:2]),labels=pretty(par('usr')[1:2])/1000000,gap.axis=-1)
box(bty='o')
dev.off()
