library(viridis)
library(RColorBrewer)
library(scales)
args<-commandArgs(TRUE)
FAI<-args[1]
FILE<-args[2]
CEN<-args[3]
OUP<-args[4]

fai<-read.table(FAI,stringsAsFactors=F,sep="\t")
fai$V1<-as.numeric(gsub("X","21",fai$V1))
file<-read.table(FILE,stringsAsFactors=F,sep="\t",fill=TRUE)
file$V1<-as.numeric(gsub("X","21",file$V1))
file<-file[file$V4>=0.3,]


cen<-read.table(CEN,stringsAsFactors=F,sep="\t")
cen$V1<-as.numeric(gsub("X","21",cen$V1))


########round_loci
wid<-0.3
r<-8000000
z<-seq(0,2*pi,length=1000)
y<-sin(z)*r
x<-cos(z)*wid

###########plot
pdf(OUP,w=8,h=3)
par(xaxs='i',yaxs='i',lend=2,mar=c(1.5,4,1,1))
plot(1:10,axes=F,ylab="Length (Mb)",xlab="Chromosome",ylim=c(-r,max(fai$V2)),xlim=rev(c(21.5,0)),type='n',xpd=NA)

################
#pal<-viridis(5)
#pal<-c('grey','#EBB23E')
pal<-c('white','cyan3')
#pal<-brewer.pal(9,"Reds")[1:9]

#cols<-col_numeric(pal,1:5)(1:5)
#leg<-data.frame(V5=c('L2','L1','H1','H2','H3'),cols=cols,stringsAsFactors=F)
#file<-file[file$V1%in%1:18,]
#file<-file[file$V5!="gap",]
#file<-merge(file,leg,by='V5')
#print(head(file))
#rect(as.numeric(file$V1)-wid,file$V2,as.numeric(file$V1)+wid,file$V3,border=NA,col=file$cols)

file$cols<-col_numeric(pal,file$V4)(file$V4)
rect(as.numeric(file$V1)-wid,file$V2,as.numeric(file$V1)+wid,file$V3,border=NA,col=file$cols)

###########
#legend("top",legend=leg$V5,lwd=6,nc=3,col=leg$cols,bty='n',cex=1.5)
########rect

for (chr in fai$V1){
        clen<-fai$V2[fai$V1==chr]
        ratio<-cen$V2[cen$V1==chr]
        if (ratio>0){
                cpos<-as.integer(clen*ratio)
		rect(chr-wid,cpos-r,chr+wid,cpos+r,border='white',col='white')
                cfai<-fai[fai$V1==chr,]
                points(x=chr,y=cpos,pch=21,col='grey',bg='black',lwd=0.5,cex=0.7)
                px<-c(chr+wid,chr+x[1:500],chr-wid,chr+x[501:1000],chr+wid)
		ncex<-1.2
                py1<-c(cfai$V2,cfai$V2+y[1:500],cpos+r*ncex,cpos+r*ncex+y[501:1000],cpos+r*ncex)
                py2<-c(cpos-r*ncex,cpos-r*ncex+y[1:500],0,0+y[501:1000],0)
                polygon(px,py1,border='black',lwd=0.5,col=adjustcolor(pal[9],alpha.f=0.1),xpd=NA)
                polygon(px,py2,border='black',lwd=0.5,col=adjustcolor(pal[9],alpha.f=0.1),xpd=NA)
        }else{
                cfai<-fai[fai$V1==chr,]
                px<-c(chr+wid,chr+x[1:500],chr-wid,chr+x[501:1000],chr+wid)
                py<-c(cfai$V2,cfai$V2+y[1:500],0,0+y[501:1000],0)
                polygon(px,py,border='black',lwd=0.5,col=adjustcolor(pal[9],alpha.f=0.1),xpd=NA)
        }
}
text(x=1:21,y=-r,labels=c(1:20,"X"),pos=1,xpd=NA)
axis(2,at=pretty(par('usr')[3:4],n=5),labels=pretty(par('usr')[3:4],n=5)/1000000,las=2)

par(fig=c(0.7,0.9,0.8,0.85),new=TRUE,mar=c(0,0,0,0),mgp=c(1.5,0.5,0),tcl=-0.2)
xseq<-seq(from=range(file$V4)[1],to=range(file$V4)[2],len=200)
image(x=xseq,y=1:2,z=t(t(xseq)),col=col_numeric(pal,xseq)(xseq),axes=F)
axis(1,xpd=NA,cex.axis=0.9,gap.axis=0)
mtext(3,text="GC percentage",cex=1,xpd=NA)
box(bty='o')

dev.off()
