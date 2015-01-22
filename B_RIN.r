#set the working directory

dir<- "H:\\B_RIN"
setwd(dir)
myfiles<-list.files()
which(substr(myfiles,69,79)=="Results.csv")
results<-myfiles[which(substr(myfiles,69,79)=="Results.csv")]
n <- which(substr(myfiles,69,79)=="Results.csv")
targets <- read.delim(myfiles[n[1]], sep=",",header=F,skip=14,as.is=T)
length(n)
targets2<- read.delim(myfiles[n[1]], sep=",",header=F,skip=0,as.is=T)
str1<-substr(targets2[1,2],1,67)
str1
str2<-"_S"
str3<-paste(str1,str2,sep="")
str3
chip_num<-length(n)
samp_num_chip<-length(grep(str3,myfiles))
samp_num<-length(which(substr(myfiles,69,74)=="Sample"))
samples1<- matrix(NA, ncol=chip_num, nrow=samp_num)

for(i in 1:chip_num)
{
#mod chip numb
targets2<- read.delim(myfiles[n[i]], sep=",",header=F,skip=0,as.is=T)
str1<-substr(targets2[1,2],1,67)
str2<-"_S"
str3<-paste(str1,str2,sep="")
samp_num_chip<-length(grep(str3,myfiles))

targets <- read.delim(results[i], sep=",",header=F,skip=14,as.is=T)
targets[1:5,]
which(targets[,1]=="Sample Name")
length(which(targets[,1]=="Sample Name"))
samples<-as.character(targets[which(targets[,1]=="Sample Name"),2])
tbl<-which(targets[,1]=="Overall Results:")
length(which(targets[,1]=="Overall Results:"))

tbl[1]
tbl[1]
colnms<-as.character(targets[(tbl[1]+1):(tbl[1]+4),1])

#mod sample numb
str1<-substr(targets2[1,2],1,67)
str4<-"_Sample"

for(j in 1:samp_num_chip)
{
iii<-j
end.st<-".csv"
str5<-paste(str1,str4,iii,end.st,sep="")
samples1[j,i]<-str5
}
}

samples1
grep("2100",samples1)
samples1<-samples1[grep("2100",samples1)]
samples1

write.table(samples1,file=paste(dir,"\\ba_lane.txt",sep=""), sep="\t",col.names=NA)


#RIN Calculator#######refined peakfinder####################

#read bioanalyzer data into a matrix called dta
#since the total RNA and mRNA assays run differently
#skip more lines for the total RNA assay

ramp <- colorRamp(c("darkmagenta","white"))
elec <- read.delim("ba_lane.txt", sep="\t",header=T,as.is=T)
qc.mat <- matrix(NA, ncol=2, nrow=nrow(elec))

dta <- matrix(NA, nrow=1060, ncol=nrow(elec))
for(i in 1:nrow(elec))
{
x <- read.csv(as.character(elec[i,2]), header=F, skip=18, nrows=1060)
dta[,i] <- x[,2]
}
time<- x[,1]

#plot ladder and define peak threshold in HD

par(mfrow=c(3,4))
#for(samples in 1:4)
for(samples in 1:nrow(elec))
{
#which(dta[,lad] > 5)
#draw some segments and store max fluor in max.peaks matrix
lad<-samples
ti1<-17
ti2<-65
p.row<-length(dta[which(time==ti1):which(time==ti1+1),lad])*(ti2-ti1+1)
peaks <- matrix(NA, nrow=p.row, ncol=1)
max.peaks <- matrix(NA, nrow=length(seq(ti1,ti2,.05)), ncol=1)

for (i in 1:length(seq(ti1,ti2,.05)))
{
ii<-seq(ti1,ti2,.05)[i]
x0<- ii
x0<-round(x0,digits=3)
y0<- ii
x<- ii+.05
x<-round(x,digits=3)
y<- ii
#segments(x0,y0,x,y)
max.peaks[i]<- max(dta[which(time==x0):which(time==x),lad])
}    

#human or mouse
cent<-max(max.peaks)*.3

#flatworm
#cent<-max(max.peaks)*.05

mat3 <- matrix("NA", ncol=length(colnms), nrow=length(samples1))
for(i in 1:length(results)) 
{
targets <- read.delim(results[i], sep=",",header=F,skip=14,as.is=T)
targets[1:5,]
which(targets[,1]=="Sample Name")
length(which(targets[,1]=="Sample Name"))
samples3<-as.character(targets[which(targets[,1]=="Sample Name"),2]) 
mat3[1:length(samples3),i]<-samples3
}
mat3
which(mat3!="NA")[i]

#for (n in 1:length(mat3[which(mat3!="NA")]))

length(mat3[which(mat3!="NA")])
plot(time,dta[,lad],type="h",col="bisque4",xlab="seconds",ylab="Fluorescence Units",xlim=c(17,75),ylim=c(0,(max(max.peaks)*1.2)),main=mat3[which(mat3!="NA")[samples]])
abline(h=cent,col="grey",lty="dashed")

#define the middle of the segments
max.avg<-seq(ti1+.025,ti2+.025,.05)

#plot the maximum value of each segment
max.peaks
lines(max.avg,max.peaks,type="n",col="bisque3",pch=20)















#lines(max.avg,min.peaks,type="b",col="blue")

#find time associated with peaks > 5 fu and fill seconds matrix with them

#12.5% of max fu incase chip is fu

p1<-which(max.avg==38)
p2<-which(max.avg==55)
sizes <- which(max.peaks>cent)
max.avg[sizes]
gt35<-which(max.avg[sizes]>35)
lt55<-which(max.avg[sizes]<55)
#as.numeric(duplicated(c(gt35,lt55)))
tfr<-which(duplicated(c(gt35,lt55))==TRUE)
c(gt35,lt55)[tfr]
max.avg[sizes][c(gt35,lt55)[tfr]]
max.peaks[sizes][c(gt35,lt55)[tfr]]

points(max.avg[sizes+1][c(gt35,lt55)[tfr]],max.peaks[sizes+1][c(gt35,lt55)[tfr]],col="springgreen3",pch=20)
points(max.avg[sizes-1][c(gt35,lt55)[tfr]],max.peaks[sizes-1][c(gt35,lt55)[tfr]],col="springgreen3",pch=20)
points(max.avg[sizes][c(gt35,lt55)[tfr]],max.peaks[sizes][c(gt35,lt55)[tfr]],col="tomato3",pch=20)










#colnms2<-c("max.avg[sizes]","max.peaks[sizes]","max.avg[sizes+1]","max.peaks[sizes+1]","max.avg[sizes-1]")
max.avg[sizes]
max.peaks[sizes]
max.avg[sizes+1]
max.peaks[sizes+1]
max.avg[sizes-1]
max.peaks[sizes-1]
area<-cbind(max.peaks[sizes-1][c(gt35,lt55)[tfr]],max.peaks[sizes][c(gt35,lt55)[tfr]],max.peaks[sizes+1][c(gt35,lt55)[tfr]])
time2<-cbind(max.avg[sizes-1][c(gt35,lt55)[tfr]],max.avg[sizes][c(gt35,lt55)[tfr]],max.avg[sizes+1][c(gt35,lt55)[tfr]])
p1<-area[which(area[,1]<cent),1][1]
t1<-time2[which(area[,1]<cent),1][1]
p2<-area[which(area[,1]<cent),1][2]
t2<-time2[which(area[,1]<cent),1][2]
#area[which(area[,2]<cent),2]
#area[which(area[,2]<cent),2]
p3<-area[which(area[,3]<cent),3][1]
t3<-time2[which(area[,3]<cent),3][1]
p4<-area[which(area[,3]<cent),3][2]
t4<-time2[which(area[,3]<cent),3][2]
#abline(v=area[which(area[,3]<cent),3])

#abline(h=11)
#abline(v=42.75)

v1<-t1
v2<-t2
v3<-t3
v4<-t4

abline(v=v1,col="grey",lty="dashed")
abline(v=v2,col="grey",lty="dashed")

abline(v=v3,col="grey",lty="dashed")
abline(v=v4,col="grey",lty="dashed")

#time[area[which(area[,1]<cent),1][1]==dta[,lad]]
#time[area[which(area[,3]<cent),3][1]==dta[,lad]]
#time[area[which(area[,1]<cent),1][2]==dta[,lad]]
#time[area[which(area[,3]<cent),3][2]==dta[,lad]]
#abline(v=time[area[which(area[,1]<cent),1][1]==dta[,lad]])
#abline(v=time[area[which(area[,3]<cent),3][1]==dta[,lad]])
#abline(v=time[area[which(area[,1]<cent),1][2]==dta[,lad]])
#abline(v=time[area[which(area[,3]<cent),3][2]==dta[,lad]])
#abline(h=area[which(area[,3]<cent),3][1])
#abline(h=area[which(area[,3]<cent),3][2]+3)

s18a<-round(t1,digits=1)
s18b<-round(t3,digits=1)
s28a<-round(t2,digits=1)
s28b<-round(t4,digits=1)

if(is.na(s18a)){s18a<-30}
if(is.na(s18b)){s18b<-35.05}
if(is.na(s28a)){s28a<-30}
if(is.na(s28b)){s28b<-30.05}

eta<-which(time==s18a)
etb<-which(time==s18b)

tea<-which(time==s28a)
teb<-which(time==s28b)

#dta[eta:etb,lad]
#dta[tea:teb,lad]
s18<-sum(dta[eta:etb,lad])
s28<-sum(dta[tea:teb,lad])
s28/s18
qc.mat[samples,2]<-  round((-1*exp((s28/s18)*-1.5)+1)*10,digits=2)
qc.mat[samples,1]<-  mat3[which(mat3!="NA")[samples]]
text(60,cent+1*cent,c(format(qc.mat[samples,2],digits=2),"\n\nB-RIN"))
Sys.sleep(1)
}
qc.mat
###end of ratio calc
###
###
#barplot(qc.mat,beside=T,ylim=c(0,11))
write.table(qc.mat,file=paste(dir,"\\B-RIN_dat.txt",sep=""), sep="\t",col.names=NA)



                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                






myfiles<-list.files()
which(substr(myfiles[1:100],69,79)=="Results.csv")
results<-myfiles[which(substr(myfiles[1:100],69,79)=="Results.csv")]
n <- which(substr(myfiles[1:100],69,79)=="Results.csv")
targets <- read.delim(myfiles[n[1]], sep=",",header=F,skip=14,as.is=T)
targets[1:5,]
tbl<-which(targets[,1]=="Overall Results:")
which(targets[,1]=="Sample Name")
length(which(targets[,1]=="Sample Name"))
samples<-as.character(targets[which(targets[,1]=="Sample Name"),2])
#read in the targets file                               
colnms<-as.character(targets[(tbl[1]+1):(tbl[1]+4),1])



mat1 <- matrix(NA, ncol=length(colnms), nrow=nrow(qc.mat))
colnames(mat1)<-colnms
rownames(mat1)<-qc.mat[,1]
#read.delim("samples02.txt", sep="\t",header=T,as.is=T)
for(i in 1:length(results))
{
targets <- read.delim(results[i], sep=",",header=F,skip=14,as.is=T)
targets[1:5,]
which(targets[,1]=="Sample Name")
length(which(targets[,1]=="Sample Name"))
samples<-as.character(targets[which(targets[,1]=="Sample Name"),2])
tbl<-which(targets[,1]=="Overall Results:")
length(which(targets[,1]=="Overall Results:"))
tbl[1]
tbl[1]
#colnms<-as.character(targets[(tbl[1]+1):(tbl[1]+4),1])





for (j in 1:length(samples))
{
mat1[((12*(i-1))+j),]<-as.matrix(targets[(tbl[j]+1):(tbl[j]+4),2])
}
}
mat1
qc.mat





cbind(mat1,qc.mat)
colnames(qc.mat)<-c("samples_name","B-RIN")
write.table(mat1,file=paste(dir,"\\conc.txt",sep=""), sep="\t",col.names=NA)
write.table(cbind(mat1,qc.mat),file=paste(dir,"\\conc.txt",sep=""), sep="\t",col.names=NA)
##