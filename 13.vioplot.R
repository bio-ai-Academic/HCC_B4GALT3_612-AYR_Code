#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("vioplot")


#引用包
library(limma)
library(vioplot)

immuneFile="CIBERSORT-Results.txt"         #免疫浸润结果文件
expFile="symbol.txt"                       #表达输入文件
gene="CXCR4"                               #基因名称
pFilter=0.05                               #CIBERSORT结果过滤条件
setwd("C:\\Users\\lexb4\\Desktop\\TMEimmune\\33.vioplot")       #设置工作目录

#读取免疫结果文件，并对数据进行整理
immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#去除正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]

#按基因表达对样品分组
lowName=colnames(data)[data[gene,]<=median(data[gene,])]
highName=colnames(data)[data[gene,]>median(data[gene,])]

#提取高低表达组免疫细胞含量
lowImm=intersect(row.names(immune),lowName)
highImm=intersect(row.names(immune),highName)
rt=rbind(immune[lowImm,],immune[highImm,])
lowNum=length(lowImm)
highNum=length(highImm)

#输出小提琴图
outTab=data.frame()
pdf("vioplot.pdf",height=8,width=13)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，低表达用绿色表示，高表达用红色表示
for(i in 1:ncol(rt)){
	  if(sd(rt[1:lowNum,i])==0){
	    rt[1,i]=0.001
	  }
	  if(sd(rt[(lowNum+1):(lowNum+highNum),i])==0){
	    rt[(lowNum+1),i]=0.001
	  }
	  lowData=rt[1:lowNum,i]
	  highData=rt[(lowNum+1):(lowNum+highNum),i]
	  vioplot(lowData,at=3*(i-1),lty=1,add = T,col = 'green')
	  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
	  wilcoxTest=wilcox.test(lowData,highData)
	  p=wilcoxTest$p.value
	  if(p<pFilter){
	      cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
		  outTab=rbind(outTab,cellPvalue)
	  }
	  mx=max(c(lowData,highData))
	  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
	  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

#输出免疫细胞和p值表格文件
write.table(outTab,file="diff.result.txt",sep="\t",row.names=F,quote=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
