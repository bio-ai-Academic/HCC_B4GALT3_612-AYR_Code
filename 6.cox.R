#install.packages("survival")

library(survival)             #引用包
inputFile="expTime.txt"       #输入文件
pFilter=0.05                  #显著性过滤条件
setwd("C:\\Users\\lexb4\\Desktop\\TMEimmune\\22.cox")       #设置工作目录

#读取输入文件
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
rt$futime=rt$futime/365           #以年为单位，除以365

outTab=data.frame()
for(gene in colnames(rt[,3:ncol(rt)])){
	#单因素cox分析
	cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	
	#KM检验
	group=ifelse(rt[,gene]<=median(rt[,gene]),"Low","High")
	diff=survdiff(Surv(futime, fustat) ~ group,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	
	#保存满足条件的基因
	if((pValue<pFilter) & (coxP<pFilter)){
		outTab=rbind(outTab,
		         cbind(gene=gene,
		               KM=pValue,
		               HR=coxSummary$conf.int[,"exp(coef)"],
		               HR.95L=coxSummary$conf.int[,"lower .95"],
		               HR.95H=coxSummary$conf.int[,"upper .95"],
				       pvalue=coxP) )		
	}
}
#输出基因和p值表格文件
write.table(outTab,file="cox.result.txt",sep="\t",row.names=F,quote=F)

######绘制森林图######
#读取输入文件
rt <- read.table("cox.result.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#输出图形
pdf(file="forest.pdf", width = 7,height = 6)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

#绘制森林图左边的基因信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
