#install.packages("survival")
#install.packages("survminer")

#引用包
library(survival)
library(survminer)

inputFile="expTime.txt"      #输入文件名称
gene="CXCR4"                 #基因名字
setwd("C:\\Users\\lexb4\\Desktop\\TMEimmune\\26.geneSur")    #设置工作目录

#读取输入文件
rt=read.table(inputFile,header=T,sep="\t",check.names=F)
rt$futime=rt$futime/365        #除以365，单位改成年

#根据gene中位值对样品分组
a=ifelse(rt[,gene]<=median(rt[,gene]),"Low","High")
#高低表达组生存差异
diff=survdiff(Surv(futime, fustat) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
#计算每个时间节点病人数目
fit=survfit(Surv(futime, fustat) ~ a, data = rt)
#绘制生存曲线
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
titleName=gene
surPlot=ggsurvplot(fit, 
				data=rt,
				conf.int=TRUE,
				pval=pValue,
				pval.size=6,
				risk.table=T,
				#ncensor.plot = TRUE,
				legend.labs=c("high","low"),
				legend.title=titleName,
				xlab="Time(years)",
				break.time.by = 1,
				risk.table.title="",
				palette=c("red", "blue"),
				risk.table.height=.25)          
pdf(file=paste0("sur.",gene,".pdf"), width = 6.5, height = 5.5,onefile = FALSE)
print(surPlot)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
