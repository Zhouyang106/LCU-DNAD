#OPLS-DA
library(ropls)
library(ggplot2)
data <- read.table("Input.csv", header =T, sep = ",",row.name=1)
group<- read.table("group.txt", header = TRUE, sep = ",", quote="\"",  stringsAsFactors = TRUE)
data1 <- t(data)
class(group)
group <- factor(group)
data.plsda <- opls(x=data1,y=group1)
data.plsda
sacurine.oplsda<-opls(x=data1,group1,predI=1,orthoI=NA)
#Functional Enrichment Analysis
colnames(pathway) <- c("Description","pvalue","zuobiao")
pathway$Description <- factor(pathway$Description, levels=unique(pathway$Description)) ##注释成为因子
ggplot(pathway,aes(zuobiao,Description))+geom_point(aes(size=-log10(pvalue), color=pvalue,inherit.aes = T))+
  scale_color_gradient2(low = "#FC5C7D",mid = "#FC5C7D", high = "#396afc",midpoint =0.00001)+theme(plot.title=element_text(hjust=5))+ggtitle("pathway") 
#Metabolite scatter plots
ggplot(data, aes(x = rtmed, y = mzmed, color = Pvalue, size = Fold_change)) +
  geom_point() +
  scale_size(range=c(1,12))+
  scale_color_gradient(low = "#EAEAEA", high = "#6A82FB") +
  labs(x = "rtmed", y = "mzmed") +
  theme_minimal()
#Spearman’s correlation analysis
library(corrplot)
library(psych) 
spearman1 <- corr.test(data, method = 'spearman') 
spearman1[[11]]
p <- spearman1$p
r <- spearman1$r
col4 <- colorRampPalette(c("#a8c0ff", "#FBD786"))
corrplot(spearman1$r, method = 'square', type = 'lower', 
         p.mat = spearman1$p, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), # P值小于0.001,0.01,0.05显示*，*越多越显著
         pch.cex = 1.2, pch.col = "black", 
         diag = FALSE, col = col4(20), tl.col = '#000000',
         tl.srt = 90) # 添加此参数使横坐标标签旋转45度
dev.off()
