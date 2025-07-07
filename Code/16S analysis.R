#Bray-Curtis principal coordinate analysis
library(vegan)
library(ggplot2)
library(plyr)
otu <- read.delim('data.csv', row.names = 1, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))
group <- read.csv(file = "group.csv",header = T,sep = ",")
distance <- vegdist(otu, method = 'bray')
pcoa <- cmdscale(distance, k = (nrow(otu) - 1), eig = TRUE)
ordiplot(scores(pcoa)[ ,c(1, 2)], type = 't')
pcoa$eig
point <- data.frame(pcoa$point)
write.csv(point, 'pcoa.sample.csv')
species <- wascores(pcoa$points[,1:2], otu)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample<- data.frame({pcoa$point})[1:2]
sample$names <- rownames(sample)
names(sample)[1:2] <- c('PCoA1', 'PCoA2')
sample$treat <- factor(sample$treat, levels = c('COB', 'OB'))
group_border <- ddply(sample, 'treat', function(df) df[chull(df[[2]], df[[3]]), ])
pcoa_plot <- ggplot(sample, aes(PCoA1, PCoA2, group = treat)) +
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.key = element_rect(fill = 'transparent')) + #去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_polygon(data = group_border, aes(fill =treat )) + 
  geom_point(aes(color = treat, shape = treat), size =1.5, alpha =0.8) + 
  scale_shape_manual(values = c(15,16,17)) + 
  scale_color_manual(values = c('#6699CC', '#FF3399', '#3300CC')) + 
  scale_fill_manual(values = c('#bce6eb', '#FF99FF', '#3366CC')) + 
  guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2), color = guide_legend(order = 3)) + #设置图例展示顺序
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%')) 
pcoa_plot
#Barplot
par(xpd = TRUE, mar = par()$mar + c(1, 3, 1, 16))
barplot(as.matrix(100 *phylum))
barplot(as.matrix(100* phylum),
        col = c('#3366CC', '#89c9b8', '#ffd571', '#bbd196', '#c3aed6', '#fbe2e5', '#2bb2bb','#ff847c', '#faf0af', '#ff8ba7', 'gray'),
        legend = rownames(phylum), 
        cex.axis = 2, cex.names = 2, ylim = c(0, 100), las = 1, width = 0.9, space = 0.6, beside = FALSE,
        args.legend = list(x = 4,y=100, 
                           bty = 'n', 
                           inset = -0.18, 
                           cex = 0.5, 
                           y.intersp = 0.7, 
                           x.intersp = 0.8, 
                           text.width = 0.05
        ))
mtext('Relative Abundance(%)', cex = 2, side = 2, line = 4)
#Heatmap
data <- read.delim('heatmap.csv', row.names = 1, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)
exprSet <- data
str(exprSet)
exprSet[, 11] <- as.numeric(exprSet[, 11]) 
qx <- as.numeric(quantile(exprSet, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
LogC
exprSet <- log2(exprSet[,]+1)
pheatmap(exprSet)
my_col<-colorRampPalette(c(c("#1a4fb3","#396afc","white","#D64161","#D7263D")))(100)
pheatmap(exprSet,scale = "row",cluster_cols=F,cluster_rows=T,show_rownames = F, color = my_col,cellwidth =5,show_colnames = T,
         cellheight =0.1, legend = T, fontsize = 6, border_color = 8 )
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
         tl.srt = 90) 
dev.off()
