cepath="C:/Users/stu/Desktop/exam"
if(getwd()!=cepath)setwd(cepath)
# -------------1-----------------------
#temp=read.csv2(file = "metadata.csv",header=TRUE,seq="\t")
pattern_str_C="Control Neurones"
pattern_str_D="Diabetic Neurones"
library(stringr)
judge_Cneu=!is.na(str_extract(pattern = pattern_str_C,metadata$Sample_organism))
judge_Dneu=!is.na(str_extract(pattern = pattern_str_D,metadata$Sample_organism))
C_sample_names=metadata$Sample_status[judge_Cneu]
D_sample_names=metadata$Sample_status[judge_Dneu]
C_sample=genematrix[C_sample_names]
D_sample=genematrix[D_sample_names]
T1_result=rbind(D_sample,C_sample)
rownames(T1_result)=genematrix$ENSEMBL
save(T1_result,file="C_D.csv")
# -------------2-----------------------
C_avg=apply(C_sample,1,mean)
D_avg=apply(D_sample,1,mean)
Fold_change=D_avg/C_avg
p_list=c()
for(i in 1:nrow(C_sample)){
  p_list=append(p_list,t.test(C_sample[i,],D_sample[i,],var.equal=TRUE)$p.value)
}
T2_result=data.frame('Control_avg'=C_avg,'Diabetes_avg'=D_avg,'Fold_change'=Fold_change,'p_list'=p_list)
rownames(T2_result)=genematrix$ENSEMBL
save(T2_result,file="T2_result.csv")
# -------------3-----------------------
symbol_ID=bitr(genematrix$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",
               OrgDb = org.Hs.eg.db,drop=FALSE)
T2_result$SYMBOL=symbol_ID$SYMBOL
# -------------4-----------------------
T2_result$Fold_change=log2(T2_result$Fold_change)
library(ggplot2)
cut_off_pvalue = 0.05
cut_off_logFC = log2(1/1.2)
# get meaningful expressions
T2_result$change = ifelse(T2_result$p_list < cut_off_pvalue & abs(T2_result$Fold_change) >= cut_off_logFC,
                        ifelse(T2_result$Fold_change> -cut_off_logFC ,'Up',
                               ifelse(T2_result$Fold_change < cut_off_logFC,'Down','Stable')),
                        'Stable')

# annotate gene labels
T2_result$delabel <- NA
library(stringr)
library(ggrepel)
for (i in 1:nrow(T2_result)){
  if(T2_result$Fold_change[i]!="Stable"){
    T2_result$delabel[i] = 
      str_extract_all(T2_result$SYMBOL[i],"([A-Z0-9]+)")[[1]][1]
  }
} 

# draw volcano graph
p=ggplot(
  data=T2_result, aes(x=Fold_change, y=-log10(p_list),col=change)) +
  geom_point(alpha=0.4, size=2.5) +
  # set color
  scale_color_manual(values=c("blue", "grey","red"))+
  theme_minimal()+
  # add annotation
  geom_label_repel(
    data = subset(T2_result, T2_result$p_list < 0.05 & abs(T2_result$Fold_change) >= cut_off_logFC),
    aes(label = delabel),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  # draw lines
  geom_vline(xintercept=c(-cut_off_logFC, cut_off_logFC),lty=4, col="grey",lwd=0.8) +
  geom_hline(yintercept=-log10(cut_off_pvalue),lty=4, col="grey",lwd=0.8)+
  # draw labels
  labs(x="log2(fold change)",
       y="-log10(p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank())
# save file
ggsave(p,filename = 'result.jpg',dpi = 300, width = 12,height = 9)
# -------------5-----------------------
dif_gene=T2_result$SYMBOL[T2_result$change!="Stable"]
if(!requireNamespace("clusterProfiler",quietly = TRUE))
  install.packages("clusterProfiler")#this might be slow to be loaded
library(clusterProfiler)
if(!requireNamespace("org.Hs.eg.db",quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")#this might be slow to be loaded
library(org.Hs.eg.db)

ego=enrichGO(gene = dif_gene,
             OrgDb='org.Hs.eg.db',
             keyType = "SYMBOL",
             ont = 'BP',
             universe=T2_result$SYMBOL,
             pvalueCutoff = 1,
             pAdjustMethod = "none",
             minGSSize = 10,
             maxGSSize = 20, 
             readable = TRUE)
if(!requireNamespace("GOplot",quietly = TRUE))
  install.packages("GOplot")
library(GOplot)
if(!requireNamespace("enrichplot",quietly = TRUE))
  install.packages("enrichplot")
library(enrichplot)
library(DOSE)
library(ggnewscale)
# draw fancy graphs that I can't depict if I did it again
g1=dotplot(ego,showCategory=30)
ggsave(g1,filename = "dotplot.png",dpi = 600, width = 12,height = 9)
g2=heatplot(ego)# no fold change here or generate randomly
ggsave(g2,filename = "heatplot.png",dpi = 600, width = 12,height = 9)
g3=cnetplot(ego,categorySize="pvalue",circular = TRUE, colorEdge = TRUE)
ggsave(g2,filename = "cneplot.png",dpi = 600, width = 12,height = 9)
g4=emapplot(ego)
# cut redundancy
ego2=simplify(ego,cutoff=0.7,by="p.adjust",select_fun=min)