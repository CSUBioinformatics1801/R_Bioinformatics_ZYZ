# ------------encode utf-8-------------------
#datetime: 20210426
#author: Sean Peldom Zhang
#R.version()=4.0.4
Exp6path="D:/Remain/9.R/EXP/chapter6"
if(getwd()!=Exp6path)setwd(Exp6path)
# ---------------Practice(1)-----------------
# rewrite practice(2) from Exp5.r
dist_t=function(gene_name){
  if(is.character(gene_name)==FALSE){
    print("Gene name must be a character!")
    return(NULL)
  }
  if(!requireNamespace("openxlsx",quietly = TRUE))
    install.packages("openxlsx")
  library(openxlsx)
  if(exists("ADproteomic")==FALSE){
    ADproteomic=read.xlsx("D:/Remain/9.R/EXP/chapter2/data/ADproteomic.xlsx")
  }# to save time from reading
  library(stringr)
  # detect disease
  patstr="(?<=\\.)[a-z]+(?=[0-9]+)"
  disease_names=str_extract(names(ADproteomic),pattern = patstr)
  
  # detect strange names I don't know
  doknstr="[A-Za-z]+(?=\\.)"
  dokn_names=str_extract(names(ADproteomic),pattern = doknstr)
  
  # detect gene name's row location
  judgeACTB=str_extract(ADproteomic$Protein.IDs,"H7BXI1")
  if(length(na.omit(judgeACTB))==0){
    cat("Cannot find",gene_name,"in ADproteomic, please try again!")
    return(NULL)
  }
  cat("Calculating welch 2 t-test between CTL & AD of",gene_name,".\nBe patient, it may takes few centuries....")
  for(ACTBnum in 1:nrow(ADproteomic)){
    if(!is.na(judgeACTB[ACTBnum]))break
  }
  #TODO: need to be optimized if there's meta rows in ADproteomic
  
  # fill data for df_disease
  df_disease=as.data.frame(matrix(nrow=ncol(ADproteomic)-1,ncol=3))
  colnames(df_disease)=c("Expressions","Disease","Don't known")
  ADrownams=names(ADproteomic)
  for(i in 2:length(ADrownams)){
    if(ADproteomic[[ADrownams[i]]][ACTBnum]==0){
      df_disease$Expressions[i-1]=NA
    }else{
      df_disease$Expressions[i-1]=log10(ADproteomic[[ADrownams[i]]][ACTBnum])
    }
    df_disease$Disease[i-1]=disease_names[i]
    df_disease$Dontknown[i-1]=dokn_names[i]
  }
  t.test(df_disease$Expressions[df_disease$Disease=="ctl" & df_disease$Dontknown=="Intensity"],
         df_disease$Expressions[df_disease$Disease=="ad" & df_disease$Dontknown=="Intensity"])
}
dist_t("A7E261")
# ---------------Practice(2)-----------------
# init env
rm(list = ls(all = TRUE))# be careful to run this row
if(!requireNamespace("clusterProfiler",quietly = TRUE))
  install.packages("clusterProfiler")#this might be slow to be loaded
library(clusterProfiler)
if(!requireNamespace("org.Hs.eg.db",quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")#this might be slow to be loaded
library(org.Hs.eg.db)
if(!requireNamespace("openxlsx",quietly = TRUE))
  install.packages("openxlsx")
library(openxlsx)
ADproteomic=read.xlsx("D:/Remain/9.R/EXP/chapter2/data/ADproteomic.xlsx")
if(!requireNamespace("stringr",quietly = TRUE))
  install.packages("stringr")
library(stringr)

# transform ID type
df_up2multi=as.data.frame(matrix(nrow=0,ncol=1))
colnames(df_up2multi)=c("Uniprot")
pattern="(?<=sp\\|)([A-Z0-9]+)(?=\\|)"
for(i in 1:nrow(ADproteomic)){
  df_up2multi=rbind(df_up2multi,data.frame(Uniprot=str_extract_all(ADproteomic$Protein.IDs[i],pattern)[[1]]))
}
Uni_Sym=bitr(df_up2multi$Uniprot,fromType = "UNIPROT",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = FALSE)
Uni_Ent=bitr(df_up2multi$Uniprot,fromType = "UNIPROT",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop = FALSE)
df_up2multi=rbind(Uni_Sym,Uni_Ent$ENTREZID)
# ---------------Practice(3)-----------------
# load result as function from Exp5 HW(1):
get_difprot=function(){
  library(openxlsx)
  library(stringr)
  if(exists("ADproteomic")==FALSE){
    ADproteomic=read.xlsx("D:/Remain/9.R/EXP/chapter2/data/ADproteomic.xlsx")
  }# to save time from reading
  ADrownams=names(ADproteomic)
  # filter data
  ad_subset=subset(ADproteomic, select=na.omit(str_extract(ADrownams,".*LFQ.*ad.*")))
  ctl_subset=subset(ADproteomic, select=na.omit(str_extract(ADrownams,".*LFQ.*ctl.*")))
  judge_3NA=function(x){return (colSums(as.matrix(x) ==0)<3)}# x is a row
  del_judge=apply(ad_subset, 1, judge_3NA)&apply(ctl_subset, 1, judge_3NA)
  ADrownams_del=ADproteomic$Protein.IDs[del_judge]
  ad_rm_3NA=ad_subset[del_judge,]
  ctl_rm_3NA=ctl_subset[del_judge,]
  # t test for each protein
  p_list=c()
  for(i in 1:nrow(ad_rm_3NA)){
    p_list=append(p_list,t.test(ad_rm_3NA[i,],ctl_rm_3NA[i,],var.equal = FALSE)$p.value)
  }
  return(data.frame(ProteinID = ADrownams_del, 'p-value' = p_list))
}
dif_prot=get_difprot()
dif_prot=subset(dif_prot,!is.na(dif_prot$Uniprot)&dif_prot$p.value<0.05)

# extract uniprot p<0.05
pattern="(?<=sp\\|)([A-Z0-9]+)(?=\\|)"
for(i in 1:nrow(dif_prot)){
  dif_prot$Uniprot[i]=str_extract(dif_prot$ProteinID[i],pattern)[[1]]
}
all_prot=as.data.frame(matrix(ncol = 1,nrow = 0))
colnames(all_prot)=c("Uniprot")
for(i in 1:nrow(ADproteomic)){
  all_prot=rbind(all_prot,as.data.frame(str_extract_all(ADproteomic$Protein.IDs[i],pattern)[[1]]))
}
print("Be patient, it may takes few centuries...")
ego=enrichGO(gene = dif_prot$Uniprot,
             OrgDb='org.Hs.eg.db',
             keyType = "UNIPROT",
             ont = 'BP',
             universe=all_prot,
             pvalueCutoff = 0.05,
             pAdjustMethod = "none",
             minGSSize = 10,
             maxGSSize = 500, 
             readable = TRUE)
write.csv(as.data.frame(ego),"G-enrich.csv",row.names =F)


# ---------------Practice(4)-----------------
# load env
library(ggplot2)
library(clusterProfiler)
if(!requireNamespace("GOplot",quietly = TRUE))
  install.packages("GOplot")
library(GOplot)
if(!requireNamespace("enrichplot",quietly = TRUE))
  install.packages("enrichplot")
library(enrichplot)
library(DOSE)
library(ggnewscale)
# draw fancy graphs that I can't depict if I did it again
g1=dotplot(ego)
g2=heatplot(ego)
g3=cnetplot(ego,categorySize="pvalue",circular = TRUE, colorEdge = TRUE)
# cut redundancy
ego2=simplify(ego,cutoff=0.7,by="p.adjust",select_fun=min)
# ego3=data.frame(ego2)
# GO=ego3[1:9,c(1,2,8,6)]
# GO$geneID=str_replace_all(GO$geneID,"/",",")
# names(GO)=c("ID","Term","Genes","adj_pval")
# GO$Category="BP"
# # construct gene matrix, and randomly generate logFC 
# genedata = data.frame(ID=dif_prot$Uniprot,logFC=rnorm(length(dif_prot$Uniprot),mean = 0,sd=2)) 
# circ = circle_dat(GO,genedata)
# ------------encode utf-8-------------------
#datetime: 20210426
#author: Sean Peldom Zhang
#R.version()=4.0.4
Exp6path="D:/Remain/9.R/EXP/chapter6"
if(getwd()!=Exp6path)setwd(Exp6path)
# ---------------Practice(1)-----------------
# rewrite practice(2) from Exp5.r
dist_t=function(gene_name){
  if(is.character(gene_name)==FALSE){
    print("Gene name must be a character!")
    return(NULL)
  }
  if(!requireNamespace("openxlsx",quietly = TRUE))
    install.packages("openxlsx")
  library(openxlsx)
  if(exists("ADproteomic")==FALSE){
    ADproteomic=read.xlsx("D:/Remain/9.R/EXP/chapter2/data/ADproteomic.xlsx")
  }# to save time from reading
  library(stringr)
  # detect disease
  patstr="(?<=\\.)[a-z]+(?=[0-9]+)"
  disease_names=str_extract(names(ADproteomic),pattern = patstr)
  
  # detect strange names I don't know
  doknstr="[A-Za-z]+(?=\\.)"
  dokn_names=str_extract(names(ADproteomic),pattern = doknstr)
  
  # detect gene name's row location
  judgeACTB=str_extract(ADproteomic$Protein.IDs,"H7BXI1")
  if(length(na.omit(judgeACTB))==0){
    cat("Cannot find",gene_name,"in ADproteomic, please try again!")
    return(NULL)
  }
  cat("Calculating welch 2 t-test between CTL & AD of",gene_name,".\nBe patient, it may takes few centuries....")
  for(ACTBnum in 1:nrow(ADproteomic)){
    if(!is.na(judgeACTB[ACTBnum]))break
  }
  #TODO: need to be optimized if there's meta rows in ADproteomic
  
  # fill data for df_disease
  df_disease=as.data.frame(matrix(nrow=ncol(ADproteomic)-1,ncol=3))
  colnames(df_disease)=c("Expressions","Disease","Don't known")
  ADrownams=names(ADproteomic)
  for(i in 2:length(ADrownams)){
    if(ADproteomic[[ADrownams[i]]][ACTBnum]==0){
      df_disease$Expressions[i-1]=NA
    }else{
      df_disease$Expressions[i-1]=log10(ADproteomic[[ADrownams[i]]][ACTBnum])
    }
    df_disease$Disease[i-1]=disease_names[i]
    df_disease$Dontknown[i-1]=dokn_names[i]
  }
  t.test(df_disease$Expressions[df_disease$Disease=="ctl" & df_disease$Dontknown=="Intensity"],
         df_disease$Expressions[df_disease$Disease=="ad" & df_disease$Dontknown=="Intensity"])
}
dist_t("A7E261")
# ---------------Practice(2)-----------------
# init env
rm(list = ls(all = TRUE))# be careful to run this row
if(!requireNamespace("clusterProfiler",quietly = TRUE))
  install.packages("clusterProfiler")#this might be slow to be loaded
library(clusterProfiler)
if(!requireNamespace("org.Hs.eg.db",quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")#this might be slow to be loaded
library(org.Hs.eg.db)
if(!requireNamespace("openxlsx",quietly = TRUE))
  install.packages("openxlsx")
library(openxlsx)
ADproteomic=read.xlsx("D:/Remain/9.R/EXP/chapter2/data/ADproteomic.xlsx")
if(!requireNamespace("stringr",quietly = TRUE))
  install.packages("stringr")
library(stringr)

# transform ID type
df_up2multi=as.data.frame(matrix(nrow=0,ncol=1))
colnames(df_up2multi)=c("Uniprot")
pattern="(?<=sp\\|)([A-Z0-9]+)(?=\\|)"
for(i in 1:nrow(ADproteomic)){
  df_up2multi=rbind(df_up2multi,data.frame(Uniprot=str_extract_all(ADproteomic$Protein.IDs[i],pattern)[[1]]))
}
Uni_Sym=bitr(df_up2multi$Uniprot,fromType = "UNIPROT",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = FALSE)
Uni_Ent=bitr(df_up2multi$Uniprot,fromType = "UNIPROT",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop = FALSE)
df_up2multi=rbind(Uni_Sym,Uni_Ent$ENTREZID)
# ---------------Practice(3)-----------------
# load result as function from Exp5 HW(1):
get_difprot=function(){
  library(openxlsx)
  library(stringr)
  if(exists("ADproteomic")==FALSE){
    ADproteomic=read.xlsx("D:/Remain/9.R/EXP/chapter2/data/ADproteomic.xlsx")
  }# to save time from reading
  ADrownams=names(ADproteomic)
  # filter data
  ad_subset=subset(ADproteomic, select=na.omit(str_extract(ADrownams,".*LFQ.*ad.*")))
  ctl_subset=subset(ADproteomic, select=na.omit(str_extract(ADrownams,".*LFQ.*ctl.*")))
  judge_3NA=function(x){return (colSums(as.matrix(x) ==0)<3)}# x is a row
  del_judge=apply(ad_subset, 1, judge_3NA)&apply(ctl_subset, 1, judge_3NA)
  ADrownams_del=ADproteomic$Protein.IDs[del_judge]
  ad_rm_3NA=ad_subset[del_judge,]
  ctl_rm_3NA=ctl_subset[del_judge,]
  # t test for each protein
  p_list=c()
  for(i in 1:nrow(ad_rm_3NA)){
    p_list=append(p_list,t.test(ad_rm_3NA[i,],ctl_rm_3NA[i,],var.equal = FALSE)$p.value)
  }
  return(data.frame(ProteinID = ADrownams_del, 'p-value' = p_list))
}
dif_prot=get_difprot()
dif_prot=subset(dif_prot,!is.na(dif_prot$Uniprot)&dif_prot$p.value<0.05)

# extract uniprot p<0.05
pattern="(?<=sp\\|)([A-Z0-9]+)(?=\\|)"
for(i in 1:nrow(dif_prot)){
  dif_prot$Uniprot[i]=str_extract(dif_prot$ProteinID[i],pattern)[[1]]
}
all_prot=as.data.frame(matrix(ncol = 1,nrow = 0))
colnames(all_prot)=c("Uniprot")
for(i in 1:nrow(ADproteomic)){
  all_prot=rbind(all_prot,as.data.frame(str_extract_all(ADproteomic$Protein.IDs[i],pattern)[[1]]))
}
print("Be patient, it may takes few centuries...")
ego=enrichGO(gene = dif_prot$Uniprot,
             OrgDb='org.Hs.eg.db',
             keyType = "UNIPROT",
             ont = 'BP',
             universe=all_prot,
             pvalueCutoff = 0.05,
             pAdjustMethod = "none",
             minGSSize = 10,
             maxGSSize = 500, 
             readable = TRUE)
write.csv(as.data.frame(ego),"G-enrich.csv",row.names =F)


# ---------------Practice(4)-----------------
# load env
library(ggplot2)
library(clusterProfiler)
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
# ego3=data.frame(ego2)
# GO=ego3[1:9,c(1,2,8,6)]
# GO$geneID=str_replace_all(GO$geneID,"/",",")
# names(GO)=c("ID","Term","Genes","adj_pval")
# GO$Category="BP"
# # construct gene matrix, and randomly generate logFC 
# genedata = data.frame(ID=dif_prot$Uniprot,logFC=rnorm(length(dif_prot$Uniprot),mean = 0,sd=2)) 
# circ = circle_dat(GO,genedata)
# TODO:bug failed, this is the end 
# # 01. ?????????
# chord &lt;- chord_dat(data=circ, genes = genedata)  # ???????????????????????????????????????
# chord &lt;- chord_dat(data=circ, process = GO$Term) #??????????????????GO term???????????????
# 
# chord &lt;- chord_dat(data=circ, genes=genedata,process = GO$Term) # ????????????
# GOChord (data=chord,
#          title = "GOChord plot", 
#          space = 0.02, # go term???????????????
#          limit = c(3,5), #?????????????????????????????????????????????go term??????,???????????????????????????????????????Go term????????????
#          gene.order ='logFC',gene.space=0.25,gene.size=10,
#          lfc.col = c('firebrick3','white','royalblue3'), #?????????????????????
#          ribbon.col=brewer.pal(length(GO$Term),"set3") ,#GO term??????
#          process.label = 18 # Go terms????????????
# )
# 
# 
# # 02. ?????????
# GOBar(circ,display = 'multiple')
# 
# # 03. ?????????
# # ???????????????,????????????????????????,??????????????????,?????????????????????,?????????????????????:
# GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)
# # ???????????????,??????????????????bg.col?????????TRUE,?????????????????????????????????????????????:
# GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)
# # ????????????????????????????????????????????????reduce_overlap ,???????????????????????????. ???????????????,??????????????????????????????+??????,?????????????????????. ??????,??????????????????????????????,???????????????????????????????????????????????????. ?????????????????????????????????????????????????????????????????????. ?????????????????????GO????????????????????????,????????????????????????:
# # Reduce redundant terms with a gene overlap &gt;= 0.75...
# reduced_circ &lt;- reduce_overlap(circ, overlap = 0.75)
# # ...and plot it
# GOBubble(reduced_circ, labels = 2.8)
# 
# # 04. ??????
# GOCircle(circ)
# 
# # ??????????????????????????????logFC????????????????????????. ???????????????,????????????????????????,????????????????????????. ??????????????????lfc.col????????????. ??????,???????????????,???????????????????????????,????????????????????????z??????????????????. Z?????????????????????????????????????????????. ????????????,???????????????????????????. ???????????????z??????????????????????????????,????????????????????????????????????????????????????????????????????????????????????. ????????????????????????????????????????????????,???????????? GOCirlce.nsub?????????????????????????????????????????????. ??????,?????????????????????????????????. ??????????????????????????????,?????????????????????????????????ID???????????????(???????????????)???
# IDs &lt;- c('GO:0007507', 'GO:0001568', 'GO:0001944', 'GO:0048729', 'GO:0048514', 'GO:0005886', 'GO:0008092', 'GO:0008047')
# GOCircle(circ, nsub = IDs)
# 
# # 05. ??????
# GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))
# 
# # 06. ?????????
# GOCluster(circ, EC$process, clust.by = 'logFC', term.width = 2)
# OCluster(circ, EC$process, clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))
# 
# # 07. ?????????
# l1 &lt;- subset(circ, term == 'heart development', c(genes,logFC))
# l2 &lt;- subset(circ, term == 'plasma membrane', c(genes,logFC))
# l3 &lt;- subset(circ, term == 'tissue morphogenesis', c(genes,logFC))
# GOVenn(l1,l2,l3, label = c('heart development', 'plasma membrane', 'tissue morphogenesis'))
