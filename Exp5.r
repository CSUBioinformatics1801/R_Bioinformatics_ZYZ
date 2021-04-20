# ------------encode utf-8-------------------
#datetime: 20210419
#author: Sean Peldom Zhang
#R.version()=4.0.4
# ---------------Practice(1)-----------------
# set workspace
Exp5path="D:/Remain/9.R/EXP/chapter5"
if(getwd()!=Exp5path)setwd(Exp5path)
# process data
load("diabetes.RData")
Male_GLU=subset(dat,dat$Sex=="Male")$GLU
Female_GLU=subset(dat,dat$Sex=="Female")$GLU
if(var(Male_GLU)!=var(Female_GLU)){
  print("Different variance of GLU between Male and Female!")
  boxplot(Male_GLU, Female_GLU, names=c("Male_GLU","Female_GLU"))
  t.test(Male_GLU, Female_GLU,var.equal = FALSE)
}
# ---------------Practice(2)-----------------

# load data from Exp4.r
rm(list = ls(all = TRUE))
if(!requireNamespace("openxlsx",quietly = TRUE))
  install.packages("openxlsx")
library(openxlsx)
ADproteomic=read.xlsx("D:/Remain/9.R/EXP/chapter2/data/ADproteomic.xlsx")
library(stringr)
# detect disease
patstr="(?<=\\.)[a-z]+(?=[0-9]+)"
disease_names=str_extract(names(ADproteomic),pattern = patstr)

# detect strange names I don't know
doknstr="[A-Za-z]+(?=\\.)"
dokn_names=str_extract(names(ADproteomic),pattern = doknstr)

# detect H7BXI1's row location
judgeACTB=str_extract(ADproteomic$Protein.IDs,"H7BXI1")
for(ACTBnum in 1:nrow(ADproteomic)){
  if(!is.na(judgeACTB[ACTBnum]))break
}

# fill data for df_disease
df_disease=as.data.frame(matrix(nrow=length(ADrownams)-1,ncol=3))
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
# ---------------Practice(3)-----------------
rm(list = ls(all = TRUE))
load("D:/Remain/9.R/EXP/chapter5/diabetes.RData")
aov(GLU~Sex, dat)

# ---------------Homework(1)-----------------
# load data
rm(list = ls(all = TRUE))
if(!requireNamespace("openxlsx",quietly = TRUE))
  install.packages("openxlsx")
library(openxlsx)
library(stringr)
ADproteomic=read.xlsx("D:/Remain/9.R/EXP/chapter2/data/ADproteomic.xlsx")
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
result=data.frame(ProteinID = ADrownams_del, 'p-value' = p_list)
# primary abandoned
# for(i in 1:nrow(ADproteomic)){
#   if(rowSums(as.matrix(ad_subset[i,]) == 0)<3){
#     del_judge[i]=TRUE
#     next
#   }else if(rowSums(as.matrix(ctl_subset[i,]) == 0)<3){
#     del_judge[i]=TRUE
#     next
#   }
# }