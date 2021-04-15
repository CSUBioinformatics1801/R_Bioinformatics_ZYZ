# ------------encode utf-8-------------------
#datetime: 20210330
#author: Sean Peldom Zhang
#R.version()=4.0.4 
# ---------------HOMEWORK(1)-----------------
# initiate .dat from Exp1.r
b <- factor(c(1, 1, 2, 2, 3, 3, 3, 1, 2));
levels(b) <- c("Male", "Female", "unknown");
dat <- data.frame(Sex = b,
                  GLU = c(8.8, 6.5, 5.4, 5.6, 6.7, 7.8, 4.5, 9.7, 5.0),
                  GHb = c(5.5, 3.4, 5.5, 4.3, 3.5, 3.4, 5.5, 7.0, 4.3)
                  );
# add random patients
for (i in 1:20){
  random_patients=data.frame(
    Sex=c(sample(levels(b),1)),
    GLU=c(round(runif(1,3.89,6.1),1)),
    GHb=c(round(runif(1,4,6),1))
  );
  dat=rbind(dat,random_patients)
}
Exp3path="D:/Remain/9.R/EXP/chapter3"
if(getwd()!=Exp3path)setwd(Exp3path)
save(dat,file = paste(c(getwd(),"/output.RData"),collapse =""))

# ---------------HOMEWORK(2)-----------------
library(openxlsx)
ADproteomic=read.xlsx("D:/Remain/9.R/EXP/chapter2/data/ADproteomic.xlsx")
AD_rnames=names(ADproteomic)

# define empty data frames for storing
Intensity_rnames=as.data.frame(matrix(nrow=nrow(ADproteomic),ncol=1))
iBAQ_rnames=as.data.frame(matrix(nrow=nrow(ADproteomic),ncol=1))
LFQ_rnames=as.data.frame(matrix(nrow=nrow(ADproteomic),ncol=1))
Intensity_rnames$V1=NULL
iBAQ_rnames$V1=NULL
LFQ_rnames$V1=NULL

# detect columns in ADproteomic
library(stringr)
Intensity_judge=str_detect(AD_rnames, 'Intensity')
iBAQ_judge=str_detect(AD_rnames, "iBAQ")
LFQ_judge=str_detect(AD_rnames, "LFQ")
for (i in 1:length(AD_rnames)){
  if (Intensity_judge[i]==TRUE){
    Intensity_rnames=cbind(Intensity_rnames,ADproteomic[AD_rnames[i]])
  }
  else if(iBAQ_judge[i]==TRUE){
    iBAQ_rnames=cbind(iBAQ_rnames,ADproteomic[AD_rnames[i]])
  }
  else if(LFQ_judge[i]==TRUE){
    LFQ_rnames=cbind(LFQ_rnames,ADproteomic[AD_rnames[i]])
  }
}
save(Intensity_rnames,
     iBAQ_rnames,
     LFQ_rnames,
     file = paste(c(getwd(),"/ADproteomic.RData"),collapse ="")
     )

# ---------------HOMEWORK(3)-----------------
