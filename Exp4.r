# ------------encode utf-8-------------------
#datetime: 20210412
#author: Sean Peldom Zhang
#R.version()=4.0.4
# ---------------Practice(1)-----------------
windows()
head(cars)
plot(cars[,1],
     cars[,2],
     pch=16,
     col="red",
     xlab = "speed",
     ylab = "dist",
     type = "p",
     font.axis=2,
     font.lab=2,
     cex.lab=1.5)
x=seq(5,25,by=0.1)
y=10*sin(x)+50
lines(x,y,col="green",lwd=2)
a=5:25
b=2*a+10*cos(a)
points(a,b,pch=15,col="blue")
legend("topleft",
       pch=c(16,-1,15),
       lty=c(-1,1,-1),
       col=c("red","green","blue"),
       legend = c("plot with point","lines","points")
       )
#
# ---------------Practice(2)-----------------
ggplot(data=mpg,mapping=aes(x=unlist(mpg[,"cty"]),y=unlist(mpg[,"hwy"])))+
  geom_point(alpha=0.4, size=2.5)+
  geom_line(color="blue",size=1)+
  labs(x="cty",y="hwy")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
          legend.position="right",
          legend.title = element_blank())
# ---------------HOMEWORK(1)-----------------
Exp4path="D:/Remain/9.R/EXP/chapter4"
if(getwd()!=Exp4path)setwd(Exp4path)




# ---------------HOMEWORK(2)-----------------
datapath="D:/Remain/9.R/EXP/chapter4/D_volcano.RData"
load(datapath)
if(!requireNamespace("ggplot2",quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)
cut_off_pvalue = 0.05
cut_off_logFC = log2(1/1.2)
# get meaningful expressions
prostat$change = ifelse(prostat$P < cut_off_pvalue & abs(prostat$FC) >= cut_off_logFC,
                        ifelse(prostat$FC> -cut_off_logFC ,'Up',
                               ifelse(prostat$FC < cut_off_logFC,'Down','Stable')),
                        'Stable')

# annotate gene labels
prostat$delabel <- NA
if(!requireNamespace("stringr",quietly = TRUE))
  install.packages("stringr")
library(stringr)

if(!requireNamespace("ggrepel",quietly = TRUE))
  install.packages("ggrepel")
library(ggrepel)

for (i in 1:nrow(prostat)){
  if(prostat$change[i]!="Stable"){
    prostat$delabel[i] = 
      str_extract_all(prostat$ID[i],"(?<=\\|)([A-Z0-9]+)(?=\\|)")[[1]][1]
  }
} 

# draw volcano graph
p=ggplot(
  data=prostat, aes(x=FC, y=-log10(P),col=change)) +
  geom_point(alpha=0.4, size=2.5) +
# set color
  scale_color_manual(values=c("blue", "grey","red"))+
  theme_minimal()+
# add annotation
  geom_label_repel(
    data = subset(prostat, prostat$P < 0.05 & abs(prostat$FC) >= cut_off_logFC),
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
saveloc=paste(c(substr(datapath,0,str_locate_all(datapath,'\\.')[[1]][2]-1),'.png'),collapse ="")
ggsave(p,filename = saveloc,dpi = 600, width = 12,height = 9)

# ---------------HOMEWORK(3)-----------------
library(openxlsx)
library(stringr)
datapath="D:/Remain/9.R/EXP/chapter2/data/ADproteomic.xlsx"
ADproteomic=read.xlsx(datapath)
# detect disease
patstr="(?<=\\.)[a-z]+(?=[0-9]+)"
disease_names=str_extract(names(ADproteomic),pattern = patstr)
# disease_names=na.omit(unique(str_extract(names(ADproteomic),pattern = patstr)))

# detect strange names I don't know
doknstr="[A-Za-z]+(?=\\.)"
dokn_names=str_extract(names(ADproteomic),pattern = doknstr)
# dokn_names=na.omit(unique(str_extract(names(ADproteomic),pattern = doknstr)))
# dokn_names=dokn_names[dokn_names!="Protein"]

# detect ACTB's row location
judgeACTB=str_extract(ADproteomic$Protein.IDs,"ACTB")
for(ACTBnum in 1:nrow(ADproteomic)){
  if(!is.na(judgeACTB[ACTBnum]))break
}

# fill data for df_disease
df_disease=as.data.frame(matrix(nrow=length(ADrownams)-1,ncol=3))
colnames(df_disease)=c("Expressions","Disease","Don't known")
for(i in 2:length(ADrownams)){
  if(ADproteomic[[ADrownams[i]]][ACTBnum]==0){
    df_disease$Expressions[i-1]=NA
  }else{
    df_disease$Expressions[i-1]=log10(ADproteomic[[ADrownams[i]]][ACTBnum])
  }
  # =log10(ADproteomic[[ADrownams[i]]][ACTBnum])
  df_disease$Disease[i-1]=disease_names[i]
  df_disease$`Don't known`[i-1]=dokn_names[i]
}

# some abandoned codes in list:
# disease_x_ACTB=list()
# for(i in 1:length(disease_names)+1){
#   disease_x_ACTB[i]=c()
# }
# ADrownams=names(ADproteomic)
# for(i in 1:length(ADrownams)){
#   for(t in 1:length(disease_names)){
#     if(!is.na(str_extract(ADrownams[i],disease_names[t]))){
#       disease_x_ACTB[[t]]=append(disease_x_ACTB[[t]],ADproteomic[[ADrownams[i]]][ACTBnum])
#     }
#   }
# }
# # change list to df
# for(i in 1:length(disease_x_ACTB)){
#   df_disease[[1]][i]=disease_names[i]
#   df_disease[[2]][i]=mean(disease_x_ACTB[[i]])
# }

# draw hist graph
library(ggplot2)
p=ggplot(data = subset(df_disease,!is.na(Expressions)),
         mapping = aes(x=`Don't known`,y=Expressions,fill=factor(Disease)))+
  geom_boxplot() +
  xlab("Labels that I don't know")+
  ylab("log10(expressions)")+
  geom_point(position = position_jitterdodge())+
  labs(title =  bquote('ACTB'))+
  theme(plot.title = element_text(hjust = 0.5))
p

# save graph
saveloc=paste(c(substr(datapath,0,str_locate_all(datapath,'\\.')[[1]][2]-1),'ACTB.png'),collapse ="")
ggsave(p,filename = saveloc,dpi = 600, width = 12,height = 9)