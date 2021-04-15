# ------------encode utf-8-------------------
#datetime: 20210330
#author: Sean Peldom Zhang
#R.version()=4.0.3
# ---------------Experiments-----------------
# 1.install packages
if(!requireNamespace("BiocManager",quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version="3.12")

if(!requireNamespace("Biostrings",quietly = TRUE))
  BiocManager::install("Biostrings")

if(!requireNamespace("openxlsx",quietly = TRUE))
  install.packages("openxlsx")

# 2. path 
getwd()
setwd("D:/")
dir.exists("Remain")
dir.create("D:/???????????????")
setwd("D:/???????????????")

# 2.2 read files
c2path="D:/Remain/9.R/EXP/chapter2/data"
if(getwd()!=c2path)setwd(c2path)
x1=read.delim("ADdata.txt")
x2=read.csv("ADdata.csv")
x4=read.table("ADdata2.txt",header = TRUE)
load("ADinf.RData")
library(openxlsx)
x3=read.xlsx("ADdata.xlsx")

# 2.3 save files
save.image(file = "class2_1.RData")
save(x1,file =  "class2_1.RData")
write.xlsx(x1,file = "x1.xlsx")
write.table(x2,file = "x2.txt",row.names = FALSE)
write.table(x3,file = "x3.txt",row.names = FALSE,seq="\t",quote = FALSE)
write.table(x4,file = "x4.csv",row.names = FALSE)

# 2.4 browse files
list.files()

# 2.5 read fasta
library(Biostrings)
protein=readAAStringSet("Protein.fasta")
print(protein)

# 3.1 regulation expression
ID=ADdatainf$Majority.protein.IDs
seqID=strsplit(ID[1:5],";")
seqID=unlist(seqID)
ID1=sapply(ID, "[",1)

# 3.2 gsub()/sub()
IDstr<-strsplit(ID,";")
t1=gsub("sp\\|(.*)\\|NUD11.*","\\1",ID[1])
t2=sub("sp\\|(.*)\\|NUD11.*","\\1",ID[1])
t3=gsub("(sp|tr)\\|([A-Z0-9]+).*","\\1",ID[1:5])
t4=sub("(sp|tr)\\|([A-Z0-9]+).*","\\2",ID[1:5])
x_nosort=merge(x1,x2,by="Protein.IDs")
x_sort=merge(x1,x2,by="Protein.IDs",sort = FALSE)

# 4 data extract and consolidate
library(openxlsx)
sheets = list("sheet1" = x1,"sheet2" = x2)
write.xlsx(sheets,paste(c(getwd(),"/output.xlsx"),collapse =""))
pos=match(x3$Protein.IDs,x4$Protein.IDs)
x=cbind(x3,x4[pos,-1])

# ---------------HOMEWORK-----------------
# install stringr
if(!requireNamespace("stringr",quietly = TRUE))
  install.packages("stringr")
library(stringr) 

# install purrr
if(!requireNamespace("purrr",quietly = TRUE))
  install.packages("purrr")
library(purrr) 

# homework (1)
replaceNA=function(IDip,i){
  if(is.na(str_extract_all(IDip[i],"(?<=\\|)([A-Z0-9]+)(?=\\|)")[[1]][2])){
    return("na")
  }else{
    return(IDip[i])
  }
}
ID2=ID
for (i in 1:length(nchar(ID))) {
  ID2[i] <- replaceNA(ID,i)
}

# homework (2)
extract_u3=function(IDip,i){
  uniprot3=str_extract_all(IDip[i],"(?<=\\|)([A-Z0-9]+)(?=\\|)")[[1]][3]
  if(is.na(uniprot3)){
    return("na")
  }else{
    return(uniprot3)
  }
}
ID3=ID
for (i in 1:length(nchar(ID))) {
  ID3[i] <- extract_u3(ID,i)
}

# homework (3)
bind_all=function(x1,x2){
  pos=match(x1$Protein.IDs,x2$Protein.IDs)
  return(cbind(x1,x2[pos,-1]))
}
sum_all_x=bind_all(bind_all(x1,x2),bind_all(x3,x4))
library(openxlsx)
write.xlsx(sum_all_x,file = "ADproteomic.xlsx")

# homework (4)
setwd("D:/Remain/9.R/EXP/chapter2/data/GSE67835")
file_list=list.files()
file_list=lapply(file_list, function(x) x[!grepl("[^csv]$", x)])
if(length(file_list)==0){
  stop("Empty file!")
} else if (sum(lengths(file_list))==0){
  stop("No csv files!")
}
# initiate sum_csv
i=1
while(lengths(file_list)[i]==0){
  # find the first csv index
  i=i+1;
}
sum_data=matrix(unlist(read.csv(file_list[[i]], sep="\t",header = FALSE)[2]))
column_name=matrix(unlist(read.csv(file_list[[i]], sep="\t",header = FALSE)[1]))
# keep files that end with "csv", starts from the second csv
initial_i=i;
for (t in 1:length(file_list)){
  if(length(file_list[[t]])!=0 && t>initial_i){
    # which means a name of a csv
    sum_data=cbind(sum_data,matrix(unlist(read.csv(file_list[[t]],sep="\t",header = FALSE)[2])))
  }
}
raw_name=t(matrix(file_list))
raw_name[1]="gene_name"
# if(length(raw_name)!=139){
#   raw_name[length(raw_name)]=NULL
# }
sum_csv=rbind(raw_name,cbind(column_name,sum_data))
# write.csv(sum_csv,file = "sum_csv.csv",row.names = F,col.names = NA)
save(sum_csv,file = "D:/Remain/9.R/EXP/chapter3/GSE67835.RData")
