#######################
## 安装本章需要的R包 ##
#######################
rm(list = ls())
setwd("c:\\workingdirectory")
#例5-0
#安装并升级Bioconductor
source('http://Bioconductor.org/biocLite.R');
biocLite() 
#已安装的R包
pkgs <- rownames(installed.packages())
#本章中需要的R包
packages_1=c("Biobase","CLL","simpleaffy", "affyPLM", "RColorBrewer","affy",
             "gcrma","graph", "GenomicRanges","affycoretools","limma","annotate",
             "hgu95av2.db","GOstats","GeneAnswers","pheatmap","Rgraphviz","GEOquery")
#将要安装的R包
packages_2=setdiff(packages_1,pkgs) 
#安装所需R包
biocLite(packages_2) ;

#############
##  例5-1  ##
#############
#加载所需R包
#CLL包会自动调用affy包，affy包含有后面需要的rma()函数
library (CLL); 
#读入数据（CLL包中附带的示例数据集）
data (CLLbatch) ; 
#调 用RMA算法来对数据预处理（详见5.3.3）
CLLrma <- rma(CLLbatch) ; 
#读取预处理后所有样品的基因（实际上是探针组）表达值
e <- exprs(CLLrma); 
#查 看 部分数据
e[1:5, 1:5]  

#############
##  例5-2  ##
#############
#加载所需R包, 前面已安装
library (CLL); 
#读入数据（CLL包中附带的示例数据集）
data (CLLbatch) ; 
#查看 数 据类型，结果不再显示
data.class (CLLbatch) ; 
#读入所有样品的状态信息
data(disease); 
#查看所有样品的状态信息
disease; 
#查看"AffyBatch"的详细介绍
help (AffyBatch) ; 

#############
##  例5-3  ##
#############
#查看第一张芯片的灰度图像
image(CLLbatch[, 1]);

#############
##  例5-4  ##
#############
#加载所需R包
library(simpleaffy);
library (CLL); 
data (CLLbatch) ; 
#获 取 质 量 分 析 报 告
Data.qc <- qc(CLLbatch);
#图 型 化 显 示 报 告 
plot(Data.qc);  

#############
##  例5-5  ##
#############
#加载所需R包
library (affyPLM);
library (CLL); 
#读入数据（CLL包中附带的示例数据集）
data (CLLbatch);
#对数据集做回归计算，结果是一个"PLMset"类型的对象
Pset<- fitPLM(CLLbatch);
#画第一个芯片数据的原始图
image(CLLbatch[, 1]);
#根据计算结果，画权重图
image(Pset, type = "weights", which = 1, main = "Weights");
#根据计算结果，画残差图
image(Pset, type = "resids", which = 1, main = "Residuals");
#根据计算结果，画残差符号图
image(Pset, type = "sign.resids", which = 1, main = "Residuals.sign");

#############
##  例5-6  ##
#############
#加载所需R包，RColorBrewer 包包含多种预设的颜色集
library(affyPLM);
library(RColorBrewer) ;
library (CLL); 
#读入数据（CLL包中附带的示例数据集）
data (CLLbatch);
#对数据集做回归计算，结果是一个"PLMset"类型的对象
Pset <-fitPLM(CLLbatch);
#载入一组颜色
colors <- brewer.pal(12, "Set3");
#绘 制RLE箱 线 图
Mbox(Pset, ylim = c(-1, 1), col = colors, main = "RLE", las = 3) ;
# 绘 制NUSE箱线 图
boxplot(Pset, ylim = c(0.95, 1.22), col = colors, main = "NUSE", las = 3) ;

#############
##  例5-7  ##
#############
#加载所需R包
library(affy);
library(RColorBrewer);
library (CLL); 
#读入数据（CLL包中附带的示例数据集）
data (CLLbatch);
#获 取 降 解 数 据
data.deg <- AffyRNAdeg(CLLbatch); 
#载入一组颜色
colors <- brewer.pal(12, "Set3");
#绘 制RNA降 解 图
plotAffyRNAdeg(data.deg, col = colors) ; 
#在左上部位加 注 图 注
legend("topleft", rownames(pData(CLLbatch)), col = colors, lwd = 1, inset = 0.05, cex = 0.5);
#从CLL数据集中去除样品CLL1, CLL10和 CLL13
CLLbatch <- CLLbatch[, -match(c("CLL10.CEL", "CLL1.CEL", "CLL13.CEL"), sampleNames(CLLbatch))];

#############
##  例5-8  ##
#############
#加载所需R包
library (CLL); 
library (gcrma); 
library (graph); 
library(affycoretools);
#读入数据（CLL包中附带的示例数据集）
data (CLLbatch);
data (disease);
#使 用gcrma算法来预处理数据
CLLgcrma <- gcrma(CLLbatch) ; 
#提取基因表达矩阵
eset <- exprs(CLLgcrma) ;  
#计算样品两两之间的Pearson相关系数
pearson_cor <- cor(eset) ; 
#得到Pearson距离的下三角矩阵
dist.lower <- as.dist(1 - pearson_cor) ;  
#聚类分析 
hc <- hclust(dist.lower, "ave") ; 
#根据聚类结果画图 
plot(hc) ; 

#PCA分析
samplenames <- sub(pattern="\\.CEL", replacement="",colnames(eset))
groups <- factor(disease[,2])
plotPCA(eset,addtext=samplenames,groups=groups,groupnames=levels(groups))


#############
##  例5-9  ##
#############
#加载所需R包
library(affy);
library (CLL); 
data (CLLbatch) ; 
eset.mas <- expresso(CLLbatch,  bgcorrect.method="mas",  normalize.method="constant",   pmcorrect.method="mas", summary.method="mas");
bgcorrect.methods();
normalize.methods(CLLbatch);
pmcorrect.methods();
express.summary.stat.methods();

#############
## 例5-10  ##
#############
#加载所需R包
library(affy);
library (CLL); 
data(CLLbatch); 
#使用mas方法做背景校正
CLLmas5 <- bg.correct(CLLbatch, method="mas");  
#使用constant方法标准化
data_mas5 <- normalize(CLLmas5, method = "constant");
#查看每个样品的缩放倍数
head(pm(data_mas5)/pm(CLLmas5), 5);
#查看第二个样品的缩放倍数是怎么计算来的
mean(intensity(CLLmas5)[,1])/mean(intensity(CLLmas5)[,2]);

#############
## 例5-11  ##
#############
#加载所需R包
library(affy);
library (gcrma); 
library(affyPLM) ;
library(RColorBrewer);
library (CLL); 
data (CLLbatch) ; 
colors <- brewer.pal(12, "Set3");
#使用MAS5算 法 来 预 处 理 数 据
CLLmas5 <- mas5(CLLbatch) ;  
#使 用rma算 法 来 预 处 理 数 据
CLLrma <- rma(CLLbatch) ; 
##使用gcrma算 法 来 预 处 理 数 据
CLLgcrma <- gcrma(CLLbatch) ;  
#直 方 图
hist(CLLbatch, main = "original", col = colors) ;
legend("topright", rownames(pData(CLLbatch)), col = colors, lwd = 1, inset = 0.05, cex = 0.5, ncol = 3) ;
hist(CLLmas5, main = "MAS 5.0", xlim = c(-150, 2^10), col = colors) ;
hist(CLLrma, main = "RMA", col = colors) ;
hist(CLLgcrma, main = "gcRMA", col = colors) ;
#箱 线 图
boxplot(CLLbatch, col = colors, las = 3, main = "original") ;
boxplot(CLLmas5, col = colors,las = 3, ylim=c(0, 1000), main = "MAS 5.0") ;
boxplot(CLLrma, col = colors, las = 3, main = "RMA") ;
boxplot(CLLgcrma, col = colors, las = 3, main = "gcRMA");


#############
## 例5-12  ##
#############
#加载所需R包
library (gcrma); 
library(RColorBrewer);
library (CLL); 
library (affy); 
data (CLLbatch) ; 
colors <- brewer.pal(12, "Set3");
#使 用gcrma算 法 来 预 处 理 数 据
CLLgcrma <- gcrma(CLLbatch);  
MAplot(CLLbatch[, 1:4], pairs = TRUE, plot.method = "smoothScatter", cex = 0.8, main = "original MA plot");
MAplot(CLLgcrma[, 1:4], pairs = TRUE, plot.method = "smoothScatter", cex = 0.8,main = "gcRMA MA plot");

#############
## 例5-13  ##
#############
#加载所需R包
library(limma) ;
library (gcrma); 
library (CLL); 
data (CLLbatch) ; 
data(disease); 
#从CLL数据集中去除样品CLL1, CLL10和 CLL13
CLLbatch <- CLLbatch[, -match(c("CLL10.CEL", "CLL1.CEL", "CLL13.CEL"), sampleNames(CLLbatch))];
#使 用gcrma算 法 来 预 处 理 数 据
CLLgcrma <- gcrma(CLLbatch);  
# 去除CLLgcrma样品名中的".CEL"
sampleNames(CLLgcrma) <-  gsub(".CEL$", "", sampleNames(CLLgcrma));
# 去除disease中对应样品CLL1, CLL10和 CLL13的记录
disease <- disease[match(sampleNames(CLLgcrma), disease[,"SampleID"]),];
# 构建余下21个样品的基因表达矩阵
eset <- exprs(CLLgcrma) ;
# 提取实验条件信息
disease <- factor(disease[, "Disease"]);
# 构 建 实 验 设 计 矩 阵
design <- model.matrix(~-1+disease);
# 构 建 对比 模 型，比较两个实验条件下表达数据
contrast.matrix <- makeContrasts (contrasts = "diseaseprogres.  - diseasestable", levels = design);
#线性模 型拟 合
fit <- lmFit(eset, design) ;  
# 根据对比 模 型进行差值计算
fit1 <- contrasts.fit(fit, contrast.matrix) ;   
# 贝 叶 斯 检验
fit2 <- eBayes(fit1) ; 
# 生成所有 基因的检验结果报表
dif <- topTable(fit2, coef = "diseaseprogres.  - diseasestable",  n = nrow(fit2), lfc = log2(1.5)) ; 
# 根据P.Value对结果进行筛选，得到全部差异表达基因
dif <- dif[dif[, "P.Value"] < 0.01, ] ; 
# 显 示 结 果 的 前 六 行
head(dif) ;       
design

#############
## 例5-14  ##
#############
# 加载注释工具包
library(annotate) ;
# 获得基因芯片注释包名称
affydb <- annPkgName(CLLbatch@annotation, type = "db");
# 查看基因芯片注释包名称
affydb
## [1] "hgu95av2.db"
# 依据注释包名称hgu95av2.db，加载注释包，必须设定character.only
library(affydb, character.only = TRUE) ;

# 根据每个探针组的ID获取对应的基因Gene Symbol，并作为一个新的列，加到数据框dif最后
dif$symbols <- getSYMBOL(rownames(dif), affydb) ;
# 根据每个探针组的ID获取对应的基因EntrezID，并作为一个新的列，加到数据框dif最后
dif$EntrezID <- getEG(rownames(dif), affydb) ;
# 显 示 结 果 的 前 六 行
head(dif) ; 

#############
## 例5-15  ##
#############
# 加载所需R包
library(GOstats);
# 提取HG_U95Av2芯片中所有探针组对应的EntrezID，注意保证uniq
entrezUniverse <- unique(unlist(mget(rownames(eset), hgu95av2ENTREZID)));
# 提取所有差异表达基因及其对应的EntrezID，注意保证uniq
entrezSelected <- unique(dif[!is.na(dif$EntrezID), "EntrezID"]);
# 设置GO富集分析的所有参数
params <- new("GOHyperGParams", geneIds = entrezSelected, universeGeneIds = entrezUniverse, 
              annotation = affydb, ontology = "BP", pvalueCutoff = 0.001, conditional = FALSE, testDirection = "over");
# 对所有的GOterm根据params参数做超几何检验
hgOver <- hyperGTest(params);
# 生成所有GOterm的检验结果报表
bp <- summary(hgOver) ;
# 同时生成所有GOterm的检验结果文件，每个GOterm都有指向官方网站的链接，可以获得其详细信息
htmlReport (hgOver,  file='ALL_go.html') ;
# 显 示 结 果 的 前 六 行
head (bp) ; 

#############
## 例5-16  ##
#############
# 加载所需R包
library(GeneAnswers) ; 
# 选取dif中的三列信息构成新的矩阵，第一列必须是EntrezID
humanGeneInput <- dif[, c("EntrezID", "logFC", "P.Value")];  
##获 得 humanGeneInput中基因 的 表 达 值
humanExpr <- eset[match(rownames(dif), rownames(eset)), ] ;  
# 前两个数据做列合并，第一列必须是EntrezID
humanExpr <- cbind(humanGeneInput[, "EntrezID"], humanExpr) ; 
# 去 除NA数 据
humanGeneInput <- humanGeneInput[!is.na(humanGeneInput[,  1]), ] ; 
humanExpr <- humanExpr[!is.na(humanExpr[, 1]), ] ;
# KEGG通路的超几何检验
y <- geneAnswersBuilder(humanGeneInput, "org.Hs.eg.db", categoryType = "KEGG", testType = "hyperG", pvalueT = 0.1, geneExpressionProfile = humanExpr, verbose = FALSE) ;
getEnrichmentInfo(y)[1:6,]

#############
## 例5-17  ##
#############
#加载所需R包
library(pheatmap);
# 从基因表达矩阵中，选取差异表达基因对应的数据
selected <- eset[rownames(dif), ] ;
# 将selected矩阵每行的名称由探针组ID转换为对应的基因symbol
rownames(selected) <- dif$symbols;
# 考虑到显示比例，我们只画前20个基因的热图
pheatmap(selected[1:20, ], color = colorRampPalette(c("green", "black", "red"))(100),  fontsize_row = 4, scale = "row", border_color = NA);

#加载所需R包
library(Rgraphviz);
# 显著富集的GO term的DAG关系图，见图5-21
ghandle <- goDag(hgOver) ; 
# 该图巨大，只能取一部分数据构建局部图
subGHandle <- subGraph (snodes=as.character(summary(hgOver)[,1]), graph=ghandle) ;
plot(subGHandle) ;

# 显著富集的KEGG通路的关系图，见图5-22
yy  <-  geneAnswersReadable (y,verbose  = FALSE);
geneAnswersConceptNet (yy, colorValueColumn= "logFC", centroidSize ="pvalue", output = "interactive");
# 显著富集的KEGG通路的热图，见图5-23
yyy  <- geneAnswersSort (yy, sortBy="pvalue");
geneAnswersHeatmap(yyy)
# 输出会话信息
sessionInfo()


#############
## 例5-18  ##
#############
#加载所需R包
library (GEOquery);
library (CLL); 
#设置当前目录
#setwd("c:\\workingdirectory"); 
#读入U133A中22215个共有探针组列表
U133Acols<-read.table("U133Acols");
#得到12个数据集的ID中的数字
numbers=c(15471,5563,3325,9844,5788,6344,10072,13911,1420,8671,5764);
#得到12个数据集在GEO数据库中的ID
GEO_IDs=paste("GSE",numbers,sep = "");
#得到12个数据集下载文件的后缀名
tars=paste(GEO_IDs,"_RAW.tar",sep = "");
#生成一个空变量，用了保存数据标准化后的结果
trainX=c();
#12次循环处理12个数据集
for(i in 1:length(GEO_IDs))
{
  #得到下载数据的目标路径＋名称
  GEO_tar <- paste(GEO_IDs[i],tars[i],sep = "/");
  #下载数据集
  getGEOSuppFiles(GEO=GEO_IDs[i],baseDir = getwd());
  #将当前数据集中的所有样品数据解压到data子目录
  untar(GEO_tar, exdir="data");
  #得到当前数据集中的所有样品对应的数据文件名称
  cels <- list.files("data/", pattern = "[gz]");
  #解压当前数据集中的所有样品对应的数据文件
  sapply(paste("data", cels, sep="/"), gunzip);
  #得到data子目录的全路径
  celpath <- paste(getwd(),"data",sep = "/");
  #转到data子目录，同时保留当前目录到oldWD
  oldWD <- setwd(celpath); 
  #读取当前目录中的所有样品对应的CEL文件
  raw_data <- ReadAffy();   
  #回到工作目录
  setwd(oldWD);
  #删除data子目录中全部文件，否则下次循环会追加写
  unlink("data", recursive=TRUE) ;
  #用RMA算法标准化数据
  rma_data <- rma(raw_data); 
  #得到基因表达矩阵eset
  eset <- exprs(rma_data);
  #基因表达矩阵转置后，选择需要的列
  x<-t(eset)[,as.vector(t(U133Acols))];
  #提取的数据集不断按行向后追加
  trainX=rbind(trainX,x);
}  
# 12个数据集处理完毕，保存到文件trainX中
write.table (trainX, file = "trainX", sep = "\t",row.names = F,col.names = F) ;


#############
## 例5-19  ##
#############
#加载所需R包
library (GEOquery);
#设置当前目录
setwd("c:\\workingdirectory")
#读入U133A中22215个共有探针组列表
U133Acols<-read.table("U133Acols");
#数据集GSE2503的基因表达数据
gds<-getGEO(GEO = "GSE2503", destdir = getwd());
#得到基因表达矩阵eset
eset <- exprs(gds[[1]]);
#基因表达矩阵转置后，选择需要的列
x<-t(eset)[,as.vector(t(U133Acols))]; 
#全部数据转化为以2为底的对数
trainX2<-log2(x) ;
# 结果保存到文件trainX2中
write.table (trainX2, file = "trainX2",sep = "\t",row.names = F,col.names = F) ;

#############
## 例5-20  ##
#############
#加载所需R包
library (Biostrings);
#设置当前目录
setwd("c:\\workingdirectory"); 
#从文件miRNA.tab读入miRNA序列，第一列是序列ID，第二列是序列内容
data1<-read.table("miRNA.tab");
#提取序列内容
seqs=as.character(data1[,2]);
#提取序列ID
names(seqs)=data1[,1];
#用生成RNAStringSet对象，保存为fasta格式文件
Biostrings::writeXStringSet(RNAStringSet(seqs, use.names=TRUE),"miRNA.fa");

#############
## 例5-21  ##
#############
#加载所需R包
library(biomaRt);
library (Biostrings);
#选中"ensembl"数据库
ensembl_mart <- useMart(biomart="ensembl");
#选中"sscrofa"数据集
dataset_pig <-useDataset(dataset="sscrofa_gene_ensembl",mart= ensembl_mart);
#从dataset_pig数据集中根据affy_porcine ID和description信息
idlist <- getBM(attributes=c("affy_porcine","description"), mart=dataset_pig); 
#从dataset_pig数据集中根据affy_porcine ID提取序列
seqs = getSequence(id=idlist["affy_porcine"], type="affy_porcine", seqType="3utr", mart = dataset_pig);
#去除没有序列内容的数据记录
seqs = seqs[!seqs[,1]=="Sequence unavailable",];
#去除没有UTR注释的数据记录
seqs = seqs[!seqs[,1]=="No UTR is annotated for this transcript",];
#提取序列的内容
x=seqs[,1];
#提取序列的ID
names(x)=seqs[,2];
#结果存入文件"UTR3seqs-2.fa"，格式为fasta
writeXStringSet(DNAStringSet(x, use.names=TRUE),"UTR3seqs-2.fa");

#############
## 例5-22  ##
#############
#加载所需R包
library(biomaRt);
library (Biostrings);
pig_affy_IDs<- read.table("pig_affy_IDs");
pig_affy_IDs<- as.character (unlist(pig_affy_IDs));
#列出"sscrofa"数据集的所有特征，才知道包括"ensembl_transcript_id"和"affy_porcine"
id_mapping <- getBM(attributes=c("ensembl_transcript_id","affy_porcine"), filters = "affy_porcine", values = pig_affy_IDs, mart=dataset_pig);
write.table (id_mapping,"emsembl-affy",sep="\t");

#############
## 例5-23  ##
#############
#加载所需R包
library(affycoretools);
library(genefilter);
library(annotate);
library(GOstats);
# 读入所有CEL文件
rawData <- read.affybatch(filenames=list.celfiles());
# 检查是否读入所有文件
sampleNames(rawData);
# 使用RMA算法预处理所有数据
eset <- rma(rawData); 
# 获得基因芯片注释包名称
annoPackage <- paste(annotation(eset), ".db", sep="");
#安装并加载所需对应的基因芯片注释包
source("http://Bioconductor.org/biocLite.R");
biocLite(annoPackage);
library(annoPackage, character.only = TRUE);
# 取出所有探针组的Affymetrix ID
affy_IDs <- featureNames(eset);
#取出数据集的注释信息
an <- annotation(eset);   
# 根据每个探针组的ID获取对应的Gene Symbol
symbols <- as.character(unlist(mget(affy_IDs,get(paste(an, "SYMBOL", sep="")))));
# 所有"NA",转换成""
symbols[is.na(symbols)] <- "";  
# 根据每个探针组的ID获取对应的Gene Name
gene_names <- as.character(unlist(mget(affy_IDs,get(paste(an,"GENENAME", sep="")))));
# 所有"NA",转换成""
gene_names[is.na(gene_names)] <- "";
# 根据每个探针组的ID获取对应的Entrez ID
entrez <- unlist(mget(affy_IDs, get(paste(an, "ENTREZID", sep=""))));
# 根据每个探针组的ID获取对应的UNIGENE ID
unigenes <- as.character(unlist(lapply(mget(affy_IDs, get(paste(an, "UNIGENE", sep=""))),paste, collapse="//")));
# 根据每个探针组的ID获取对应的Refseq  ID
refseqs <- as.character(unlist(lapply(mget(affy_IDs, get(paste(an, "REFSEQ", sep=""))),paste, collapse="//")));
# 数据基因表达矩阵，同时把注释的其它数据库ID，按注释顺序对应到探针组中
out <- data.frame(ProbeID=affy_IDs,Symbol=symbols, Name=gene_names, EntrezGene=entrez, UniGene=unigenes,RefSeq=refseqs, exprs(eset), stringsAsFactors=FALSE);
row.names(out) <- 1:length(affy_IDs);
# 将带注释的数据基因表达矩阵输出
write.table(out, "Expression_table.xls", sep='\t',row.names=F); 
# 六个样品的名称
samples <- c('HEK293_EGFP_monolayer','HEK293_EGFP_sphere','HEK293_CD147_monolayer',
             'HEK293_CD147_sphere','MIAPaCa_2_NC','MIAPaCa_2_A6');
# 五组对比的名称    
compNames = paste(samples[c(1,1,2,3,5)],samples[c(2,3,4,4,6)],sep=' vs '); 
# 设计一个过滤准则，要求每个基因在6个样品中至少表达一次
flt1 <- kOverA(1,6); 
# 五组对比中log2Fold超过1倍的基因算作差异表达基因，输出到五个不同的文件
out <- foldFilt(eset,fold=1,groups=1:6,comps=list(c(1,2),c(1,3),c(2,4),c(3,4),c(5,6)),compNames,text=T,html=F,save=T,filterfun=flt1);


#############
## 例5-24  ##
#############
# GO的3个领域，这里3个领域分开注释
go_domains <- c('BP','CC','MF');
# 提取所有探针组对应的GO信息，没有GO注释的记录为FALSE，
no_go <- sapply(mget(affy_IDs, hgu133plus2GO), function(x) if(length(x) == 1 && is.na(x))TRUE else FALSE);
# 提取所有有GO注释的探针组
affy_IDs <- affy_IDs[!no_go];
# 提取所有探针组对应的"EntrezGene"，并去除重复记录
all_probes <- unique(getEG(affy_IDs, "hgu133plus2"));
# 提取对比1和对比4之间，对比2和对比3之间的差异表达探针组，以及对比5的差异表达探针组
sig_probes.temp <- list(names(which(abs(out[[2]][,1])==1 & abs(out[[2]][,4])==1)),
                        names(which(abs(out[[2]][,2])==1 & abs(out[[2]][,3])==1)),
                        names(which(abs(out[[2]][,5])==1)));
# 定义结果文件的名称  	       
fnames <- c('common probesets for monolayer vs sphere','common probesets for EGFP vs CD147',
            'MIAPaCa_2_NC vs MIAPaCa_2_A6');
# 每次循环处理一次超几何检验，来根据差异表达基因做GO的富集	    
for(i in 1:3){
  sig_probes = unique(getEG(sig_probes.temp[[i]][sig_probes.temp[[i]]%in%affy_IDs], "hgu133plus2"));
  for(j in 1:3){
    # 设定超几何检验参数
    params = new("GOHyperGParams", geneIds=sig_probes, universeGeneIds=all_probes, conditional = TRUE,
                 annotation = annotation(eset), ontology = go_domains[j],pvalueCutoff=.01);
    # 超几何检验
    hypt = hyperGTest(params);
    # 输出html格式的报告文件
    htmlReport(hypt, digits=8, file=paste(fnames[i],'_',go_domains[j],'.html',sep=''));
  }
}

