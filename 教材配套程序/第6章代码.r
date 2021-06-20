##############
## 准备工作 ##
##############
# 如果是第一次使用bioconductor或者主程序需更新，要先安装核心包
source("http://bioconductor.org/biocLite.R")
biocLite()
# 建立C:\workingdirectory目录


#############
##  例6-1  ##
#############

# 安装并加载所需R包
source('http://Bioconductor.org/biocLite.R');
biocLite("ShortRead") ;
library(ShortRead);
Q=20;
# 计算Q的Sanger分数 (Phred+33)
PhredQuality (as.integer(Q));
# 计算Q的Illumina分数 (Phred+64)
SolexaQuality(as.integer(Q));

#############
##  例6-2  ##
#############
# 以下为fastq文件中一条序列的信息，转存为6-2.fastq至C:\workingdirectory以便于例6-3使用
@HWUSI-EAS100R:123:C0EPYACXX:6:73:941:1973#0/1
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+ HWUSI-EAS100R:123:C0EPYACXX:6:73:941:1973#0/1
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65

#############
##  例6-3  ##
#############
# 加载所需R包，前面已经安装
library(ShortRead);
# 更换工作目录
# 需保证该目录存在，否则会报错
# R中\符号为转义符，比如\n为换行，因此表示路径时为免程序误会，需要使用双斜杠\\或者直接使用反斜杠/
setwd("c:\\workingdirectory");
# 读入FASTQ文件
reads <- readFastq("6-2.fastq");
# 得到质量分数的类型
score_sys = data.class(quality(reads));
# 得到质量分数(字符表示)
qual <- quality(quality(reads));
# 质量分数转为16进制表示	
myqual_16L <- charToRaw(as.character(unlist(qual)));
# 如果是Phred+64分数表示系统
if(score_sys=="SFastqQuality"){
# 显示分数系统类型
cat("The quality score system is Phred+64" ,"\n");	
# 输出原始分数值，16进制转10进制，再减去64
strtoi(myqual_16L, 16L)-64;			
}
# 如果是Phred+33分数表示系统
if(score_sys=="FastqQuality"){
# 显示分数系统类型
cat("The quality score system is Phred+33" ,"\n");	
# 输出原始分数值，16进制转10进制，再减去33
strtoi(myqual_16L, 16L)-33;
}

#############
##  例6-6  ##
#############
# 由于RNA-seq数据量较大，本例建议在Linux服务器（操作系统Fedora）或者4G内存以上的Windows系统运行。
# 安装并加载所需R包
source('http://Bioconductor.org/biocLite.R');
biocLite("ShortRead") ;
biocLite("SRAdb");
biocLite("R.utils");
library(ShortRead);
library(SRAdb);
library(R.utils);
# 下载需要的数据文件
getFASTQfile("SRR921344");
# 解压后，改名为"T2-1.fastq"
gunzip ("SRR921344.fastq.gz", destname="T2-1.fastq");
# 需要分析的数据文件名称
fastqfile="T2-1.fastq";
# 得到质量分析的结果
qa <- qa(dirPath=".", pattern=fastqfile, type="fastq");
# 输出html格式的分析报告
report(qa, dest="qcReport", type="html");

#############
##  例6-7  ##
#############
# 进入linux系统，下载高通量数据SRR987316.sra
# 此处当注意，wget命令的默认下载位置为当前目录，亦可使用-O人为指定位置
wget ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR987/SRR987316/SRR987316.sra
# 将SRA格式的高通量数据转换为FASTQ格式
./fastq-dump SRR987316.sra
# 调用fastqc软件（版本0.10.1）输出质量控制报告
fastqc SRR987316.fastq

#############
##  例6-8  ##
#############
# 安装并加载所需R包
source('http://Bioconductor.org/biocLite.R');
biocLite("Rsamtools");
biocLite("DESeq");
biocLite("edgeR");
library(Rsamtools);
library(DESeq);
library(edgeR);
# 示例数据来自GenomicRanges包中的文件extdata
fls <- list.files(system.file("extdata",package="GenomicRanges"),
    recursive=TRUE, pattern="*bam$", full=TRUE);
bfl <- BamFileList(fls, index=character());
features <- GRanges(
    seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)),
    ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600, 4000, 
        7500, 5000, 5400), width=c(rep(500, 3), 600, 900, 500, 300, 900, 
        300, 500, 500)), "-",
    group_id=c(rep("A", 4), rep("B", 5), rep("C", 2))) ;
olap <- summarizeOverlaps(features, bfl) ;
deseq <- newCountDataSet(assays(olap)$counts, rownames(colData(olap))) ;
edger <- DGEList(assays(olap)$counts, group=rownames(colData(olap))) ;

#############
##  例6-9  ##
#############
setwd("c:/workingdirectory");
# 导入数据，第1列是基因名称用作row.name，最后1列是基因长度，第1行是sample名称
raw.data <- read.table("raw_counts.table",row.names=1);
# 去除第1列长度信息，只保留样本的基因表达值
counts <- raw.data[, 2:dim(raw.data)[2]];
# 取出第1列长度信息
length<-raw.data[, 1];
# 读取lib_size，即每个样品可以对齐到转录组的读段总数
lib_size <- read.table("lib_size.txt");
lib_size <- unlist(lib_size);
# 计算RPKM值
rpkm <- t(t(counts/length)*10^9/lib_size);
# 输出为表的格式
write.table(rpkm,file = "rpkm.table",sep = "\t");
# 计算相关系数
cor_table = cor(rpkm);
# 输出为表的格式
write.table(cor_table,file = "correlations.table",sep = "\t");

##############
##  例6-10  ##
##############
# 安装并加载所需R包
source('http://Bioconductor.org/biocLite.R');
biocLite("DESeq");
library(DESeq);
# 只用100bp的样品进行计算，51bp的样品不选
counts <- cbind(counts[, 1:3],counts[, 7:9]);
# 设置每个样品的实验条件，3个处理，3个对照
conditions=c(rep("T", 3),rep("CK", 3)); 
# 创建CountDataSet对象(DESeq包的核心数据结构)
cds <- newCountDataSet(counts, conditions);
# 估计每个样本的Size Factor，标准化，标准化后数据存于对象cds
cds <- estimateSizeFactors(cds); 
# 显示每个样本的Size Factor，这步可以跳过去       
sizeFactors(cds);                          
# 散度估计
cds <- estimateDispersions(cds, method = "per-condition", sharingMode="maximum"); 
# 显示散度，这步可以跳过去
dispTable(cds);                           
# 检验"T-CK"差异，结果存入对象et
et <- nbinomTest(cds, "T", "CK");
# 输出表达差异分析的结果
write.table(et,file = "T-CK.table",sep = "\t");
# 查看对象cds的内容
cds

##############
##  例6-11  ##
##############
# 指定工作目录，该目录包括所有的数据文件
setwd("C:/workingdirectory");
# 安装并加载所需R包
library(ShortRead);
library(ggplot2);
# 指定转录组数据文件名称
file="Trinity.fasta";
# 读入fasta文件
seqs <- readFasta(file);
# 从网上下载程序源代码文件contigStats.R
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/contigStats.R");
# 统计长度分布等信息
N <- list(seqs=width(seqs)) ;
reflength <- sapply(N, sum) ;
contigStats(N=N, reflength=reflength, style="ggplot2");
stats <- contigStats(N=N, reflength=reflength, style="data");
# 显示统计结果
stats[["Contig_Stats"]];
# 输出统计结果
write.table(t(stats[["Contig_Stats"]]),file="unigenes.dis");
