rm(list = ls())
setwd("c:\\workingdirectory")

#安装本章用到的软件包：         
#source("http://www.bioconductor.org/biocLite.R");
#biocLite("Biostrings");
#biocLite("BSgenome.Hsapiens.UCSC.hg19");                          
#biocLite("hgu133a2probe");
#加载Biostrings包
library(Biostrings);
#加载人类基因组序列数据包
library(BSgenome.Hsapiens.UCSC.hg19);                                         
#加载HG-U133A的探针数据包
library(hgu133a2probe);                                         

#例4.1: 基本操作：互补，反向，反向互补，翻译，转录和逆转录。
#用DNAString生成一个dna对象                                                  
dna<-DNAString("TCTCCCAACCCTTGTACCAGT");
#查看这个对象
dna; 

#将对象dna由DNAString类型转为"RNAString"类型，直接查看内容
Biostrings::dna2rna(dna); 
# 将对象dna中的DNA转录，产生一个"RNAString"类型新对象rna
rna<-transcribe(dna); 
#查看rna内容
rna  

#再转为"DNAString"类型，RNA序列中的Ｕ全部替换为T
rna2dna(rna); 

#对象rna逆转录，得到新对象cD（"DNAString"类型）
cD<-cDNA(rna);

#查看rna的三连密码子
codons(rna);

# rna翻译，产生新对象AA（"AAString"类型） 
AA <-translate(rna);

# 查看AA的内容
AA;  

# dna的互补，又得到一个"DNAString"类型的对象
complement(dna);

# dna的反向互补序列，还是"DNAString"类型的对象
reverseComplement(dna);

# dna的反向序列，还是"DNAString"类型的对象
reverse(dna);

# 例4.2: 统计人类基因组数据中的碱基频率
# 将第22号染色体全序列对有N的地方遮盖，以方便后续步骤时提高工作效率
chr22NoN <-mask (Hsapiens$chr22, "N"); 
# 统计第2号染色体全序列中的所有基础碱基[ATCG]的出现次数
alphabetFrequency(Hsapiens$chr22, baseOnly =TRUE);
#再统计染色体中所有碱基的出现次数
alphabetFrequency(Hsapiens$chr22);
#看看Hsapiens$chr22是否只有基础碱基[ATCG]（字母）
hasOnlyBaseLetters(Hsapiens$chr22);
#显示Hsapiens$chr22中碱基（字母）种类（不含冗余）
uniqueLetters(Hsapiens$chr22);
#计算Hsapiens$chr22中C或G的数量，注意不是CG两连子
GC_content<-letterFrequency(Hsapiens$chr22, letters ="CG");
#查看C或G的数量
GC_content
#计算Hsapiens$chr22中C或G所占的含量（比例）
GC_pencentage<-letterFrequency(Hsapiens$chr22, letters ="CG")/letterFrequency(Hsapiens$chr22, letters ="ACGT");
#查看C或G的含量
GC_pencentage

# 例4.3: 模板匹配，在一组序列中匹配一个模板
#生成连续7个碱基组成的模板
my_pattern = "TATAAAA";
#在chr22NoN中匹配该模板，读者可自己查看结果
mT = matchPattern(my_pattern, chr22NoN);
#计算chr22NoN中匹配该模板的数量
countPattern(my_pattern, chr22NoN);
#在chr22NoN中匹配该模板且允许一个错配
mmT = matchPattern(my_pattern, chr22NoN, max.mismatch =1);
#另一种方法计算匹配的数量，可以看到多匹配了很多
length(mmT);
#观察前5个匹配得到的片段中错配碱基所在的位置
mismatch(my_pattern, mmT[1:5]);
#左侧将要匹配的模板序列
Lpattern <- "CTCCGAG";
#右侧将要匹配的模板序列
Rpattern <- "GTTCACA";
#用左右模板同时匹配Hsapiens$chr22，要求中间的序列长度不能超过500bp
LRsegments<-matchLRPatterns(Lpattern, Rpattern, 500, Hsapiens$chr22);
#查看匹配到的前5条序列
LRsegments[1:5];

# 例4.4: 模板匹配，在一组序列中匹配一组模板（必须长度一样）
#提取所有探针的序列，组成一组模板，存于对象dict
dict<-hgu133a2probe$sequence;
#计算所有探针（序列）的数量
length(dict);
#查看探针的长度nchar(dict)有多少种，只有一种是25
unique(nchar(dict));
#查看dict的前三项内容（探针序列）
dict[1:3];
#用这组探针序列构建DNA词典（模板），允许最大有一个错配
pdict<-PDict(dict, max.mismatch =1);
#用词典匹配Hsapiens$chr22序列
vindex<-matchPDict(pdict, Hsapiens$chr22);
#每个模板（探针序列）对Hsapiens$chr22匹配的个数
count_index<-countIndex(vindex);
#计算所有模板匹配的总数
sum(count_index);
#看看前3个模板的匹配情况， 结果全是0，看来匹配不多
count_index[1:3];
#统计一下匹配数量的分步，可以看到大部分（243903）模板的匹配数都是0
table(count_index);
#从探针序列中提取匹配数最多的模版对应的序列
dict[count_index == max(count_index)];
#用这个序列在Hsapiens$chr22匹配，看看匹配数量，结果和上面统计分布最后一个相同
countPattern("CTGTAATCCCAGCACTTTGGGAGGC", Hsapiens$chr22);

# 例4.5: 搜索回文结构
#计算chr22_pals长度，限定间隔至少40bp
chr22_pals <-findPalindromes(chr22NoN, min.armlength =40, max.looplength =20);
#计算chr22_pals长度                                                                 
nchar(chr22_pals);
#查看前5个找到的回文结构
chr22_pals[1:5];
#查看回文结构序列中的间隔长度
palindromeArmLength(chr22_pals);
#统计回文结构中的所有基础碱基[ATCG]的出现次数
ans<-alphabetFrequency(chr22_pals, baseOnly =TRUE);
#查看基础碱基的频率
ans;

# 例4.6: 序列比对
#用AAString函数生成一个"AAString"对象aa1
aa1 <-AAString("HXBLVYMGCHFDCXVBEHIKQZ");
#用AAString函数生成一个"AAString"对象aa2
aa2 <-AAString("QRNYMYCFQCISGNEYKQN");

#序列全局比对，应用矩阵"BLOSUM62"打分，要求gap open罚分为3
needwunsQS(aa1, aa2, "BLOSUM62", gappen =3);
#用DNAString函数生成一个"DNAString"对象dna1
dna1 <-DNAString("CTCCGAGGGTTTGAATGAT");

#用DNAString函数生成一个"DNAString"对象dna2
dna2 <-DNAString("CTCCGAGTAGCTGGGATTA");

#构建4x4的矩阵DNA打分矩阵，
mat <-matrix(-5L, nrow =4, ncol =4); 
for (i in seq_len(4)) mat[i, i] <-0L;
rownames(mat) <-colnames(mat) <-DNA_ALPHABET[1:4];

#序列全局比对，应用矩阵mat打分，要求gap open罚分为0
needwunsQS(dna1, dna2, mat, gappen =0); 

# 例4.7: 读写序列文件(Fasta和Fastq格式)
# 指定文件的目录(Biostrings安装目录中的extdata子目录)和文件名(someORF.fa)
filepath<-system.file("extdata", "someORF.fa", package ="Biostrings");  
# 显示上面FASTA文件中的数据信息
fasta.info(filepath);
# 读取FASTA文件
x <-readDNAStringSet(filepath); 

# 查看FASTA文件的内容
x;  
# 命名输出文件
out1 <- 'example1.fasta';
# 把序列输出到文件out1，格式还是FASTA
writeXStringSet(x, out1);

# 指定文件的目录和文件名(s_1_sequence.txt)，这个是FASTQ文件
filepath<-system.file("extdata", "s_1_sequence.txt", package ="Biostrings") ; 

# 显示上面FASTQ文件中的数据信息
fastq.geometry(filepath); 
#读取FASTQ文件
x <-readDNAStringSet(filepath, format ="fastq");  
# 查看FASTQ文件，结果不在这里显示
x;  
# 从第1号染色体上按照固定长度（50bp）依次取短序列（read）的起点（向量）
sw_start<-seq.int(1, length(Hsapiens$chr1) -50, by =50);
# 从起点开始取read，每个长度为10bp，注意sw的格式是"XStringViews"
sw<-Views(Hsapiens$chr1, start =sw_start, width =10);
# 变量sw的格式从"XStringViews"转换为"XStringSet"
my_fake_shortreads<-as(sw, "XStringSet");
# 按照"ID"加开头6个数字的格式，得到一组新名称
my_fake_ids<-sprintf("ID%06d", seq_len(length(my_fake_shortreads)));  
# 用新名称替换旧名称
names(my_fake_shortreads) <-my_fake_ids;
# 查看第500000到500005条数据，结果不在这里显示
my_fake_shortreads[500000:500005];
# 命名输出文件
out2 <- 'example2.fastq';
# 把序列输出到文件out2，格式是FASTQ， 但是缺少质量信息
writeXStringSet(my_fake_shortreads, out2, format ="fastq");
# 产生质量信息
my_fake_quals <- rep.int(BStringSet("DCBA@?>=<;"), length(my_fake_shortreads));
#查看my_fake_quals内容，结果不在这里显示
my_fake_quals;
# 命名输出文件
out3 <- 'example3.fastq';
# 把序列输出到文件out3，格式还是FASTQ， 这次含有质量信息
writeXStringSet(my_fake_shortreads, out3, format ="fastq", qualities =my_fake_quals);

# 例4-8
# 安装CLL数据包
#source("http://www.bioconductor.org/biocLite.R");
#biocLite("CLL"); 
# 载入CLL数据包
library(CLL) ; 
# 载入数据（库文件中附带的示例数据）
data(CLLbatch);

# 查看数据内容与结构
phenoData(CLLbatch);


# 例4-9
# 安装biomaRt包
#source("http://www.bioconductor.org/biocLite.R");
#biocLite("biomaRt"); 
# 载入biomaRt包
library(biomaRt) ; 
# 获取当前可用的数据源，一个数据源叫做一个mart
marts <- listMarts(); 
#只查看前几个
head(marts); 
#使用ensembl数据源，如果知道用这个，前面没必要查看所有数据源
ensembl_mart <- useMart(biomart="ensembl"); 
#获取ensembl_mart中可用数据集
datasets <- listDatasets(ensembl_mart); 
#查看前10个
datasets[1:10, ]; 
#使用猪基因组数据集
dataset_pig <- useDataset("sscrofa_gene_ensembl", mart= ensembl_mart);
#获取dataset_pig数据集上可用的筛选器
filters <- listFilters(dataset_pig); 
#只查看前几个，后面没用到任何筛选器
head(filters); 
#获取可选择的属性（列）
attributes <- listAttributes(dataset_pig); 
#只查看前几个
head(attributes); 
#从dataset_pig数据集中提取ensembl_transcript_id和description信息
idlist <- getBM(attributes= c("ensembl_transcript_id", "description"), mart= dataset_pig);
#从dataset_pig数据集中根据ensembl_transcript_id提取序列
seqs = getSequence(id=idlist["ensembl_transcript_id"], type="ensembl_transcript_id", seqType="3utr", mart = dataset_pig);
#去除没有序列内容的数据记录
seqs = seqs[!seqs[, 1]=="Sequence unavailable", ];
#去除没有UTR注释的数据记录
seqs = seqs[!seqs[ ,1]=="No UTR is annotated for this transcript", ];
#提取序列的内容
x=seqs[ ,1];
#提取序列的ID
names(x)=seqs[ ,2];
#结果存入文件"UTR3seqs-1.fa"，格式为fasta
writeXStringSet(DNAStringSet(x, use.names=TRUE),"UTR3seqs-1.fa");
#同时提取对应3'UTR序列的cDNA序列
cDNAseqs = getSequence(id=idlist["ensembl_transcript_id"], type="ensembl_transcript_id", seqType="cdna", mart = dataset_pig);
x=cDNAseqs[ ,1];
names(x)=cDNAseqs[ ,2];
#结果存入文件"UTR3seqs-1.fa"，格式为fasta
writeXStringSet(DNAStringSet(x,  use.names=TRUE), " transcriptom.fasta");


# 例4.10
# 读入解压后的注释文件
probeset <- read.csv("PrimeView.na32.annot.csv", comment.char ="#");
# 查看列名，可以看到每个探针可以对应到多少种公开数据库的ID上
colnames(probeset);
# 只保留probeset两列，探针ID（probeset id）和Entrez数据库的ID(Entrez.Gene)
Id_mapping <- probeset[,c("Probe.Set.ID", "Entrez.Gene")];
# 将Entrez.Gene一列信息中的空白字符（"\\s+"），转为空字符（""）
Id_mapping$Entrez.Gene <- gsub("\\s+","", Id_mapping$Entrez.Gene);
# Entrez.Gene一列信息不为空的数据保留下来
Id_mapping <- Id_mapping[Id_mapping$Entrez.Gene!="---", ];
# 提取所有的Entrez.Gene，
l<-as.character(Id_mapping$Entrez.Gene);
# 通过查看，发现有些探针ID对应多个Entrez.Gene，
head(l[!grepl("^\\d+$", l)]);
# 先将以"///"分割的多个Entrez.Gene分割开
entrez<-strsplit(as.character(Id_mapping$Entrez.Gene),"///");
# entrez中所有对象的名称赋值为探针ID
names(entrez)<-as.character(Id_mapping$Probe.Set.ID);
# 由于无法直接将对象entrez从list格式转换成matrix格式，先将list内的元素转成以
# 探针ID及Entrez.Gene为列名的matrix，然后再合并成一个长表。为了提高效率， 
# 要用到yapply函数并行运算，下面先定义一个yapply函数
# yapply函数的作用是对list操作时可以同时调用其名称以及索引值，
# 该函数表达简洁、功能强大，但语法复杂，初学者可以跳过下面2句
yapply<-function(X,FUN, ...) { 
  index <- seq(length.out=length(X))
  namesX <- names(X) 
  if(is.null(namesX)) namesX <- rep(NA,length(X))
  FUN <- match.fun(FUN) 
  fnames <- names(formals(FUN)) 
  if( ! "INDEX" %in% fnames ){ 
    formals(FUN) <- append( formals(FUN), alist(INDEX=) )   } 
  if( ! "NAMES" %in% fnames ){ 
    formals(FUN) <- append( formals(FUN), alist(NAMES=) )   } 
  mapply(FUN,X,INDEX=index, NAMES=namesX, MoreArgs=list(...)) };
# 调用上面定义好的yapply函数，实现需要的功能
entrez<-yapply(entrez, function(.ele){cbind(rep(NAMES, length(.ele)), gsub(" ","",.ele))});
# 转换成两列的matrix
entrez<-do.call(rbind, entrez);
# 去除重复记录
entrez<-unique(entrez);
# 结果输出，最后得到一个ID映射文件
write.table(entrez, file="primeviewHumanGeneExprs.txt", sep="\t", col.names=F, row.names=F, quote=F);

# 例4.11
# 安装并加载相应的R包
source("http://bioconductor.org/biocLite.R");
biocLite("AnnotationDbi");
library(AnnotationDbi);
# 查看所有的可用模版
available.chipdbschemas();
# 在当前目录下新建叫做primeview的文件夹
dir.create("primeview");
# 根据之前生成的primeviewHumanGeneExprs.txt，生成一个SQLite数据库，注释信息来源自human.db0，模版采用"HUMANCHIP_DB"。这一步耗时相当长
# 这步结束后，可以在primeview的文件夹中看到一个结果文件primeview.sqlite
populateDB("HUMANCHIP_DB", affy=F, prefix="primeview", fileName="primeviewHumanGeneExprs.txt", metaDataSrc=c("DBSCHEMA"="HUMANCHIP_DB", "ORGANISM"="Homo sapiens", "SPECIES"="Human", "MANUFACTURER"="Affymetrix", "CHIPNAME"="PrimeView Human Gene Expression Array", "MANUFACTURERURL"="http://www.affymetrix.com"), baseMapType="eg", outputDir="primeview");
# 生成数据库种子，指明生成模板"HUMANCHIP.DB"，数据包的名称"primeview.db"以及版本号"1.0.0"等相关信息
seed<-new("AnnDbPkgSeed", Package = "primeview.db", Version="1.0.0", PkgTemplate="HUMANCHIP.DB", AnnObjPrefix="primeview");
# 在primeview的文件夹， 生成最后的primeview.db文件
makeAnnDbPkg(seed, file.path("primeview", "primeview.sqlite"), dest_dir="primeview");



# 例4.12
# 安装并加载相应的R包
#source("http://bioconductor.org/biocLite.R");
#biocLite("AnnotationDbi");
#biocLite("AnnotationForge");
library(AnnotationDbi);
library(AnnotationForge);
# 查看所有的可用模板
available.chipdbschemas();
# 在当前目录下新建叫做primeview的文件夹
dir.create("primeviewdb");
# 直接使用affy的注释文件，通过自动处理，在primeviewdb文件夹中生成sqlite文件
# 在primeviewdb的文件夹，生成最后的primeview.db文件
makeDBPackage("HUMANCHIP_DB", affy=TRUE, prefix="primeviewdb", fileName="PrimeView.na32.annot.csv", baseMapType="eg", outputDir="primeviewdb", version="1.0.0", manufacturer="affymetrix", chipName="PrimeView Human Gene Expression Array", manufacturerUrl="http://www.affymetrix.com");
# 输出本章的会话信息
sessionInfo();
