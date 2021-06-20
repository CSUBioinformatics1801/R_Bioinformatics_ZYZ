#####################################
#3.2 用R包（非Bioconductor）实现课题#
#####################################

#——————————————————————————————————————
##################################
#A. 定义序列导入函数'seq_import'#
##################################
seq_import<- function(input_file) {
  # 逐行读取数据，并存入向量my_fasta，向量每个元素对应文件input_file中的一行，
  #这样以后可以通过操作向量my_fasta，来操作对应文件的行。
  my_fasta<- readLines(input_file);
  
 # 判断my_fasta中每个元素第一个字母是否是“>”（表示一个fasta记录的注释
#行），判断结果用1和-1表示，并存入向量y。
  y <- regexpr("^>", my_fasta, perl = T);

  # 向量y中为1的元素替换为0，即序列行对应-1，注释行对应0。
  # 这行语句只是一个习惯问题，不是必须的。
  y[y == 1] <- 0;

  # 用index记录下y中全部0的在向量中的位置，对应注释行的行号。
  index<- which(y == 0) ;

  # 生成数据框distance，包括第1列start（除最后一个fasta记录外的所有注释 
#行的位置）和第2列end（除第一个fasta记录外的所有注释行的位置）。
  distance <- data.frame(start = index[1:(length(index) - 1)], end = index[2:length(index)]);

  # 在数据框distance最后增加一行（两个元素），第1个是最后一个fasta记录的
  #注释行位置，第2个是为所有行的行数+1）。
  distance<- rbind(distance, c(distance[length(distance[, 1]), 2], length(y) + 1));

  # 在数据框distance后面加1列，其值是第2列和第1列之差，注释行之间的距离，
#实际上就是每条序列记录对应的行数。
  distance <- data.frame(distance, dist = distance[, 2] - distance[, 1]);

  # 建立从1开始的连续正整数向量，长度等于注释行的数量。
  seq_no<- 1:length(y[y == 0]);

  # 重复正整数向量seq_no中的每一个元素，重复次数为两个临近注释行之间的距离 
#（即distance[, 3]）。
  index<- rep(seq_no, as.vector(distance[, 3]));

  # 建立一个新的数据框变量，名称还是my_fasta，包括3列内容，第1列是index，
#第2列是y，第3列是旧的my_fasta。
  my_fasta<- data.frame(index, y, my_fasta);

  # 数据框my_fasta中，第2列为0的元素，对应的第1列赋值为0。
  my_fasta[my_fasta[, 2] == 0, 1] <- 0;

  # tapply函数调用paste函数的字符串连接功能，把my_fasta[, 3]中的同一类
#元素合并，my_fasta[, 3]的类别由对应my_fasta[, 1]的数据来决定，如“0”表示
#序列所有的注释行，“1”表示第一条记录的序列内容，以此类推。
  seqs <- tapply(as.vector(my_fasta[, 3]), factor(my_fasta[, 1]), paste, collapse ="", simplify = F);

  # 将变量seq由数组类型转化为字符串向量，不包括第1个元素（所有注释行），剩下
#的内容为所有记录的序列。
  seqs <- as.character(seqs [2:length(seqs)]);

  # 从my_fasta[, 3]中提取所有的注释行，存入向量Desc。
  Desc<- as.vector(my_fasta[c(grep("^>", as.character(my_fasta[, 3]), perl =TRUE)), 3]);

  # 建立一个新的数据框变量，名称还是my_fasta，每行对应一个序列记录，包括3列信息（序列的注释，长度和序列内容）。
  my_fasta<- data.frame(Desc, Length =nchar(seqs), seqs);

  # 从my_fasta第一列的注释行中提取序列的ID(Accession Number)。
  Acc<- gsub(".*gb\\|(.*)\\|.*", "\\1", as.character(my_fasta[, 1]), perl = T);

  # 将字符串向量Acc添加到数据框左边，成为一列。
  my_fasta<- data.frame(Acc, my_fasta);

  # 将my_fasta返回，这是习惯性的，R把最后出现的数据作为返回值。
  my_fasta;
}
#—————————————————————————————————————






#_____________________________________________________________________________
#######################################
#B. 定义模式匹配函数'pattern_match'#
#######################################
pattern_match<- function(pattern, sequences, hit_num) {

  # 获取正则表达式pattern表示的模序在所有序列中出现的位置（未找到匹配将返回
#-1），所有位置存入一个列表对象pos，perl=T表示兼容perl的正则表达式格式。
  pos<- gregexpr(pattern, as.character(sequences[, 4]), perl= T);

  # lapply函数调用paste函数的字符串连接功能，对pos中的每个成员的第一个元素操作，即用
#逗号连接成一个字符串，再用unlist将所得的列表转换为向量posv。
  posv<- unlist(lapply(pos, paste, collapse =", "));

  # 将向量posv中值为-1的项赋值为0，即表示该序列中未找到模序pattern。
  posv[posv == -1] <- 0;

  # lapply函数调用自定义函数function，根据pos中的每一个元素，计算
#pattern在每条序列中匹配的个数，再由unlist函数将结果转变为向量。
  hitsv<- unlist(lapply(pos, function(x) if (x[1] == -1) {0} else {length(x)})); 
 
  # 产生一个数据框类型的结果sequences，保留了原来sequences数据的第1、2、
#3、4列，又插入了2列，即匹配位点（Position）和匹配次数（Hits）。
  sequences <- data.frame(sequences[, 1:3], Position = as.vector (posv), Hits =hitsv, sequences[, 4]);

  # 找出匹配次数大于hit_num的序列，并将大写形式替换为小写，gsub中第一个参数
#[A-Z]匹配任意大写字母，“\\L\\1”表示将前面小括号中匹配的任意字母替换为其小写形式。
  tag <- gsub("([A-Z])", "\\L\\1", as.character(sequences[sequences[, 5]> hit_num, 6]), perl = T, ignore.case = T);

  # 为模序pattern加上小括号，以适合perl正则表达式格式，方便下面使用。
  pattern2 = paste("(", pattern, ")", sep ="");

  # 将tag序列中，和模序pattern匹配的部分替换为大写，原理同上，“\\U\\1”表示
#替换为大写。
  tag<- gsub(pattern2, "\\U\\1", tag, perl = T, ignore.case = T);

  # 找出匹配次数大于hit_num的序列，并将序列内容替换为tag中的序列内容，存于
#数据框export。
  export<- data.frame(sequences[sequences[, 5] > hit_num,-6],tag);

  # Acc号前添加fasta格式标识“>”，得到数据框export，第1列是Acc，第2列
#是小写字母表示的蛋白质序列（模式用大写表示）
  export<- data.frame(Acc =paste(">", export[, 1], sep =""), seq = export[,6]);   
                                                               
  # 数据框export矩阵转置输出，到文件Hit_sequences.fasta（fasta文件格式）
  write.table(as.vector(as.character(t(export))), file = "Hit_sequences.fasta", quote = F, row.names = F, col.names = F);

  # 输出提示信息
  cat("含有模序\"", pattern, "\"超过", hit_num, "个的所有蛋白质序列已写入当前工作目录下文件'Hit_sequences.fasta'", "\n", sep ="");

  # 选中匹配次数（sequences的第5列）大于hit_num的序列
  selected<- sequences[sequences[, 5] >hit_num, ];

  # 输出提示信息
  cat("极端嗜盐古菌蛋白组中以下序列含有模序\"", pattern, "\"的数量超过2个：", "\n", sep ="");

  # 输出选中序列的第1到5列到终端，第6列是序列内容太长，不显示
  print(selected[, 1:5]);

  # 返回选中序列
  selected;
}
#_____________________________________________________________________________








#_____________________________________________________________________________
############################################
#C. 定义氨基酸含量统计函数'getAApercentage'#
############################################
getAApercentage<- function(sequences) {
  # 生成一个包含20种标准氨基酸单字母简写的数据框AA。
  AA <- data.frame(AA =c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"));

  # strsplit函数将序列内容sequences[, 6])转换成字符数组，lapply函数调用 
#table函数统计每条序列中各字符（氨基酸）出现的次数。
  AAstat<- lapply(strsplit(as.character(sequences[, 6]), ""), table);

  # 下面循环每次处理一条序列，全部序列共length(AAstat)条。
  for (i in 1:length(AAstat)) { 
    # 计算每条序列中20种氨基酸出现的百分比。
AAperc<- AAstat[[i]]/sequences[, 3][i] * 100;

    # 转换为数据框，第1列是氨基酸种类，第2列是百分比含量。
AAperc<- as.data.frame(AAperc);

    # 将数据框第2列名称改为第i条序列的Acc。
names(AAperc)[2] <- as.vector(sequences[i, 1]);

    # 通过AA中的列名为“AA”的列和AAperc中列名为“Var1”的列之间的元素同名映射合并#AA和AAperc，并产生一个新的对象AA保存结果，这样做实质就是按照“AA”的20种氨基       
    #酸的顺序不断添加在每个序列中的分布数据，每个循环至少一列。
    AA <- merge(AA, AAperc, by.x ="AA", by.y ="Var1", all = T);
  }#循环结束。

  # 将AA中氨基酸种类或百分比为“NA”的项赋值为0。
  for (i in 1:length(AA[[1]])) { #外循环总次数是20（种氨基酸）。
    for (j in 1:length(AA)) {    #内循环总次数是序列总数+1
      if (is.na(AA[i, j])) {     #如果发现“NA” 。
        AA[i, j] <- 0;           #替换为0。
      }
    }
  }#循环结束。

  # 统计所有序列中每种氨基酸出现的平均百分比，放入AA最后一列。
  AApercentage <- data.frame(AA, Mean =apply(AA[, 2:length(AA)], 1, mean, na.rm = T));

  # 将对象AApercentage输出到同名的csv文件。
  write.csv(AApercentage, file ="AApercentage.csv", row.names = F, quote = F) ;

  # 提示计算完成。
  cat("氨基酸分布数据已经写入当前工作目录下的文件'AApercentage.csv'", "\n");

  # 返回AApercentage。
  AApercentage;
}
#_____________________________________________________________________________







#_____________________________________________________________________________
####################################
#D. 定义两两比对函数'seq_alignment'#
####################################
seq_alignment<- function(sequences) {
      # shell可调用操作系统命令，命令以字符串形式给出，del为windows系统上的删除命令，/f选项表示强制删除只读文件，my_needle_file为所要删除的文件名，这样做的目的是删除上次程序运行的结果文件，否则本次运行结果会追加写入上次的结果文件。
	shell("del /f my_needle_file");

      # 下面循环每次写一条序列存入file1，另一条存入file2，然后调用needle程序做比对，这样每次都是对比两条序列，结果追加写入结果文件。
	for (i in 1:length(sequences[, 1])) {

      # 第1条序列写入file1（fasta格式）。
		cat(as.character(paste(">", as.vector(sequences[i, 1]), sep ="")), as.character(as.vector(sequences[i, 6])), file ="file1", sep ="\n");

		for (j in 1:length(sequences[, 1])) {
                    # 第2条序列写入file2（fasta格式）。
			cat(as.character(paste(">", as.vector (sequences[j, 1]), sep ="")), as.character(as.vector(sequences[j, 6])), file ="file2", sep ="\n");

                    # 调用needle程序对比file1和file2中的序列，结果追加写入文件“my_needle_file”。
			shell("needle file1 file2 stdout -gapopen 10.0 -gapextend 0.5 >> my_needle_file");
		}
	}
      # 提示结果
	cat("Needle程序完成所有序列的两两比对，结果存入文件\"my_needle_file\"\n");
}   
#_____________________________________________________________________________




#____________________________________________________________________________
#########################################
#E1. 定义函数'getScoreMatrix'求得分矩阵#
########################################
getScoreMatrix<- function(sequences) {
    # 读取my_needel_file中的所有行，存入向量score。
score <- readLines("my_needle_file");

    # 查找以“# Score”开头的行（如# Score: 290.5），存入向量score。
score <- score[grep("^# Score", score, perl = T)];

    # 将任意结尾带空格的字符串替换为空，只保留score后面的数字得分。
score <- gsub(".* ", "", as.character(score), perl = T);

    # 将字符向量转为数值向量。
score <- as.numeric(score);

    # 将score转换为n*n的数值矩阵，n为序列条数length(sequences[, 1])
scorem<- matrix(score, length(sequences[, 1]), length (sequences[, 1]),dimnames =list(as.vector(sequences[, 1]), as.vector(sequences[, 1])));

    # 得分矩阵求倒数，得到普通距离矩阵，用as.dist函数转换为下三角距离矩阵。
scorem.dist<- as.dist(1/scorem);

    # 根据距离矩阵，调用层次聚类函数hclust对所有序列聚类。
hc<- hclust(scorem.dist, method ="complete");

    # 绘制层次聚类结果。
    plot(hc, hang = -1, main ="Distance Tree Based on Needle All-Against-All Comparison", xlab =" sequence name", ylab ="distance");

    # 返回比对得分矩阵。
    scorem;
}
#____________________________________________________________________________









#_____________________________________________________________________________

##############################
#E2. 定义函数'infile_produce'#
##############################
infile_produce<- function(scorem) {
    # 求得分矩阵倒数，作为距离矩阵。
z <- 1/scorem;

    # 计算scorem行或列的长度
len = sqrt(length(scorem)) ;

    # 将距离矩阵中对角线赋值为0，seq生成一个从1到length(scorem)的向量，by为步长，向量中的值即对角线元素在z中的位置
z[seq(1, length(scorem), by = (len + 1))] <- 0;

    # 利用round函数将z中数值保留到小数点后7位
z <- round(z, 7) ;

    # 输出len长度信息到infile
write.table(len, file ="infile", quote = F, row.names = F, col.names = F) ;

    # 追加（append=T）输出距离矩阵z到infile
write.table(as.data.frame(z), file ="infile", append = T, quote = F, col.names = F, sep ="\t");

#结果提示
    cat("Phylip格式的距离矩阵已经输出到工作目录下名为'infile'的文件，以便使用Phylip软件继续进行分析。",  "\n");
}














################
#3.2.2 课题实现#
################
#_______________________________________________
#################################
#A. 调用函数'seq_import'导入数据#
################################
setwd("C:/workingdirectory");
my_file<- "AE004437.faa";
my_sequences<- seq_import(input_file = my_file);

#B. 调用函数'pattern_match'寻找模序
hit_sequences<- pattern_match(pattern ="H..H{1,2}", sequences = my_sequences, hit_num =2)
AA_percentage<- getAApercentage(sequences = hit_sequences)

#D. 调用函数'seq_alignment'进行序列两两比对
seq_alignment(sequences = hit_sequences)

#E1. 调用函数'getScoreMatrix'得到得分矩阵
score_matrix<- getScoreMatrix(sequences = hit_sequences)

#E2. 调用函数'infile_produce'生成PHYLIP软件的输入文件infile
infile_produce(scorem = score_matrix)
#_________________________________________________________







#############################################
#3.3 用R包（Bioconductor）再实现课题（方法一）
#############################################
#_________________________________________________________
#####################################
#A. 定义序列导入函数'bio_seq_import'#
#####################################
bio_seq_import<- function(input_file) {
      # 读入fasta文件，存入对象my_fasta
	my_fasta <- read.AAStringSet(input_file);

      #从my_fasta第一列的注释行中提取序列的ID(Accession Number)。
	Acc<- gsub(".*gb\\|(.*)\\|.*", "\\1", as.character(my_fasta[, 1]), perl = T);
     
      #修改my_fasta对象的names属性
	names(my_fasta)<-Acc;

	my_fasta;
}






#B. 定义模式匹配函数'bio_pattern_match'
bio_pattern_match<- function(pattern, sequences, hit_num) {
  # 从sequences对象中获取蛋白质序列的内容。
  seqs = as.character(sequences);

  # 获取正则表达式pattern表示的模序在所有序列中出现的位置（未找到匹配将返回
#-1），所有位置存入一个列表对象pos，perl=T表示兼容perl的正则表达式格式。
  pos<- gregexpr(pattern, seqs, perl = T);

  # lapply函数调用自定义函数function,根据pos中的每一个元素，计算
#pattern在每条序列中匹配的个数，再由unlist函数将结果转变为向量。
  hitsv<- unlist(lapply(pos, function(x) if (x[1] == -1) {0} else {length(x)})); 
 
  # 找出匹配次数大于hit_num的序列，并将大写形式替换为小写，gsub中第一个参数
#[A-Z]匹配任意大写字母，“\\L\\1”表示将前面小括号中匹配的任意字母替换为其小写形式。
  tag <- gsub("([A-Z])", "\\L\\1", as.character(sequences[hitsv > hit_num]), perl = T, ignore.case = T);

  # 为模序pattern加上小括号，以适合perl正则表达式格式，方便下面使用。
  pattern2 = paste("(", pattern, ")", sep ="");

  # 将tag序列中和模序pattern匹配的部分替换为大写，原理同上，“\\U\\1”表示替换为大写。
  tag<- gsub(pattern2, "\\U\\1", tag, perl = T, ignore.case = T);

  # 生成新的AAStringSet对象export保存选中的序列信息（大小写混合表示氨基酸）
  export <- AAStringSet(tag);

  # 输出对象export到文件Hit_sequences.fasta（fasta文件格式）
  write.XStringSet(export, file="Hit_sequences1.fasta");

  # 输出提示信息
  cat("含有模序\"", pattern, "\"超过", hit_num, "个的所有蛋白质序列已写入当前工作目录下文件'Hit_sequences.fasta'", "\n", sep ="");

  # 输出提示信息
  cat("极端嗜盐古菌蛋白组中以下序列含有模序\"", pattern, "\"的数量超过2个：", "\n", sep ="");

  # 输出选中的序列
  print(export);

  # 生成新的AAStringSet对象export保存选中的序列信息（大写表示氨基酸）
  selected <- AAStringSet(as.character(sequences[hitsv > hit_num]));

  # 返回选中序列
  selected;
}


#C. 定义氨基酸含量统计函数'bio_getAApercentage'#

bio_getAApercentage<- function(sequences) {
# 生成一个包含20种标准氨基酸单字母简写的数据框AA。
AA <- data.frame(AA =c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"));

# strsplit函数将序列内容sequences[, 6])转换成字符数组，lapply函数调用#table函数统计每条序列中各字符（氨基酸）出现的次数。
AAstat<- lapply(strsplit(as.character(sequences[, 6]), ""), table);
bio_getAApercentage<- function(sequences) {
      # 生成一个包含20种标准氨基酸单字母简写的数据框AA。
	AA <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y");
      # 得到每条序列中20种氨基酸出现的次数
	AApercentage <- letterFrequency(sequences, AA);

      # 次数/序列长度，得到每条序列中20种氨基酸的百分含量
      AApercentage <- t(AApercentage/width(sequences)* 100);
      
      # 修改数据框AApercentage各列的名字
	colnames(AApercentage) <- names(sequences);	
      
      # 下面代码与函数getAApercentage中相同，参看例3-1中的注释
	AApercentage <- data.frame(AApercentage, Mean =apply(AApercentage[, 1:dim(AApercentage)[2]], 1, mean, na.rm = T));
	write.csv(AApercentage, file ="AApercentage.csv", row.names = F, quote = F) ;
	cat("氨基酸百分比含量已经写入当前工作目录下的文件'AApercentage.csv'", "\n");
	AApercentage;
}

##########################################
#D. 定义两两比对和函数'bio_alignAndScore'#
##########################################
bio_alignAndScore<- function(sequences) {
      # 定义一个空的得分矩阵，初始值都为0
     scorem=matrix(rep(0,length(sequences)*length(sequences)),nrow = length(sequences), ncol = length(sequences));
      # 下面循环每次都是对比两条序列，结果存入得分矩阵scorem
	for (i in 1:length(sequences)) {
		for (j in 1:length(sequences)) {
                     #调用pairwiseAlignment函数两两对比
			scorem[i,j]=pairwiseAlignment(sequences[[i]], sequences[[j]], type = "overlap", substitutionMatrix = "BLOSUM62", gapOpening = 9.5, gapExtension = 0.5,scoreOnly=T)
		}
	}
	cat("程序完成所有序列的两两比对\n");
}
#_________________________________________________________




###############
#3.3.2 课题实现
###############
#____________________________________________________________

#首先，安装Bioconductor扩展包“Biostrings”：#
source("http://bioconductor.org/biocLite.R") ;
biocLite("Biostrings") ;
biocLite("ape") ;

#A. 利用函数'bio_seq_import'导入极端嗜盐古菌的蛋白质组#
library("Biostrings")
my_file <- "AE004437.faa" 
my_sequences <- bio_seq_import(input_file = my_file)

#B. 调用函数'bio_pattern_match'寻找模序
hit_sequences <- bio_pattern_match(pattern ="H..H{1,2}", sequences = my_sequences, hit_num =2)

#C. 调用函数' bio_getAApercentage'统计氨基酸百分含量
AA_percentage<- bio_getAApercentage(sequences = hit_sequences)

#D. 调用函数'alignAndScore'得到两两比对得分矩阵
score_matrix<- bio_alignAndScore(sequences = hit_sequences)
#____________________________________________________________________



######################################
#3.4 用R包（Bioconductor）再实现课题#
######################################
#_____________________________________________________
#E1. 定义函数'read.alignment'和'printMultipleAlignment'
printMultipleAlignment <- function(alignment, chunksize=60)
{
	#此函数需要Biostrings扩展包
	require("Biostrings") ;
	# 从对象alignment （由read.alignment读入）中得到序列总数
	numseqs <- alignment$nb ;
	# 得到所有序列比对在一起时的序列长度（=原始序列+indel）
	alignmentlen <- nchar(alignment$seq[[1]]) ;
	# 设定显示时，每个新行起始位置对应序列中的实际位置，如果设定每行最多不能超过#60bp（chunksize=60），则起始位置是1，61，121 …。这样，一个chunk包括了所有序列的一段 #不超过60bp的序列
	starts <- seq(1, alignmentlen, by=chunksize) ;
             # 得到chunk总数
	n <- length(starts) ;
	# 定义两个空向量
	aln <- vector();  #每个元素是一条序列
	lettersprinted <- vector(); #对应新行中最后一位在原始行中的位置
             # 这个循环完全可以用向量运算解决
	for (j in 1:numseqs)
	{
                         # 每次提取一条序列
		aln[j]  <- alignment$seq[[j]] ;        
                         # 初始化向量                  
		lettersprinted[j] <- 0 ;
	}
	# 每次循环处理一个chunk
	for (i in 1:n) { 
                          # 每次循环处理一条序列
		for (j in 1:numseqs)
		{
                                      # 每次取一条序列
			alnj <- aln[j] ;
                                      # 取该序列第i个chunk长度的片段
			chunkseqjaln <- substring(alnj, starts[i], starts[i]+chunksize-1) ;
                                       # 序列中的字母全部转换为大写
			chunkseqjaln <- toupper(chunkseqjaln) ;
			# 统计有多少“-”
			gapsj <- countPattern("-",chunkseqjaln) ;	
		             # chunk的长度减去“-”就是DNA字母数量，逐渐累加，即可对应到原始序#列中的位置
			lettersprinted[j] <- lettersprinted[j] + chunksize - gapsj;
			# 打印chunk中的片段，然后加上最后一位碱基在原始序列中的位置                         
			print(paste(chunkseqjaln,lettersprinted[j])) ;
		}
  			# 每个chunk打印完毕，需要空一行，分割多个chunk
print(paste(' ')) ;
	}
}



#E2. 定义函数'cleaned_aln'
cleanAlignment <- function(alignment, minpcnongap, minpcid)
{
     # 保留一份变量alignment的copy，用于存储更新后的信息，并作为返回值
     newalignment <- alignment;
     # 从对象alignment （由read.alignment读入）中得到序列总数
     numseqs <- alignment$nb;
     # 得到所有序列比对在一起时的序列长度（=原始序列+indel） 
     alignmentlen <- nchar(alignment$seq[[1]]) ;

     # 把newalignment对象中所有序列置空
     for (j in 1:numseqs) { newalignment$seq[[j]] <- "" };

     # 循环1开始，每次循环处理对齐序列中的一个位置
     for (i in 1:alignmentlen)
     {
        # 定义变量nongap记录该位置gap总数，并初始化:
        nongap <- 0;
        # 每次循环处理一条序列中的该位置所有对齐的残基
        for (j in 1:numseqs)
        {
          # 取第j条序列的所有残基
           seqj <- alignment$seq[[j]];
          # 只截取第j条序列第i个位置的残基
           letterij <- substr(seqj,i,i);
          # 如果出现不是“-”，nongap总数就加1
           if (letterij != "-") { nongap <- nongap + 1};
        }
         # 第i个位置的nongap总数除以序列数量，得到nongap的百分比
        pcnongap <- (nongap*100)/numseqs;
        # 如果某个位置的nongap 百分比含量大于等于阈值minpcnongap
        # 条件判断1开始
        if (pcnongap >= minpcnongap)
        {
           # 满足第一个条件，则还需要看是否满足第2个条件
           # 定义两个变量，第1个记录两两残基对总数，第2个记录相同残基对的数量，并初始化。
           numpairs <- 0; numid <- 0;
           # 下面两重循环用于第i个位置上所有残基，两两比较
           for (j in 1:(numseqs-1))
           {
               # 只截取第j条序列第i个位置的残基
              seqj <- alignment$seq[[j]];
               # 只截取第j条序列第i个位置的残基
              letterij <- substr(seqj,i,i);
               # 再取第k条序列第i个位置的残基，逐个与第j条i位置比较
              for (k in (j+1):numseqs)
              {
                 seqk <- alignment$seq[[k]];
                 letterkj <- substr(seqk,i,i);
               # 再取第k条和j条都不是gap的
                 if (letterij != "-" && letterkj != "-")
                 {
               # 两两残基对总数计数增加1次
                    numpairs <- numpairs + 1;
               # 如果这对残基相同，相同残基对的计数加1
                    if (letterij == letterkj) { numid <- numid + 1};
                 }
              }
           }
           # 相同残基对除以两两残基对总数，得到比例
           pcid <- (numid*100)/(numpairs);
           # 条件判断2开始，如果上面的比例大于阈值
           if (pcid >= minpcid)
           {
               # 就把第i位置上，所有序列相应的残基依次写入，否则就丢弃该位点
               for (j in 1:numseqs)
               {
                  seqj <- alignment$seq[[j]];
                  letterij <- substr(seqj,i,i);
                  newalignmentj <- newalignment$seq[[j]];
                  newalignmentj <- paste(newalignmentj,letterij,sep="");
                  newalignment$seq[[j]] <- newalignmentj;
               }
           } # 条件判断2结束
        } # 条件判断1结束
     }  # 循环1结束
     return(newalignment);
}


#E3. 定义函数'unrootedNJtree'
unrootedNJtree <- function(alignment,type)
{
     # 这个函数需要 ape 和 seqinR 扩展包:
     require("ape");
     require("seqinr");
     # 定义一个函数，注意这个新知识点，是函数内定义函数
     makemytree <- function(alignmentmat)
     {
     # as开头的函数都是格式转换
        alignment <- ape::as.alignment(alignmentmat);
     # 如果序列类型是蛋白质
        if      (type == "protein")
        {
     # 从比对结果对象中得到一个两两距离矩阵
           mydist <- dist.alignment(alignment);
        }
     # 如果序列类型是DNA
        else if (type == "DNA")
        {
     # as开头的函数都用于格式转换
           alignmentbin <- as.DNAbin(alignment);
     # 从比对结果对象中得到一个两两距离矩阵
           mydist <- dist.dna(alignmentbin);
        }
     # 用邻位相连法（Neighbor-joining）构建进化树对象
        mytree <- nj(mydist);
     # 将构建的进化树对象返回
        return(mytree);
     }
     # as开头的函数都用于格式转换
     mymat  <- as.matrix.alignment(alignment);
     # 调用上面定义的函数构建进化树
     mytree <- makemytree(mymat);
     # 对构建的进化树做自举（bootstrap）分析
     myboot <- boot.phylo(mytree, mymat, makemytree);
     # 画进化树，类型是无根树
     plot.phylo(mytree,type="u");
     # 在画好的进化树上显示自举值
     nodelabels(myboot,cex=0.7);   
     # 把自举值设定为节点的标签
     mytree$node.label <- myboot;  
     # 返回构建的进化树（对象）
     return(mytree);
}


#E4. 定义函数'rootedNJtree'
rootedNJtree <- function (alignment, theoutgroup, type)
{
     # 本函数大部分代码与unrootedNJtree函数相同，因此只注释不同的语句
     require("ape");
     require("seqinr");
     # 定义一个函数，多了一个参数，用于指定哪条序列来充当根
     makemytree <- function (alignmentmat, outgroup=`theoutgroup`)
     {
        alignment <- ape::as.alignment(alignmentmat);
        if      (type == "protein")
        {
           mydist <- dist.alignment(alignment);
        }
        else if (type == "DNA")
        {
           alignmentbin <- as.DNAbin(alignment);
           mydist <- dist.dna(alignmentbin);
        }
        mytree <- nj(mydist);
     # 这里需要指明哪条序列作为根
        myrootedtree <- root(mytree, outgroup, r=TRUE);
        return(myrootedtree);
     }
     mymat  <- as.matrix.alignment(alignment);
     # 这里函数调用，也多了一个参数
     myrootedtree <- makemytree(mymat, outgroup=theoutgroup);     
     myboot <- boot.phylo(myrootedtree, mymat, makemytree);
     # 画进化树，类型是有根树
     plot.phylo(myrootedtree,type="p");
     nodelabels(myboot,cex=0.7);   
     myrootedtree$node.label <- myboot;  
     return(myrootedtree);
}
#____________________________________________________________________________


################
#3.4.2 课题实现#
################
#___________________________
install.packages('seqinr');
install.packages('ape');

#E1. 读入多重对比结果并显示，供人工查看
library('seqinr')
hit_aln  <- read.alignment(file = "Hit_sequences.clustalw", format = "clustal");
printMultipleAlignment(hit_aln, 60)


#E2. 去除多重比对中的低质量区
cleaned_aln <- cleanAlignment(hit_aln, 25, 25)

#E3. 根据多重比对结果绘制无根树
install.packages('ape')
library('ape')
unrooted_tree <- unrootedNJtree(cleaned_aln,"protein")

#E4. 根据多重比对结果绘制有根树
rooted_tree <- rootedNJtree(cleaned_aln, "AAG19157.1", "protein")
#_____________________________________________________________________


