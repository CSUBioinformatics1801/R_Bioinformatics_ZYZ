#######################
## 安装本章需要的R包 ##
#######################
#已安装的R包
pkgs <- rownames(installed.packages())
#本章所需包
packages_1=c("ggmap","maps","XML","ggplot2","maproj","igraph","quantmod","scales")
#将要安装的R包
packages_2=setdiff(packages_1,pkgs) 
#安装所需R包 
install.packages(packages_2)


##简单的语法知识##
#赋值变量
x <- c(1,2,3,4)
x

#加法运算
c(1,2,3,4) + c(3,4,5,6)

#字符向量
c("hello world", "I am a R user")

#调用函数exp计算0的自然指数
exp(0)

#产生一个向量 (1,2,3,4)
x<-1:4
#计算x中每个值的指数
exp(x)

#计算圆面积的函数，计算半径为4的圆面积
area<-function(x){  
 result<-pi*x^2
 return(result)
}
area(4)

#计算两个数值的和与差
add.diff<-function(x,y){
 add<-x+y
 diff<-x-y
 return(c(add,diff))
}
#调用函数add.diff来计算5和3的和与差
add.diff(x=5,y=3)  


#############
##  例1-1  ##
#############

library(ggplot2)  #加载ggplot2包
library(ggmap)  #加载ggmap包
library(XML)   #加载XML包
library(maps)    #加载maps包
library(mapproj)  #加载mapproj包

# 将数据所在网址作为字符串存入url变量中
url <- 'http://data.earthquake.cn/datashare/globeEarthquake_csn.html';
#windows下可能要加下面这句解码
url= htmlParse(url,encoding="UTF-8")
# 对该网页内容进行解析，读取其中的所有表格，并存入tables变量中
tables <- readHTMLTable(url,stringsAsFactors = FALSE); 
# 取出我们所需的第6个表格，存入变量raw。
raw <- tables[[6]] ;  
#查看表格的第一行数据，这里不显示结果。
raw[1,]
# 只保留时间、经度、纬度这三列数据，并存入变量data。
data <- raw[ ,c(1,3,4)] ;  
# 修改data包含的表格各列的名称为'date'、'lan'和'lon'。
names(data) <- c('date','lan','lon') ;
# 将经度(data$lan)和纬度(data$lon)的数据类型用函数as.numeric（）转换为数值类型
data$lan <- as.numeric(data$lan) ;
data$lon <- as.numeric(data$lon) ;
# 将时间(data$date)的数据类型用函数as.Date（）转换为时间类型("%Y-%m-%d"))
data$date <- as.Date(data$date, "%Y-%m-%d");
# 用ggmap包读取地图（该地图中心为 'china'，放大4倍，地图类型为地形图，范围为整个图形设备）+ 叠加散点图（散点数据来源于data数据框，以数据经纬度作为坐标值，散点颜色为红色，透明度为0.7，图解位置无）。
ggmap(get_googlemap(center='china', zoom=4, maptype='terrain'), extend='device')+
  geom_point (data=data, aes(x=lon,y=lan), colour='red', alpha=0.4)+
  opts(legend.position="none");

#############
##  例1-2  ##
#############

#加载igraph包
library (igraph);
#根据Barabasi-Albert模型生成一个网络，该网络包括100个节点，每次生成时出现一条边 #（m=1）。变量g是一个对象（参见第七章），包含了这个网络的所有信息
g <- barabasi.game (100, m=1);
plot( g,                        #画图的对象是g
 vertex.size=4, # 顶点的大小设置为4
 vertex.label=NA, # 顶点的标签设置为无
 edge.arrow.size=0.1, # 点之间连线的箭头大小设置为0.1
 edge.color="grey40", # 点之间连线的颜色设置为灰度
 layout=layout.fruchterman.reingold, # 设置整体的布局方式
 vertex.color="red", # 顶点的颜色设置为红色
 frame=TRUE);   #绘图包括整体边框

#############
##  例1-3  ##
#############

#加载quantmod包
library(quantmod)  
# getSymbols函数自动连接Yahoo数据源，必须保证网络连接正常。读取的数据源是SSEC #（上证指数的Yahoo代码），时间范围是从 '2011-01-01' 到 '2012-07-13'
getSymbols('^SSEC', from = '2011-01-01', to='2012-07-13')
#得到的数据保存在SSEC这个对象中，查看表格的第一行数据，这里不显示结果
SSEC[1,]
#绘图K线图，时间是最近的4个月，主题使用白色蜡烛图，不使用其它的技术指标(TA)
candleChart (last(SSEC,'4 months'), theme=chartTheme('white'), TA=NULL)
# 添加MACD指标，使用默认参数
addMACD()

#############
##  例1-4  ##
#############

library(ggplot2);  #加载ggplot2包。
library(scales);   #加载scales包。
revigo.names <- c("term_ID", "description", "frequency_%", "plot_X", "plot_Y", "plot_size", "log10_p_value", "uniqueness", "dispensability");
revigo.data <- rbind(c("GO:0009698", "phenylpropanoid metabolic process",  0.024,  1.380, -8.143,  3.438, -6.7645, 0.675, 0.000), 
c("GO:0021722", "superior olivary nucleus maturation",  0.000, -4.773,  3.914,  0.845, -3.5935, 0.874, 0.000), 
c("GO:0019748", "secondary metabolic process",  0.081, -6.414,  0.228,  3.966, -5.0550, 0.909, 0.007), 
c("GO:0042440", "pigment metabolic process",  0.324, -0.765,  0.844,  4.568, -4.3893, 0.909, 0.009), 
c("GO:0051299", "centrosome separation",  0.001,  0.901, -5.814,  1.771, -3.0482, 0.902, 0.012), 
c("GO:0052695", "cellular glucuronidation",  0.000,  6.070, -2.020,  1.398, -6.7645, 0.667, 0.023), 
c("GO:0042501", "serine phosphorylation of STAT protein",  0.001,  4.269,  4.220,  1.944, -3.5935, 0.760, 0.025), 
c("GO:0016101", "diterpenoid metabolic process",  0.005, -4.590, -4.132,  2.772, -3.5346, 0.763, 0.028), 
c("GO:0090313", "regulation of protein targeting to membrane",  0.000,  1.015,  6.202,  1.230, -3.5935, 0.822, 0.133), 
c("GO:0010225", "response to UV-C",  0.001,  6.233,  3.496,  2.121, -3.1409, 0.864, 0.171), 
c("GO:0006063", "uronic acid metabolic process",  0.025,  6.023, -2.786,  3.461, -5.4473, 0.697, 0.408), 
c("GO:0006069", "ethanol oxidation",  0.015,  5.585, -3.328,  3.222, -3.1555, 0.762, 0.427), 
c("GO:0006720", "isoprenoid metabolic process",  0.401, -3.931, -4.897,  4.660, -3.1707, 0.775, 0.481), 
c("GO:0021819", "layer formation in cerebral cortex",  0.000, -4.381,  4.333,  1.740, -3.5935, 0.845, 0.559), 
c("GO:0001523", "retinoid metabolic process",  0.003, -4.238, -4.361,  2.480, -3.6819, 0.741, 0.595), 
c("GO:0090314", "positive regulation of protein targeting to membrane",  0.000,  0.506,  6.348,  0.301, -3.5935, 0.827, 0.598),
c("GO:0052697", "xenobiotic glucuronidation",  0.000,  5.854, -0.744,  1.114, -6.9626, 0.570, 0.617));
#以下都是数据格式的转换
one.data <- data.frame(revigo.data);  #矩阵格式转换为数据框格式
names(one.data) <- revigo.names;    #改变各列名称
#只保留x、y坐标都不为null，即有数值的行
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );#factor类型转字符，再转数字
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) ); #同上
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) ); #同上
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) ); #同上
one.data$frequency <- as.numeric( as.character(one.data$frequency) ); #同上
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) ); #同上
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) ); #同上
#以下使用ggplot绘图
p1 <- ggplot( data = one.data );    #建立基本绘图对象,将绘图所需数据传递给该对象。
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_area();   #确定X轴、Y轴、颜色、大小和透明度的映射规则
p1 <- p1 + scale_colour_gradientn (colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );  #确定颜色过度的标度控制
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_area();      #添加点图对象，并设置颜色
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw();    #设置点的大小标度
ex <- one.data [ one.data$dispensability < 0.15, ];    #数据取子集（只要最后一列小于0.15的）
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );            #添加文字对象并设置颜色和大小
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");    #添加X轴和Y轴说明
p1 <- p1 + opts(legend.key = theme_blank()) ;   #添加图例说明
#以下是设置X轴和Y轴的刻度限
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);
p1; #执行绘图指令

