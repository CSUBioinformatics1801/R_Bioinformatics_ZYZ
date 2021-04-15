
# R软件更新
install.packages("installr")
library(installr)
updateR()

# 查看程序包对应哪个版本的R编译的
pkgs<-installed.packages()
plot(as.factor(pkgs[,'Built']),col=2:4,main='Packages built version',ylab='Count of packages')

# 检查包库路径
.libPaths()

# 获取旧包名称
old_packages <- installed.packages(lib.loc = "/home/gzucm04/R/x86_64-pc-linux-gnu-library/3.6")
old_packages <- as.data.frame(old_packages)
list.of.packages <- unlist(old_packages$Package)

# 删除旧R包
remove.packages( installed.packages( priority = "NA" )[,1] )

# 重新安装所有R程序包
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages,function(x){library(x,character.only=TRUE)})