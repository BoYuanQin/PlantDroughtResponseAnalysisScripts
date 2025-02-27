# GENE_PATH_NEXUS

## 简介
PPI和WGCNA分析。

## 指南
1、使用conda安装依赖
2、ProNet需要本地编译，flashClust需要使用R安装
install.packages("ProNet_1.0.0.tar.gz", repos=NULL, type="source")
install.packages("flashClust")
3、把data物种数据放入input,下载STRING数据库放入STRINGdb（aliases\info\links三个文件）
4、运行main.r
注：原始数据要去除NA，并且进行ID转换，满足stringdb数据输入要求
