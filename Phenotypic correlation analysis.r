#Analysis by traits
setwd("D:\\phenotype")
data1 <- read.table("AP3.0.txt", header = TRUE)
install.packages("GGally")  
library(GGally)
library(ggplot2)
my_cor <- function(data, mapping, ...) {
  ggally_cor(data, mapping, ..., digits = 2)  
}
ggpairs(data1,
        columns = 2:4,
        aes(color = group),
        upper = list(continuous = wrap("cor", size = 5)),
        lower = list(continuous = wrap("points", alpha = 0.6)),
        diag = list(continuous = wrap("densityDiag", alpha = 0.4)))  
#Analysis by source
#Hebei(DC)
setwd(dir="D:\\phenotype")
data1=read.table("hb.txt",head=T)
data1=data1[,2:4]
library("Hmisc")
test1 <- rcorr(as.matrix(data1))
test1=rcorr(as.matrix(data1))
head(test1)
library("PerformanceAnalytics")
chart.Correlation(data1,histogram=TRUE)
#Heilongjiang(SW)
setwd(dir="D:\\phenotype")
data1=read.table("hlj.txt",head=T)
data1=data1[,1:3]
library("Hmisc")
test1 <- rcorr(as.matrix(data1))
test1=rcorr(as.matrix(data1))
head(test1)
library("PerformanceAnalytics")
chart.Correlation(data1,histogram=TRUE)
#qingyang
setwd(dir="D:\\phenotype")
data1=read.table("qy.txt",head=T)
data1=data1[,1:3]
library("Hmisc")
test1 <- rcorr(as.matrix(data1))
test1=rcorr(as.matrix(data1))
head(test1)
library("PerformanceAnalytics")
chart.Correlation(data1,histogram=TRUE)

