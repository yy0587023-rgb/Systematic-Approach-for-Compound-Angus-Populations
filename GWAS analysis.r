#GWAS分析
setwd("D:/gwas")
file_path ="gapit_functions.txt"
source(file_path)
install.packages("Biobase_2.66.0.zip", repos = NULL, type = "win.binary")
library(Biobase)

myY=read.table("AP.txt", head = TRUE)
myG=read.table("mdp.genotype.hmp.txt", head = FALSE)
myCV=read.csv("CV.csv", head = TRUE) 
myCV=myCV[,1:3]
myGAPIT <- GAPIT(
 G=myG,
 PCA.total=3,
 )
myGD= myGAPIT$GD
myGM= myGAPIT$GM
write.csv(myGD,"D:/gwas/myGD.csv",quote=F,row.names=F)
write.csv(myGM,"D:/gwas/myGM.csv",quote=F,row.names=F)

#shut off 方法降低阈值线
setwd("D:/gwas/降低阈值线")
myGD=read.csv("myGD.csv", head = TRUE) 
myGM=read.csv("myGM.csv", head = TRUE) 
myY=read.table("AP.txt", head = TRUE)
myCV=read.csv("CV.csv", head = TRUE) 
myCV=myCV[,1:3]
Y <- myY[, c("Taxa", "CW")]  # 选择 "Taxa" 和 "CW" 列作为表型数据
n_perm <-50               # 设置进行n次 permutation
#存储每次 permutation 得到的最小 p 值
permuted_pvals <- numeric(n_perm)
#进行n次 permutation
for (i in 1:n_perm) {
#打乱表型数据的 "CW" 列
  Y_permuted <- Y
  Y_permuted[, 2] <- sample(Y[, 2])  # 打乱CW列的表型数据
#用 permutation 后的表型数据进行 GWAS 分析
  myGAPIT_perm <- GAPIT(
    Y = Y_permuted,
    GM = myGM,
    GD = myGD,
    CV = myCV,
    PCA.total = 3,
    model = c("BLINK"),
	)
#保存每次 permutation 得到的最小p值
  permuted_pvals[i] <- min(myGAPIT_perm$GWAS$P.value)
}
#循环完成后一次性打印所有最小 p 值
cat("All permutation min p-values:\n")
print(permuted_pvals)
#计算 permutation-based 的显著性阈值（分位数法）
alpha <- 0.05
threshold <- quantile(permuted_pvals, probs = alpha)
#输出新的阈值
cat("Permutation-based shut-off threshold (50 times) for CW:", threshold, "\n")

#用新的阈值线对全部个体做GWAS
myGAPIT=GAPIT(
    Y =myY[,c(1,3)] ,
    GM = myGM,
    GD = myGD,
    CV = myCV,
    PCA.total = 3,
    model = c("BLINK"),
	cutOff=0.0861)
#用新的阈值线对单独的三个群体做GWAS
myY=read.table("hb.txt", head = TRUE)
myCV=read.csv("CV.csv", head = TRUE) 
myGD=read.csv("myGD.csv", head = TRUE) 
myGM=read.csv("myGM.csv", head = TRUE) 
myGAPIT=GAPIT(
  Y = myY[,c(1,3)],
  GM=myGM,
  GD=myGD,
  PCA.total = 3,
  model = c("BLINK"),
  cutOff=0.0861
 )
 
myY=read.table("qy.txt", head = TRUE)
myGAPIT=GAPIT(
  Y = myY[,c(1,5)],
  GM=myGM,
  GD=myGD, 
  PCA.total = 3,
  model = c("BLINK"),
  cutOff=0.0861
 )

myY=read.table("hlj.txt", head = TRUE)
myGAPIT=GAPIT(
  Y = myY[,c(1,3)],
  GM=myGM,
  GD=myGD,
  PCA.total = 3,
  model = c("BLINK"),
  cutOff=0.0861
 )
