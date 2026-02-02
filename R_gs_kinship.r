setwd("D:/gs2")
myY  <- read.table("AP.txt", head = TRUE)
myGD <- read.csv("myGD.csv", head = TRUE)
myGM <- read.csv("myGM.csv", head = TRUE)
myGD_data <- myGD[, -1]
colnames(myGD)[1] <- "Taxa"
myY <- myY[!myY$Taxa %in% remove_id, ]###+
dim(myY)
remove_id1<- c("D00592537", "D00628064")
myGD<- myGD[!myGD$Taxa %in% remove_id1, ]
dim(myGD)
source("D:/gs2/10.16/gapit_functions (6).txt")
metal_path <- "D:/gs2/generic-metal/metal.exe"

nrep <- 1
nfold <- 5
cutOff <- 0.01
nSNP <- nrow(myGM)
sig_thre <- cutOff / nSNP


gd_mat <- as.matrix(myGD[,-1])
gd_var <- apply(gd_mat, 2, var, na.rm = TRUE)
gd_mat <- gd_mat[, gd_var > 0, drop = FALSE]
pca_res <- prcomp(gd_mat, scale. = TRUE)
PCs <- pca_res$x[, 1:3]
colnames(PCs) <- paste0("PC", 1:3)

results_all <- data.frame(
  trait = character(),
  rep = integer(),
  fold = integer(),
  method = character(),
  cor_pred = numeric(),
  cor_gBV  = numeric(),
  QTNs = integer(),
  stringsAsFactors = FALSE
)
traits <- c("DP")
set.seed(1000)
for (trait in traits) {
  cat("========== 正在预测性状:", trait, "==========\n")
  for (i in 1:nrep) {
      n_total <- nrow(myY)
      group_labels <- sample(rep(1:3, length.out = n_total))
      myY$Group <- group_labels
     
      population1 <- subset(myY, Group == 1)[, c("Taxa", trait)]
      population2 <- subset(myY, Group == 2)[, c("Taxa", trait)]
      population3 <- subset(myY, Group == 3)[, c("Taxa", trait)]
    
     cat(sprintf("rep=%d 随机分群完成：population1=%d population2=%d population3=%d\n",
                i, nrow(population1), nrow(population2), nrow(population3)))
    
     myCV <- data.frame(
      Taxa = myY$Taxa,
      population1 = ifelse(myY$Group == 1, 1, 0),
      population2 = ifelse(myY$Group == 2, 1, 0),
      population3 = ifelse(myY$Group == 3, 1, 0)
    )
     
     PCs1<- data.frame(Taxa = myGD$Taxa, PCs)
     myCV1<- merge(myCV[, 1:3], PCs1, by = "Taxa", all.x = TRUE)##blink
	 myCV2 <- myCV1[, -c(2, 3)]##meta
	 myY_pop1 <- population1
     myY_pop2 <- population2
     myY_pop3 <- population3
     for (j in 1:nfold) {
      cat(sprintf("rep=%d fold=%d trait=%s\n", i, j, trait))

      sets_pop1 <- sample(cut(1:nrow(myY_pop1), nfold, labels = FALSE), nrow(myY_pop1))
      sets_pop2 <- sample(cut(1:nrow(myY_pop2), nfold, labels = FALSE), nrow(myY_pop2))
      sets_pop3 <- sample(cut(1:nrow(myY_pop3), nfold, labels = FALSE), nrow(myY_pop3))

   
      train_pop1 <- myY_pop1[sets_pop1 != j, c("Taxa", trait)]
      test_pop1  <- myY_pop1[sets_pop1 == j, c("Taxa", trait)]
      train_pop2 <- myY_pop2[sets_pop2 != j, c("Taxa", trait)]
      test_pop2  <- myY_pop2[sets_pop2 == j, c("Taxa", trait)]
      train_pop3 <- myY_pop3[sets_pop3 != j, c("Taxa", trait)]
      test_pop3  <- myY_pop3[sets_pop3 == j, c("Taxa", trait)]
      
      names(train_pop1)[2] <- "Phenotype"
      names(train_pop2)[2] <- "Phenotype"
      names(train_pop3)[2] <- "Phenotype"
      names(test_pop1)[2]  <- "Phenotype"
      names(test_pop2)[2]  <- "Phenotype"
      names(test_pop3)[2]  <- "Phenotype"
    
      train_all <- rbind(train_pop1, train_pop2, train_pop3)
      colnames(train_all) <- c("Taxa", "Phenotype")
      myGAPIT1 <- GAPIT(
        Y = train_all,
        GD = myGD,
        GM = myGM,
        CV = myCV1,
        cutOff = cutOff,
        Random.model = FALSE,
        model = "Blink", 
		Multi_iter = FALSE
      )
      if (!is.null(myGAPIT1$GWAS) && "P.value" %in% names(myGAPIT1$GWAS)) {
      sig_res <- subset(myGAPIT1$GWAS, P.value < sig_thre)
      blink_sig_count <- nrow(sig_res)
      blink_sig_snps  <- intersect(sig_res$SNP, colnames(myGD))
}     else {
      blink_sig_count <- 0
      blink_sig_snps  <- character(0)
}

      cat("BLINK显著位点个数为:", blink_sig_count, "\n")

  	  if (length(blink_sig_count) > 0) {
      myGD_sig_df <- data.frame(Taxa = myGD$Taxa, myGD[, blink_sig_snps, drop=FALSE])
      myCV3<- merge(myCV1, myGD_sig_df, by = "Taxa", all.x = TRUE)
}     else {
      myCV3 <- myCV1
} 
     
      blink_gs <- GAPIT(
       Y = train_all,
       GD = myGD,
       GM = myGM,
       CV = myCV3,
	   CV.Extragenetic = 2,
       model=c("gBLUP")
)
      pred_blink <- blink_gs$Pred
	  
	  test_all   <- rbind(test_pop1,test_pop2,test_pop3)
      colnames(test_all) <- c("Taxa", "Phenotype")
	  
      merged_blink <- merge(test_all, pred_blink[, c("Taxa", "BLUE.N", "gBreedingValue")],
                           by = "Taxa", all.x = TRUE)
						   
						   merged_blink <- merge(test_all, pred_blink,
                           by = "Taxa", all.x = TRUE)
      cor_pred_blink <- cor(merged_blink$Phenotype, merged_blink$BLUE.N)
      cor_gBV_blink<- cor(merged_blink$Phenotype, merged_blink$gBreedingValue)
      
	  cor_pred_blink
	  cor_gBV_blink
	  
     results_all <- rbind(results_all,
       data.frame(trait=trait, rep=i, fold=j, method="BLINK",
             cor_pred=cor_blue.N_blink, cor_gBV=cor_gBV_blink,
             QTNs = blink_sig_count,
             stringsAsFactors = FALSE))
      # =================== (2)meta+GS ===================
	  run_gwas <- function(train_df){
        res <- GAPIT(
          Y = train_df,
          GD = myGD,
          GM = myGM,
          CV = myCV2,
          cutOff = cutOff,
          model = "Blink",
		  Multi_iter = FALSE
        )
        return(res$GWAS)
      }

      gwas_pop1 <- run_gwas(train_pop1)
      gwas_pop2 <- run_gwas(train_pop2)
      gwas_pop3 <- run_gwas(train_pop3)

      fn_pop1  <- paste0("pop1_rep", i, "_fold", j, ".txt")
      fn_pop2  <- paste0("pop2_rep", i, "_fold", j, ".txt")
      fn_pop3 <- paste0("pop3_rep", i, "_fold", j, ".txt")

      write.table(gwas_pop1[,c("SNP","effect","P.value","nobs")],  fn_pop1,  sep="\t", row.names=FALSE, quote=FALSE)
      write.table(gwas_pop2[,c("SNP","effect","P.value","nobs")],  fn_pop2,  sep="\t", row.names=FALSE, quote=FALSE)
      write.table(gwas_pop3[,c("SNP","effect","P.value","nobs")],  fn_pop3, sep="\t", row.names=FALSE, quote=FALSE)
      

      metal_script <- paste0("metal_rep", i, "_fold", j, ".txt")
      meta_out <- paste0("meta_rep", i, "_fold", j)

      ms <- c(
        "SCHEME SAMPLESIZE",
        "MARKER SNP",
        "EFFECT effect",
        "PVALUE P.value",
        "WEIGHT nobs",
        paste0("PROCESS ", fn_pop1),
        paste0("PROCESS ", fn_pop2),
        paste0("PROCESS ", fn_pop3),
        paste0("OUTFILE ", meta_out),
        "ANALYZE",
        "QUIT"
      )
      writeLines(ms, metal_script)
      
	 system(paste(metal_path, metal_script), ignore.stdout=FALSE, ignore.stderr=FALSE)
     #taxa  <- myGD$taxa
     meta_tbl <- "METAANALYSIS1.TBL"

     if(!file.exists(meta_tbl)) stop("METAANALYSIS1.TBL 不存在！")

     meta_res <- read.table(meta_tbl, header=TRUE, sep="\t", stringsAsFactors=FALSE)
     sig_res  <- subset(meta_res, P.value < sig_thre)
     cat("meta显著位点数量:", nrow(sig_res), "\n")
     meta_sig_count <- nrow(sig_res)
     sig_snps <- intersect(sig_res$MarkerName, colnames(myGD_data))
 
     if (length(sig_snps) > 0) {
         myGD_sig_df <- data.frame(Taxa = myGD$Taxa, myGD_data[, sig_snps, drop = FALSE])
         myCV4 <- merge(myCV2, myGD_sig_df, by = "Taxa", all.x = TRUE)
        }else {
         myCV4 <- myCV2
} 
     meta_gs <- GAPIT(
       Y = train_all,
       GD = myGD,
       GM = myGM,
       CV = myCV4,
      model=c("gBLUP")
)
   
      pred_meta <- meta_gs$Pred
     
      merged_meta <- merge(
	                      test_all, 
						  pred_meta[, c("Taxa", "BLUE.N", "gBreedingValue")],
                          by = "Taxa", 
						  all.x = TRUE)
    
      cor_pred_meta <- cor(merged_meta$Phenotype, merged_meta$BLUE.N)
      cor_gBV_meta  <- cor(merged_meta$Phenotype, merged_meta$gBreedingValue)
      
      results_all <- rbind(results_all,
                           data.frame(trait=trait, rep=i, fold=j, method="Meta",
                                      cor_pred=cor_pred_meta, cor_gBV=cor_gBV_meta,
                                      QTNs = meta_sig_count, stringsAsFactors = FALSE))
      
      cat(sprintf("Trait=%s | rep=%d | fold=%d | BLINK cor(Pred)=%.4f cor(gBV)=%.4f | Meta cor(Pred)=%.4f cor(gBV)=%.4f\n",
                  trait, i, j,
                  ifelse(is.na(cor_pred_blink), NaN, cor_pred_blink),
                  ifelse(is.na(cor_gBV_blink), NaN, cor_gBV_blink),
                  ifelse(is.na(cor_pred_meta), NaN, cor_pred_meta),
                  ifelse(is.na(cor_gBV_meta), NaN, cor_gBV_meta)))
}}}
   
if (!dir.exists("D:/gs results")) dir.create("D:/gs results", recursive = TRUE)
write.csv(results_all, "D:/gs results/5.0DP_blink_vs_meta_randomGroups_kinship.csv", row.names = FALSE)

