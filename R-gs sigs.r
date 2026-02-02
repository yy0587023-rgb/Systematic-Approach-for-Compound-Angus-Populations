# 随机分群significants 
setwd("D:/gs/significant sites method")
myY  <- read.table("AP.txt", head = TRUE)
myGD <- read.csv("myGD.csv", head = TRUE)
myGM <- read.csv("myGM.csv", head = TRUE)
myCV <- read.table("CV.txt", head = TRUE)[, 1:3]
source("D:/gs/significant sites method/gapit_functions (6).txt")
metal_path <- "D:/gs/generic-metal/metal.exe"
remove_id<- c("D00592531", "D00608064")
myY <- myY[!myY$Taxa %in% remove_id, ]
dim(myY)
remove_id1<- c("D00592537", "D00628064")
myGD<- myGD[!myGD$Taxa %in% remove_id1, ]
colnames(myGD)[1] <- "Taxa"
dim(myGD)
myCV<- myCV[!myCV$Taxa %in% remove_id, ]
dim(myCV)
nrep <- 20
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

traits <- c("CW")
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
     myCV1<- merge(myCV[, 1:3], PCs1, by = "Taxa", all.x = TRUE)
	 
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

      test_all <- rbind(test_pop1,test_pop2,test_pop3)
	  
      train_df1 <- merge(train_all, myCV1[, setdiff(names(myCV1), "Phenotype")], by = "Taxa")
      train_df1$Phenotype <- NULL  
      test_df1 <- merge(test_all, myCV1[, setdiff(names(myCV1), "Phenotype")], by = "Taxa")
      test_df1$Phenotype <- NULL  

  	  if (length(blink_sig_snps) > 0) {
      myGD_sig_df <- data.frame(Taxa = myGD$Taxa, myGD[, blink_sig_snps, drop=FALSE])
      train_df1 <- merge(train_df1, myGD_sig_df, by = "Taxa", all.x = TRUE)
      test_df1  <- merge(test_df1,  myGD_sig_df, by = "Taxa", all.x = TRUE)
}
      train_df1$taxa <- NULL  
      test_df1$taxa  <- NULL  
	  
      train_pheno <- train_all$Phenotype
      lm_model <- lm(train_pheno ~ ., data = train_df1[, -1])
      pred_values <- predict(lm_model, newdata = test_df1[, -1])

      cor_pred_blink <- cor(test_all$Phenotype, pred_values, use = "complete.obs")
     
      cat("BLINK预测与真实表型的相关性为:", cor_pred_blink, "\n")
####cor_gBV_blink 	 
      train_df2 <- data.frame(Taxa = train_all$Taxa)
      train_df2 <- merge(train_df2, PCs1, by = "Taxa", all.x = TRUE)
     
      test_df2 <- data.frame(Taxa = test_all$Taxa)
      test_df2 <- merge(test_df2, PCs1, by = "Taxa", all.x = TRUE)
     
      train_pheno <- train_all$Phenotype
      lm_model <- lm(train_pheno ~ ., data = train_df2[, -1])
      pred_values1 <- predict(lm_model, newdata = test_df2[, -1])

      cor_gBV_blink  <- cor(test_all$Phenotype, pred_values1, use = "complete.obs")
   
      cat("BLINK_gbv与真实表型的相关性为:",  cor_gBV_blink , "\n")
	  
      results_all <- rbind(
        results_all,
        data.frame(
          trait = trait,
          rep = i,
          fold = j,
          method = "BLINK",
          cor_pred = cor_pred_blink,
          cor_gBV = cor_gBV_blink,
          QTNs = blink_sig_count,
          stringsAsFactors = FALSE
        )
      )
      
#Meta 
      run_gwas <- function(train_df) {
        res <- GAPIT(
          Y = train_df,
          GD = myGD,
          GM = myGM,
          CV = myCV1,
          cutOff = cutOff,
		  Random.model = FALSE,
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
      
      system(paste(metal_path, metal_script), ignore.stdout = FALSE, ignore.stderr = FALSE)
      
      myGD_data <- myGD[, -1]
      meta_tbl <- "METAANALYSIS1.TBL"
      
      if (!file.exists(meta_tbl)) stop("METAANALYSIS1.TBL 不存在！")
      
      meta_res <- read.table(meta_tbl, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      sig_res  <- subset(meta_res, P.value < sig_thre)
      meta_sig_count <- nrow(sig_res)
      
      sig_snps <- intersect(sig_res$MarkerName, colnames(myGD_data))

      train_df3 <- data.frame(Taxa = train_all$Taxa)
      train_df3 <- merge(train_df3, PCs1, by = "Taxa", all.x = TRUE)
     
      test_df3 <- data.frame(Taxa = test_all$Taxa)
      test_df3 <- merge(test_df3, PCs1, by = "Taxa", all.x = TRUE)
     
  	  if (length(meta_sig_count) > 0) {
      myGD_sig_df <- data.frame(Taxa = myGD$Taxa, myGD[, sig_snps, drop=FALSE])
      train_df3 <- merge(train_df3, myGD_sig_df, by = "Taxa", all.x = TRUE)
      test_df3 <- merge(test_df3, myGD_sig_df, by = "Taxa", all.x = TRUE)
}
      lm_model <- lm(train_all$Phenotype ~ ., data = train_df3[, -1])
      pred_values <- predict(lm_model, newdata = test_df3[, -1])
      cor_pred_meta <- cor(test_all$Phenotype, pred_values, use = "complete.obs")
      cor_gBV_meta  <- cor_gBV_blink 
      cat("Meta预测与真实表型的相关性为:", cor_pred_meta, "\n")

      results_all <- rbind(
        results_all,
        data.frame(
          trait = trait,
          rep = i,
          fold = j,
          method = "Meta",
          cor_pred = cor_pred_meta,
          cor_gBV = cor_gBV_meta,
          QTNs = meta_sig_count,
          stringsAsFactors = FALSE
        )
      )
      
     cat(sprintf(
  "Trait=%s | rep=%d | fold=%d | BLINK cor(Pred)=%.4f | BLINK cor(gBV)=%.4f | Meta cor(Pred)=%.4f | Meta cor(gBV)=%.4f | SNPs=%d\n",
  trait, i, j,
  ifelse(is.na(cor_pred_blink), NaN, cor_pred_blink),
  ifelse(is.na(cor_gBV_blink), NaN, cor_gBV_blink),
  ifelse(is.na(cor_pred_meta), NaN, cor_pred_meta),
  ifelse(is.na(cor_gBV_meta), NaN, cor_gBV_meta),
  meta_sig_count
))

    }
  }
}

if (!dir.exists("D:/gs results")) dir.create("D:/gs results", recursive = TRUE)
write.csv(results_all, "D:/gs results/6.0CW_blink_vs_meta_randomGroups_Sigs.csv", row.names = FALSE)
