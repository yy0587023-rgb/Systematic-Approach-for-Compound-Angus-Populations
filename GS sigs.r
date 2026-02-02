setwd("D:/gs/R significant sites") 
myY_hb  <- read.table("hb.txt",  head = TRUE, stringsAsFactors = FALSE)[,1:4]
myY_qy  <- read.table("qy.txt",  head = TRUE, stringsAsFactors = FALSE)[,1:4]
myY_hlj <- read.table("hlj.txt", head = TRUE, stringsAsFactors = FALSE)[,1:4]
myY_hlj=myY_hlj[,1:4]
myY_hb=myY_hb[,1:4]
myY_hlj[1:10,]
myY_hb[1:10,]
myGD <- read.csv("myGD.csv", head = TRUE)
myGD_data <- myGD[, -1]
colnames(myGD)[1] <- "Taxa"
myGM <- read.csv("myGM.csv", head = TRUE)
myCV <- read.table("CV.txt", head = TRUE)[,1:3]
remove_id<- c("D00592531", "D00608064")
myY_hb <- myY_hb[!myY_hb$Taxa %in% remove_id, ]
dim(myY_hb)
remove_id1<- c("D00592537", "D00628064")
myGD<- myGD[!myGD$Taxa %in% remove_id1, ]
dim(myGD)
myCV<- myCV[!myCV$Taxa %in% remove_id, ]
dim(myCV)
source("D:/gs/R significant sites/gapit_functions (6).txt") ##newest version
metal_path <- "D:/gs/generic-metal/metal.exe"
nrep <- 1
nfold <- 5
cutOff <- 0.01
nSNP <- nrow(myGM)   
sig_thre <- cutOff / nSNP

gd_mat <- as.matrix(myGD[,-1])
gd_var <- apply(gd_mat, 2, var, na.rm=TRUE) 
gd_mat <- gd_mat[, gd_var > 0, drop=FALSE]
pca_res <- prcomp(gd_mat, scale. = TRUE)
PCs <- pca_res$x[,1:3]
colnames(PCs) <- paste0("PC", 1:3)
PCs_df <- data.frame(Taxa = myGD$Taxa, PCs)
head(PCs)

myCV1 <- data.frame(Taxa = myGD$Taxa, PCs)
head(myCV1)

PCs1<- data.frame(Taxa = myGD$Taxa, PCs)
myCV2<- merge(myCV, PCs1, by = "Taxa", all.x = TRUE)
head(myCV2)

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
##index+pca+sigs
traits <- c("CW")###
set.seed(1000)
for (trait in traits) {
  cat("========== 正在预测性状:", trait, "==========\n")
  
  for (i in 1:nrep) {
  
    sets_hb  <- sample(cut(1:nrow(myY_hb),  nfold, labels = FALSE), nrow(myY_hb))
    sets_qy  <- sample(cut(1:nrow(myY_qy),  nfold, labels = FALSE), nrow(myY_qy))
    sets_hlj <- sample(cut(1:nrow(myY_hlj), nfold, labels = FALSE), nrow(myY_hlj))
    
    for (j in 1:nfold) {
      cat("rep=", i, " fold=", j, " trait=", trait, "\n")
     
      train_hb  <- myY_hb[ sets_hb  != j, c("Taxa", trait)]
      test_hb   <- myY_hb[ sets_hb  == j, c("Taxa", trait)]
      train_qy  <- myY_qy[ sets_qy  != j, c("Taxa", trait)]
      test_qy   <- myY_qy[ sets_qy  == j, c("Taxa", trait)]
      train_hlj <- myY_hlj[sets_hlj != j, c("Taxa", trait)]
      test_hlj  <- myY_hlj[sets_hlj == j, c("Taxa", trait)]

      names(train_hb)[2]  <- "Phenotype"
      names(train_qy)[2]  <- "Phenotype"
      names(train_hlj)[2] <- "Phenotype"
      names(test_hb)[2]   <- "Phenotype"
      names(test_qy)[2]   <- "Phenotype"
      names(test_hlj)[2]  <- "Phenotype"
      
      # =================== (1) blink ===================
      train_all <- rbind(train_hb, train_qy, train_hlj)
      colnames(train_all) <- c("Taxa", "Phenotype")
      
      myGAPIT1 <- GAPIT(
        Y = train_all,
        GD = myGD,
        GM = myGM,
        CV = myCV2,
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

      test_all <- rbind(test_hb, test_qy, test_hlj)
	  
      train_df1 <- merge(train_all, myCV2[, setdiff(names(myCV2), "Phenotype")], by = "Taxa")
      train_df1$Phenotype <- NULL  
      test_df1 <- merge(test_all, myCV2[, setdiff(names(myCV2), "Phenotype")], by = "Taxa")
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
      train_df2 <- merge(train_df2, PCs_df, by = "Taxa", all.x = TRUE)
     
      test_df2 <- data.frame(Taxa = test_all$Taxa)
      test_df2 <- merge(test_df2, PCs_df, by = "Taxa", all.x = TRUE)
     
      train_pheno <- train_all$Phenotype
      lm_model <- lm(train_pheno ~ ., data = train_df2[, -1])
      pred_values <- predict(lm_model, newdata = test_df2[, -1])

      cor_gBV_blink  <- cor(test_all$Phenotype, pred_values, use = "complete.obs")
	  
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

      ### =================== (2)meta ===================
      run_gwas <- function(train_df) {
        res <- GAPIT(
          Y = train_df,
          GD = myGD,
          GM = myGM,
          CV = myCV1,
          cutOff = cutOff,
          model = "Blink",
		  Multi_iter = FALSE,
		  Random.model = FALSE
        )
        return(res$GWAS)
      }

      gwas_hb  <- run_gwas(train_hb)
      gwas_qy  <- run_gwas(train_qy)
      gwas_hlj <- run_gwas(train_hlj)

      fn_hb  <- paste0("hb_rep", i, "_fold", j, ".txt")
      fn_qy  <- paste0("qy_rep", i, "_fold", j, ".txt")
      fn_hlj <- paste0("hlj_rep", i, "_fold", j, ".txt")

      write.table(gwas_hb[, c("SNP", "effect", "P.value", "nobs")], fn_hb,  sep = "\t", row.names = FALSE, quote = FALSE)
      write.table(gwas_qy[, c("SNP", "effect", "P.value", "nobs")], fn_qy,  sep = "\t", row.names = FALSE, quote = FALSE)
      write.table(gwas_hlj[, c("SNP", "effect", "P.value", "nobs")], fn_hlj, sep = "\t", row.names = FALSE, quote = FALSE)

      metal_script <- paste0("metal_rep", i, "_fold", j, ".txt")
      meta_out <- paste0("meta_rep", i, "_fold", j)

      ms <- c(
        "SCHEME SAMPLESIZE",
        "MARKER SNP",
        "EFFECT effect",
        "PVALUE P.value",
        "WEIGHT nobs",
        paste0("PROCESS ", fn_hb),
        paste0("PROCESS ", fn_qy),
        paste0("PROCESS ", fn_hlj),
        paste0("OUTFILE ", meta_out),
        "ANALYZE",
        "QUIT"
      )
      writeLines(ms, metal_script)

      system(paste(metal_path, metal_script), ignore.stdout = FALSE, ignore.stderr = FALSE)

      meta_tbl <- "METAANALYSIS1.TBL"
      if (!file.exists(meta_tbl)) stop("METAANALYSIS1.TBL 不存在！")

      meta_res <- read.table(meta_tbl, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      sig_res  <- subset(meta_res, P.value < sig_thre)
      meta_sig_count <- nrow(sig_res)
      sig_snps <- intersect(sig_res$MarkerName, colnames(myGD_data))

      train_df3 <- data.frame(Taxa = train_all$Taxa)
      train_df3 <- merge(train_df3, PCs_df, by = "Taxa", all.x = TRUE)
     
      test_df3 <- data.frame(Taxa = test_all$Taxa)
      test_df3 <- merge(test_df3, PCs_df, by = "Taxa", all.x = TRUE)
     
  	  if (length(meta_sig_count) > 0) {
      myGD_sig_df <- data.frame(Taxa = myGD$Taxa, myGD[, sig_snps, drop=FALSE])
      train_df3 <- merge(train_df3, myGD_sig_df, by = "Taxa", all.x = TRUE)
      test_df3 <- merge(test_df3, myGD_sig_df, by = "Taxa", all.x = TRUE)
}
      lm_model <- lm(train_all$Phenotype ~ ., data = train_df3[, -1])
      pred_values <- predict(lm_model, newdata = test_df3[, -1])
      cor_pred_meta <- cor(test_all$Phenotype, pred_values, use = "complete.obs")
      cor_gBV_meta  <- cor_gBV_blink 

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

          cat(sprintf("Trait=%s | rep=%d | fold=%d | BLINK cor(Pred)=%.4f cor(gBV)=%.4f | Meta cor(Pred)=%.4f cor(gBV)=%.4f\n",      
                  trait, i, j,
                  ifelse(is.na(cor_pred_blink), NaN, cor_pred_blink),
				  ifelse(is.na(cor_gBV_blink), NaN, cor_gBV_blink),
                  ifelse(is.na(cor_pred_meta), NaN, cor_pred_meta),
				  ifelse(is.na(cor_gBV_meta), NaN, cor_gBV_meta)
				  ))
    }
  }
}

if (!dir.exists("D:/gs results")) dir.create("D:/gs results", recursive = TRUE)
write.csv(results_all, "D:/gs results/6.0CW_blink_vs_meta(significants).csv", row.names = FALSE)

