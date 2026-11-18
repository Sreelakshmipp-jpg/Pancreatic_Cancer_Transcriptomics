# 02_QC_outlier_filtering.R
# QC, outlier removal (patient 51681), low-expression filtering

library(limma)
library(ggplot2)

obj <- readRDS("data/organized_data.rds")
ex <- obj$expr_averaged
patient_info <- obj$patient_info

# 1. Scale check (2-16 = log2 scale, RMA-normalized)
range(ex)

# 2. Sample correlation: identify outliers
cor_mat <- cor(ex)
sort(rowMeans(cor_mat))[1:5]   # 51681_normal = 0.709; all others >= 0.899

# 3. Remove outlier patient (both samples, keeps paired design -> 35 pairs)
bad_patient <- "51681"
keep_samples <- patient_info$Patient != bad_patient
ex2 <- ex[, keep_samples]
info2 <- patient_info[keep_samples, ]
info2$Patient <- droplevels(factor(info2$Patient))
info2$Tissue <- factor(info2$Tissue, levels = c("normal", "tumor"))
dim(ex2)   # 54675 x 70

# 4. QC plots (35 pairs)
png("figures/QC_boxplot_35pairs.png", width = 1400, height = 700)
boxplot(ex2, las = 2, outline = FALSE, cex.axis = 0.5,
        col = ifelse(info2$Tissue == "tumor", "tomato", "steelblue"),
        main = "Expression per sample (35 pairs)")
dev.off()

png("figures/QC_density_35pairs.png", width = 900, height = 600)
plotDensities(ex2, group = info2$Tissue, legend = "topright")
dev.off()

pca2 <- prcomp(t(ex2), scale. = TRUE)
pct2 <- round(100 * summary(pca2)$importance[2, 1:2], 1)
df2 <- data.frame(pca2$x[, 1:2], Tissue = info2$Tissue)
p2 <- ggplot(df2, aes(PC1, PC2, colour = Tissue)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1 (", pct2[1], "%)"), y = paste0("PC2 (", pct2[2], "%)")) +
  theme_bw()
ggsave("figures/QC_PCA_no_outlier.png", p2, width = 7, height = 5, dpi = 300)

# 5. Low-expression filter: above global median in >= 18 of 70 samples
keep <- rowSums(ex2 > median(ex2)) >= 18
ex_f <- ex2[keep, ]
dim(ex_f)   # 30,529 probes

design2 <- model.matrix(~ Patient + Tissue, data = info2)

saveRDS(list(ex2 = ex2, ex_f = ex_f, info2 = info2, design2 = design2),
        "data/qc_objects.rds")