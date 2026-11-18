# 06_acinar_sensitivity.R
# Sensitivity analysis: adjust for acinar content (bulk-tissue composition)

library(limma)
library(ggplot2)

de_f <- readRDS("data/de_f.rds")
obj <- readRDS("data/qc_objects.rds")
ex_f <- obj$ex_f
info2 <- obj$info2

ac_probes <- c("205615_at","205912_at","205869_at","205971_s_at","210080_x_at")
info2$Acinar <- colMeans(t(scale(t(ex_f[ac_probes, ]))))

# PCA coloured by acinar score
pca2 <- prcomp(t(obj$ex2), scale. = TRUE)
df3 <- data.frame(pca2$x[, 1:2], Tissue = info2$Tissue, Acinar = info2$Acinar)
p3 <- ggplot(df3, aes(PC1, PC2, colour = Acinar, shape = Tissue)) +
  geom_point(size = 3) + scale_colour_viridis_c() + theme_bw()
ggsave("figures/PCA_acinar_score.png", p3, width = 7, height = 5, dpi = 300)

# limma with acinar score as covariate
design3 <- model.matrix(~ Patient + Acinar + Tissue, data = info2)
fit3 <- eBayes(lmFit(ex_f, design3))
de3 <- topTable(fit3, coef = "Tissuetumor", number = Inf, adjust.method = "BH")
write.csv(de3, "results/DE_results_acinar_adjusted.csv")

sig_old <- rownames(de_f)[de_f$sig != "NS"]
sig_new <- rownames(de3)[de3$adj.P.Val < 0.05 & abs(de3$logFC) >= 1]
length(sig_old)                        # 2,253
length(sig_new)                        # 1,859
length(intersect(sig_old, sig_new))    # 1,841

# Which probes were lost after adjustment?
lost <- setdiff(sig_old, sig_new)
table(de_f[lost, "sig"])
head(de_f[lost, c("Symbol","logFC","adj.P.Val")][order(de_f[lost, "logFC"]), ], 15)

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")