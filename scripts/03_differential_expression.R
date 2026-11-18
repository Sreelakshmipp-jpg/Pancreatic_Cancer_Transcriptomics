# 03_differential_expression.R
# Paired limma analysis (35 pairs), annotation, volcano plot

library(limma)
library(ggplot2)
library(hgu133plus2.db)
library(AnnotationDbi)

obj <- readRDS("data/qc_objects.rds")
ex_f <- obj$ex_f
info2 <- obj$info2
design2 <- obj$design2

fit_f <- eBayes(lmFit(ex_f, design2))
de_f <- topTable(fit_f, coef = "Tissuetumor", number = Inf, adjust.method = "BH")

sum(de_f$adj.P.Val < 0.05)                              # 22,311
sum(de_f$adj.P.Val < 0.05 & abs(de_f$logFC) >= 1)       # 2,253

# Annotation
de_f$ProbeID <- rownames(de_f)
de_f$Symbol <- mapIds(hgu133plus2.db, keys = de_f$ProbeID,
                      column = "SYMBOL", keytype = "PROBEID",
                      multiVals = "first")

de_f$sig <- ifelse(de_f$adj.P.Val < 0.05 & de_f$logFC >= 1, "Up",
                   ifelse(de_f$adj.P.Val < 0.05 & de_f$logFC <= -1, "Down", "NS"))
table(de_f$sig)   # Up 1942, Down 311

write.csv(de_f, "results/DE_results_filtered_35pairs.csv", row.names = FALSE)

# Volcano plot
v <- ggplot(de_f, aes(logFC, -log10(adj.P.Val), colour = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = c(Up = "tomato", Down = "steelblue", NS = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_bw()
ggsave("figures/volcano_35pairs.png", v, width = 7, height = 6, dpi = 300)

saveRDS(de_f, "data/de_f.rds")