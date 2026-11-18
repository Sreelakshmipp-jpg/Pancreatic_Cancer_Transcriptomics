# 05_candidate_genes_and_figures.R
# Candidate gene table, heatmap, paired plots (candidates and acinar genes)

library(ggplot2)
library(pheatmap)

de_f <- readRDS("data/de_f.rds")
obj <- readRDS("data/qc_objects.rds")
ex_f <- obj$ex_f
info2 <- obj$info2

cands <- c("AHNAK2","CDH3","CEACAM1","CEACAM5","CEACAM6","COL10A1","COL11A1",
           "FN1","KRT19","KRT16","MMP11","MMP12","VCAN","TRIM29","KRAS")

# Candidate gene results (all probes)
cand_res <- de_f[de_f$Symbol %in% cands, c("ProbeID","Symbol","logFC","adj.P.Val","sig")]
cand_res <- cand_res[order(cand_res$Symbol), ]
write.csv(cand_res, "results/candidate_genes_35pairs.csv", row.names = FALSE)

# Top 25 up + top 25 down genes (by t-statistic, one probe per gene)
top <- de_f[!is.na(de_f$Symbol) & de_f$sig != "NS", ]
top <- top[!duplicated(top$Symbol), ]
sel <- rbind(head(top[order(-top$t), ], 25), head(top[order(top$t), ], 25))

mat <- ex_f[sel$ProbeID, ]
rownames(mat) <- sel$Symbol
mat <- t(scale(t(mat)))
ann <- data.frame(Tissue = info2$Tissue, row.names = colnames(mat))

png("figures/heatmap_top50.png", width = 2000, height = 1800, res = 250)
pheatmap(mat, annotation_col = ann, show_colnames = FALSE,
         clustering_method = "ward.D2", fontsize_row = 7)
dev.off()

# Paired candidate plot (probe with highest mean expression per gene)
pick <- de_f[de_f$Symbol %in% cands, ]
pick$MeanExpr <- rowMeans(ex_f[pick$ProbeID, ])
pick <- pick[order(-pick$MeanExpr), ]
pick <- pick[!duplicated(pick$Symbol), ]
write.csv(pick[, c("Symbol","ProbeID","logFC","adj.P.Val","MeanExpr")],
          "results/candidate_probe_used.csv", row.names = FALSE)

long <- do.call(rbind, lapply(seq_len(nrow(pick)), function(i)
  data.frame(Symbol = pick$Symbol[i], Expr = ex_f[pick$ProbeID[i], ],
             Patient = info2$Patient, Tissue = info2$Tissue)))

pp <- ggplot(long, aes(Tissue, Expr, group = Patient)) +
  geom_line(colour = "grey70", alpha = 0.6) +
  geom_point(aes(colour = Tissue), size = 1.5) +
  scale_colour_manual(values = c(normal = "steelblue", tumor = "tomato")) +
  facet_wrap(~Symbol, scales = "free_y") +
  theme_bw() + theme(legend.position = "none")
ggsave("figures/candidate_genes_paired.png", pp, width = 11, height = 8, dpi = 300)

# Acinar genes (not significant; bimodal expression)
acinar <- de_f[de_f$Symbol %in% c("CPA1","PRSS1","CELA3A","PNLIP","CTRB1"),
               c("Symbol","logFC","adj.P.Val")]
print(acinar)

ac <- c("205615_at","205912_at","205869_at","205971_s_at","210080_x_at")
long_ac <- do.call(rbind, lapply(ac, function(p)
  data.frame(Probe = p, Symbol = de_f[p, "Symbol"], Expr = ex_f[p, ],
             Patient = info2$Patient, Tissue = info2$Tissue)))

pa <- ggplot(long_ac, aes(Tissue, Expr, group = Patient)) +
  geom_line(colour = "grey70", alpha = 0.6) +
  geom_point(aes(colour = Tissue), size = 1.5) +
  scale_colour_manual(values = c(normal = "steelblue", tumor = "tomato")) +
  facet_wrap(~Symbol, scales = "free_y") +
  theme_bw() + theme(legend.position = "none")
ggsave("figures/acinar_genes_paired.png", pa, width = 9, height = 6, dpi = 300)