# 04_enrichment_GO_KEGG_GSEA.R
# GO BP over-representation (up/down), KEGG, and GSEA

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

de_f <- readRDS("data/de_f.rds")

# Symbol -> Entrez
syms <- unique(na.omit(de_f$Symbol))
sym2ent <- mapIds(org.Hs.eg.db, keys = syms, column = "ENTREZID",
                  keytype = "SYMBOL", multiVals = "first")
de_f$Entrez <- sym2ent[de_f$Symbol]

universe <- unique(na.omit(de_f$Entrez))   # genes actually tested (15,035)
up   <- unique(na.omit(de_f$Entrez[de_f$sig == "Up"]))     # 1,315
down <- unique(na.omit(de_f$Entrez[de_f$sig == "Down"]))   # 205
all_sig <- union(up, down)

# GO Biological Process
run_go <- function(g) enrichGO(gene = g, universe = universe, OrgDb = org.Hs.eg.db,
                               ont = "BP", pAdjustMethod = "BH",
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                               readable = TRUE)
ego_all  <- run_go(all_sig)
ego_up   <- run_go(up)
ego_down <- run_go(down)

write.csv(as.data.frame(ego_all),  "results/GO_BP_all.csv",  row.names = FALSE)
write.csv(as.data.frame(ego_up),   "results/GO_BP_up.csv",   row.names = FALSE)
write.csv(as.data.frame(ego_down), "results/GO_BP_down.csv", row.names = FALSE)

ggsave("figures/GO_dotplot_up.png",
       dotplot(ego_up, showCategory = 15) + ggtitle("GO BP: upregulated in tumour"),
       width = 8, height = 7, dpi = 300)
ggsave("figures/GO_dotplot_down.png",
       dotplot(ego_down, showCategory = 15) + ggtitle("GO BP: downregulated in tumour"),
       width = 8, height = 7, dpi = 300)

# Simplified (remove redundant terms)
ego_up_s   <- simplify(ego_up,   cutoff = 0.7, by = "p.adjust", select_fun = min)
ego_down_s <- simplify(ego_down, cutoff = 0.7, by = "p.adjust", select_fun = min)
write.csv(as.data.frame(ego_up_s),   "results/GO_BP_up_simplified.csv",   row.names = FALSE)
write.csv(as.data.frame(ego_down_s), "results/GO_BP_down_simplified.csv", row.names = FALSE)
ggsave("figures/GO_dotplot_up_simplified.png",
       dotplot(ego_up_s, showCategory = 15) + ggtitle("GO BP up (simplified)"),
       width = 8, height = 7, dpi = 300)

# KEGG (needs internet)
kegg_up   <- enrichKEGG(gene = up,   universe = universe, organism = "hsa", pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene = down, universe = universe, organism = "hsa", pvalueCutoff = 0.05)
kegg_up   <- setReadable(kegg_up,   org.Hs.eg.db, keyType = "ENTREZID")
kegg_down <- setReadable(kegg_down, org.Hs.eg.db, keyType = "ENTREZID")
write.csv(as.data.frame(kegg_up),   "results/KEGG_up.csv",   row.names = FALSE)
write.csv(as.data.frame(kegg_down), "results/KEGG_down.csv", row.names = FALSE)

# GSEA on full ranked list (moderated t, one probe per gene with strongest |t|)
rank_df <- de_f[!is.na(de_f$Entrez), ]
rank_df <- rank_df[order(-abs(rank_df$t)), ]
rank_df <- rank_df[!duplicated(rank_df$Entrez), ]
geneList <- sort(setNames(rank_df$t, rank_df$Entrez), decreasing = TRUE)

set.seed(123)
gsea_go <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, ont = "BP",
                 minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05,
                 eps = 0, verbose = FALSE)
gsea_go <- setReadable(gsea_go, org.Hs.eg.db, keyType = "ENTREZID")
write.csv(as.data.frame(gsea_go), "results/GSEA_GO_BP.csv", row.names = FALSE)

ggsave("figures/GSEA_GO_dotplot.png",
       dotplot(gsea_go, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign),
       width = 10, height = 8, dpi = 300)

g <- as.data.frame(gsea_go)
table(sign(g$NES))   # -1: 80, +1: 1305