# 01_download_and_organize.R
# Load the patient-level expression matrix (replicates already averaged)
# and rebuild the patient table from the column names.

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# --- Record of how the averaged matrix was originally made -----------------
# library(GEOquery)
# geo_data <- GEOquery::getGEO("GSE15471", GSEMatrix = TRUE)
# eset <- geo_data[[1]]            # 54,675 probes x 78 hybridizations
# Replicate hybridizations (3 patients) were checked by Pearson correlation
# (0.987-0.994; results/replicate_correlation.csv) and averaged within
# patient x tissue -> 72 patient-level profiles (36 normal + 36 tumour).
# The result was saved as data/pancreatic_expression_averaged.rds
# ---------------------------------------------------------------------------

expr_averaged <- readRDS("data/pancreatic_expression_averaged.rds")
dim(expr_averaged)   # 54675 x 72

# Column names look like "51681_normal" / "51681_tumor"
patient_info <- data.frame(
  Patient = sub("_.*$", "", colnames(expr_averaged)),
  Tissue  = sub("^.*_", "", colnames(expr_averaged)),
  stringsAsFactors = FALSE
)

stopifnot(all(patient_info$Tissue %in% c("normal", "tumor")))
table(patient_info$Tissue)   # 36 normal, 36 tumor

saveRDS(list(expr_averaged = expr_averaged, patient_info = patient_info),
        "data/organized_data.rds")