#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else
  "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(project_dir, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(scuttle)
  library(scDblFinder)
  library(BiocParallel)
  library(readxl)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

raw_dir <- file.path(project_dir, "data", "raw", "scrna", "GSE183904")
matrix_dir <- file.path(raw_dir, "extracted")
supp_xlsx <- file.path(raw_dir, "paper_supplement", "cd-21-0683_supp.xlsx")
soft_file <- file.path(raw_dir, "GSE183904_family.soft.gz")
processed_dir <- file.path(project_dir, "data", "processed", "scrna", "GSE183904")
intermediate_dir <- file.path(project_dir, "data", "intermediate", "scrna", "GSE183904", "sample_sce")
result_dir <- file.path(project_dir, "results", "scrna", "GSE183904")
figure_dir <- file.path(project_dir, "figures", "supplementary")
audit_dir <- file.path(project_dir, "audit")
source_dir <- file.path(project_dir, "source_data")
for (d in c(processed_dir, intermediate_dir, result_dir, figure_dir, audit_dir, source_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Parse the exact patient-to-file mapping embedded in the official GEO SOFT.
soft <- readLines(gzfile(soft_file), warn = FALSE)
mapping_lines <- sub("^!Series_summary = ", "", soft[grepl("^!Series_summary = NGCII", soft)])
mapping_parts <- strsplit(mapping_lines, "\t", fixed = TRUE)
mapping_cols <- c("Patient_ID", "Primary_Tumor", "Primary_Normal", "Peritoneal_Tumor", "Peritoneal_Normal")
mapping <- rbindlist(lapply(mapping_parts, function(x) {
  length(x) <- length(mapping_cols)
  as.list(x)
}))
setnames(mapping, mapping_cols)
mapping[, Patient_ID := trimws(Patient_ID)]

mapping_long <- melt(mapping, id.vars = "Patient_ID", variable.name = "specimen_class",
                     value.name = "submitted_name")
mapping_long <- mapping_long[!is.na(submitted_name) & submitted_name != "-"]
mapping_long[, sample_id := sub("\\.csv$", "", submitted_name)]
mapping_long[, sample_number := as.integer(sub("^sample", "", sample_id))]
mapping_long[, tissue_site := fifelse(grepl("^Primary", specimen_class), "stomach", "peritoneum")]
mapping_long[, disease_status := fifelse(grepl("Tumor$", specimen_class), "Tumor", "Normal")]

# Read-only extraction of author Supplementary Tables 1 and 2. The original
# workbook is preserved unchanged; flat TSVs make joins auditable.
clinical <- as.data.table(read_excel(supp_xlsx, sheet = "Supp Table 1"))
setnames(clinical, 1L, "author_row")
setnames(clinical, trimws(names(clinical)))
clinical[, `Patient ID` := trimws(`Patient ID`)]
markers <- as.data.table(read_excel(supp_xlsx, sheet = "Supp Table 2", skip = 1L))
setnames(markers, trimws(names(markers)))
fwrite(clinical, file.path(processed_dir, "GSE183904_author_clinical_supp_table1.tsv"), sep = "\t", na = "")
fwrite(markers, file.path(processed_dir, "GSE183904_author_marker_supp_table2.tsv"), sep = "\t", na = "")

files <- data.table(file_path = sort(list.files(matrix_dir, pattern = "^GSM[0-9]+_sample[0-9]+\\.csv\\.gz$",
                                                full.names = TRUE)))
files[, file_name := basename(file_path)]
files[, GSM := sub("_.*$", "", file_name)]
files[, sample_id := sub("^GSM[0-9]+_(sample[0-9]+)\\.csv\\.gz$", "\\1", file_name)]
files[, sample_number := as.integer(sub("^sample", "", sample_id))]
if (nrow(files) != 40L || uniqueN(files$sample_id) != 40L) stop("Expected 40 unique GSE183904 matrices")

manifest <- merge(files, mapping_long[, .(Patient_ID, specimen_class, submitted_name, sample_id,
                                           tissue_site, disease_status)], by = "sample_id", all.x = TRUE)
manifest <- merge(manifest, clinical, by.x = "Patient_ID", by.y = "Patient ID", all.x = TRUE)
if (manifest[, anyNA(Patient_ID)] || uniqueN(manifest$sample_id) != 40L) stop("Incomplete GEO patient mapping")

count_header_cells <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  length(strsplit(readLines(con, n = 1L), ",", fixed = TRUE)[[1L]]) - 1L
}
manifest[, author_retained_cells := vapply(file_path, count_header_cells, integer(1))]
manifest[, primary_gastric_tumor := specimen_class == "Primary_Tumor"]
manifest[, include_primary_replication := primary_gastric_tumor]
setorder(manifest, sample_number)
fwrite(manifest[, !"file_path"], file.path(processed_dir, "GSE183904_sample_patient_clinical_manifest.tsv"),
       sep = "\t", na = "")

robust_thresholds <- function(sum_counts, detected, mito_percent) {
  log_sum <- log10(sum_counts + 1)
  log_det <- log10(detected + 1)
  med_sum <- median(log_sum); mad_sum <- mad(log_sum, constant = 1.4826)
  med_det <- median(log_det); mad_det <- mad(log_det, constant = 1.4826)
  med_mt <- median(mito_percent); mad_mt <- mad(mito_percent, constant = 1.4826)
  if (!is.finite(mad_sum) || mad_sum == 0) mad_sum <- IQR(log_sum) / 1.349
  if (!is.finite(mad_det) || mad_det == 0) mad_det <- IQR(log_det) / 1.349
  if (!is.finite(mad_mt) || mad_mt == 0) mad_mt <- IQR(mito_percent) / 1.349
  if (!is.finite(mad_sum)) mad_sum <- 0
  if (!is.finite(mad_det)) mad_det <- 0
  if (!is.finite(mad_mt)) mad_mt <- 0
  list(
    min_sum = max(0, 10^(med_sum - 4 * mad_sum) - 1),
    min_detected = max(0, 10^(med_det - 4 * mad_det) - 1),
    review_max_sum = min(10^(med_sum + 5 * mad_sum) - 1,
                         as.numeric(quantile(sum_counts, 0.999))),
    review_max_detected = min(10^(med_det + 5 * mad_det) - 1,
                              as.numeric(quantile(detected, 0.999))),
    max_mito_percent = min(100, med_mt + 5 * mad_mt)
  )
}

qc_rows <- list()
qc_summaries <- list()
threshold_rows <- list()

for (i in seq_len(nrow(manifest))) {
  info <- manifest[i]
  checkpoint <- file.path(intermediate_dir, paste0(info$sample_id, "_secondary_qc_sce.rds"))
  message(sprintf("[%d/40] %s (%s; %d author-retained cells)", i, info$sample_id,
                  info$Patient_ID, info$author_retained_cells))
  if (file.exists(checkpoint)) {
    sce <- readRDS(checkpoint)
    cd0 <- as.data.table(as.data.frame(colData(sce)))
    thr <- robust_thresholds(cd0$sum_counts, cd0$detected_genes, cd0$mito_percent)
    passes_low_qc <- cd0$sum_counts >= thr$min_sum & cd0$detected_genes >= thr$min_detected &
      cd0$mito_percent <= thr$max_mito_percent
    upper_review <- cd0$sum_counts > thr$review_max_sum | cd0$detected_genes > thr$review_max_detected
    colData(sce)$passes_secondary_low_quality_qc <- passes_low_qc
    colData(sce)$upper_count_review_flag <- upper_review
    colData(sce)$passes_secondary_qc_singlet <- passes_low_qc &
      colData(sce)$secondary_scDblFinder_class == "singlet"
    metadata(sce)$secondary_qc_thresholds <- thr
    saveRDS(sce, checkpoint, compress = "gzip")
    rm(cd0)
  } else {
    dt <- fread(info$file_path, check.names = FALSE, showProgress = FALSE)
    genes <- as.character(dt[[1L]])
    if (anyDuplicated(genes)) stop("Duplicated genes in ", info$file_name)
    dense <- as.matrix(dt[, -1L])
    storage.mode(dense) <- "integer"
    counts <- Matrix(dense, sparse = TRUE)
    rownames(counts) <- genes
    colnames(counts) <- paste0(colnames(dense), "__", info$sample_id)
    rm(dt, dense); gc(verbose = FALSE)

    sce <- SingleCellExperiment(list(counts = counts))
    is_mito <- grepl("^MT-", rownames(sce))
    qc <- perCellQCMetrics(sce, subsets = list(Mito = is_mito))
    thr <- robust_thresholds(qc$sum, qc$detected, qc$subsets_Mito_percent)
    passes_low_qc <- qc$sum >= thr$min_sum & qc$detected >= thr$min_detected &
      qc$subsets_Mito_percent <= thr$max_mito_percent
    upper_review <- qc$sum > thr$review_max_sum | qc$detected > thr$review_max_detected

    dbl_score <- rep(NA_real_, ncol(sce))
    dbl_class <- rep("not_tested_low_quality", ncol(sce))
    dbl_error <- NA_character_
    if (sum(passes_low_qc) >= 100L) {
      set.seed(183904L + info$sample_number)
      tested <- tryCatch(
        scDblFinder(sce[, passes_low_qc, drop = FALSE], clusters = FALSE, verbose = FALSE,
                    BPPARAM = SerialParam()),
        error = function(e) e
      )
      if (inherits(tested, "error")) {
        dbl_error <- conditionMessage(tested)
        dbl_class[passes_low_qc] <- "scDblFinder_error"
      } else {
        dbl_score[passes_low_qc] <- colData(tested)$scDblFinder.score
        dbl_class[passes_low_qc] <- as.character(colData(tested)$scDblFinder.class)
      }
    } else {
      dbl_error <- "Fewer than 100 cells passed secondary low-quality thresholds"
    }

    colData(sce)$sample_id <- info$sample_id
    colData(sce)$GSM <- info$GSM
    colData(sce)$Patient_ID <- info$Patient_ID
    colData(sce)$specimen_class <- info$specimen_class
    colData(sce)$tissue_site <- info$tissue_site
    colData(sce)$disease_status <- info$disease_status
    colData(sce)$sum_counts <- qc$sum
    colData(sce)$detected_genes <- qc$detected
    colData(sce)$mito_percent <- qc$subsets_Mito_percent
    colData(sce)$passes_secondary_low_quality_qc <- passes_low_qc
    colData(sce)$upper_count_review_flag <- upper_review
    colData(sce)$secondary_scDblFinder_score <- dbl_score
    colData(sce)$secondary_scDblFinder_class <- dbl_class
    colData(sce)$passes_secondary_qc_singlet <- passes_low_qc & dbl_class == "singlet"
    metadata(sce)$secondary_qc_thresholds <- thr
    metadata(sce)$scDblFinder_error <- dbl_error
    metadata(sce)$input_note <- paste(
      "Official GEO matrix contains author-retained cells after the publication QC/DoubletFinder pipeline;",
      "this object adds data-derived low-quality flags and a secondary scDblFinder diagnostic."
    )
    saveRDS(sce, checkpoint, compress = "gzip")
  }

  cd <- as.data.table(as.data.frame(colData(sce)), keep.rownames = "cell_id")
  qc_rows[[i]] <- cd
  thr <- metadata(sce)$secondary_qc_thresholds
  threshold_rows[[i]] <- data.table(
    sample_id = info$sample_id,
    min_sum = thr$min_sum,
    min_detected = thr$min_detected,
    review_max_sum = thr$review_max_sum,
    review_max_detected = thr$review_max_detected,
    max_mito_percent = thr$max_mito_percent
  )
  dbl_err <- metadata(sce)$scDblFinder_error
  if (is.null(dbl_err) || !length(dbl_err)) dbl_err <- NA_character_
  qc_summaries[[i]] <- data.table(
    sample_id = info$sample_id, GSM = info$GSM, Patient_ID = info$Patient_ID,
    specimen_class = info$specimen_class,
    input_author_retained_cells = ncol(sce),
    secondary_low_quality_fail = sum(!cd$passes_secondary_low_quality_qc),
    secondary_upper_review_flag = sum(cd$upper_count_review_flag),
    secondary_scDblFinder_doublet = sum(cd$secondary_scDblFinder_class == "doublet"),
    secondary_scDblFinder_singlet = sum(cd$secondary_scDblFinder_class == "singlet"),
    final_secondary_qc_singlet = sum(cd$passes_secondary_qc_singlet),
    scDblFinder_error = dbl_err
  )
  rm(sce, cd); gc(verbose = FALSE)
}

qc_all <- rbindlist(qc_rows, fill = TRUE)
thresholds <- rbindlist(threshold_rows)
qc_summary <- rbindlist(qc_summaries, fill = TRUE)
fwrite(qc_all, file.path(processed_dir, "GSE183904_secondary_cell_qc_metadata.tsv.gz"), sep = "\t", na = "")
fwrite(thresholds, file.path(result_dir, "GSE183904_secondary_qc_thresholds.tsv"), sep = "\t")
fwrite(qc_summary, file.path(result_dir, "GSE183904_secondary_qc_summary.tsv"), sep = "\t", na = "")

plot_stats <- qc_all[, .(
  q25_sum = as.numeric(quantile(sum_counts, 0.25)), median_sum = as.numeric(median(sum_counts)),
  q75_sum = as.numeric(quantile(sum_counts, 0.75)),
  q25_detected = as.numeric(quantile(detected_genes, 0.25)),
  median_detected = as.numeric(median(detected_genes)),
  q75_detected = as.numeric(quantile(detected_genes, 0.75)),
  q25_mito = as.numeric(quantile(mito_percent, 0.25)), median_mito = as.numeric(median(mito_percent)),
  q75_mito = as.numeric(quantile(mito_percent, 0.75))
), by = sample_id]
plot_stats <- merge(plot_stats, thresholds, by = "sample_id")
plot_long <- rbindlist(list(
  plot_stats[, .(sample_id, metric = "Total UMIs", q25 = q25_sum, median = median_sum, q75 = q75_sum,
                 lower = min_sum, upper = review_max_sum)],
  plot_stats[, .(sample_id, metric = "Detected genes", q25 = q25_detected, median = median_detected,
                 q75 = q75_detected, lower = min_detected, upper = review_max_detected)],
  plot_stats[, .(sample_id, metric = "Mitochondrial fraction (%)", q25 = q25_mito, median = median_mito,
                 q75 = q75_mito, lower = NA_real_, upper = max_mito_percent)]
))
plot_long[, sample_id := factor(sample_id, levels = paste0("sample", 1:40))]
p_qc <- ggplot(plot_long, aes(sample_id, median)) +
  geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 0.5, color = "#4C78A8") +
  geom_point(size = 1, color = "#1F4E79") +
  geom_point(aes(y = lower), shape = 25, fill = "#E45756", color = "#E45756", size = 1.4, na.rm = TRUE) +
  geom_point(aes(y = upper), shape = 24, fill = "#F2CF5B", color = "#9C7C00", size = 1.4, na.rm = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_y_continuous(trans = "log1p") +
  labs(x = NULL, y = "Observed distribution and data-derived thresholds (log1p scale)",
       title = "GSE183904 secondary QC across 40 author-filtered samples",
       subtitle = "Lines show IQR; circles show medians; triangles show lower/upper review thresholds") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

save_all <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = 300, device = agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         dpi = 600, device = agg_tiff, compression = "lzw", bg = "white")
}
save_all(p_qc, "FigS_GSE183904_secondary_QC", 12, 8.5)
fwrite(plot_long, file.path(source_dir, "FigS_GSE183904_secondary_QC.tsv"), sep = "\t", na = "")

audit <- c(
  "# GSE183904 manifest and secondary QC audit",
  "",
  sprintf("- Official matrices: 40 samples and %s author-retained cells.", format(sum(manifest$author_retained_cells), big.mark = ",")),
  sprintf("- Patient mapping: %d distinct manuscript patient IDs.", uniqueN(manifest$Patient_ID)),
  sprintf("- Primary gastric tumor replication set: %d independent tumors.", sum(manifest$include_primary_replication)),
  sprintf("- Peritoneal tumor sensitivity set: %d samples.", sum(manifest$specimen_class == "Peritoneal_Tumor")),
  sprintf("- Secondary low-quality failures: %s.", format(sum(qc_summary$secondary_low_quality_fail), big.mark = ",")),
  sprintf("- Secondary scDblFinder doublet flags: %s.", format(sum(qc_summary$secondary_scDblFinder_doublet), big.mark = ",")),
  sprintf("- Secondary-QC singlets retained: %s.", format(sum(qc_summary$final_secondary_qc_singlet), big.mark = ",")),
  "",
  "The GEO matrices already contain cells retained after the publication QC and DoubletFinder workflow. The present scDblFinder calls are secondary diagnostic flags run separately per sample; both author-retained and secondary-singlet analyses remain distinguishable.",
  "",
  "Only detected genes, total UMIs, and mitochondrial percentage drive the secondary low-quality rule. Upper UMI/gene outliers are review flags, while scDblFinder supplies the explicit doublet decision.",
  "",
  "Empty droplets were not submitted, so empty-droplet calling and ambient-RNA estimation cannot be reconstructed. This is a source-data limitation, not evidence that ambient RNA is absent.",
  "",
  "QuickGO GO:0030198 was not downloaded or accessed by this script."
)
writeLines(audit, file.path(audit_dir, "GSE183904_manifest_secondary_QC.md"), useBytes = TRUE)

session <- capture.output(sessionInfo())
writeLines(session, file.path(audit_dir, "GSE183904_secondary_QC_sessionInfo.txt"))
cat(paste(audit, collapse = "\n"), "\n")
