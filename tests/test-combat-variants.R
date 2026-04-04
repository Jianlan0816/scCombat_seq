library(testthat)
library(Seurat)

source(here::here("R", "core", "evaluation.R"))
source(here::here("R", "core", "combat_variants.R"))

make_test_seurat <- function() {
  data("pbmc_small")
  pbmc_small$batch    <- rep(c("A", "B"), length.out = ncol(pbmc_small))
  pbmc_small$celltype <- sample(c("T", "B", "NK"), ncol(pbmc_small), replace = TRUE)
  pbmc_small <- NormalizeData(pbmc_small, verbose = FALSE) |>
    FindVariableFeatures(verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = 10, verbose = FALSE) |>
    RunUMAP(dims = 1:10, verbose = FALSE)
  pbmc_small
}

test_that("ComBat_combo returns correct structure for combat_seq", {
  seu <- make_test_seurat()
  res <- ComBat_combo(seu, combat_seq = TRUE, save_path = NULL)
  expect_named(res, "combat_seq")
  expect_s4_class(res$combat_seq$seurat, "Seurat")
  expect_type(res$combat_seq$kBET_acceptance, "double")
  expect_type(res$combat_seq$batch_LISI,      "double")
  expect_type(res$combat_seq$celltype_LISI,   "double")
})

test_that("ComBat_combo returns correct structure for combat_scseq", {
  seu <- make_test_seurat()
  res <- ComBat_combo(seu, combat_scseq = TRUE, save_path = NULL)
  expect_named(res, "combat_scseq")
  expect_s4_class(res$combat_scseq$seurat, "Seurat")
  expect_true(is.numeric(res$combat_scseq$kBET_acceptance) || is.na(res$combat_scseq$kBET_acceptance))
})

test_that("ComBat_combo combat_ind handles per-celltype correction", {
  seu <- make_test_seurat()
  res <- ComBat_combo(seu, combat_ind = TRUE, save_path = NULL)
  expect_named(res, "combat_ind")
  expect_s4_class(res$combat_ind$seurat, "Seurat")
})

test_that("ComBat_combo combat_pcseq returns corrected PCs without Seurat object", {
  seu <- make_test_seurat()
  res <- ComBat_combo(seu, combat_pcseq = TRUE, save_path = NULL)
  expect_named(res, "combat_pcseq")
  expect_null(res$combat_pcseq$seurat)
  expect_true(is.matrix(res$combat_pcseq$corrected_pcs))
})

test_that("summarise_metrics returns a data.frame with correct columns", {
  seu <- make_test_seurat()
  res <- ComBat_combo(seu, combat_seq = TRUE, save_path = NULL)
  df  <- summarise_metrics(res)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("method", "kBET_acceptance", "batch_LISI", "celltype_LISI") %in% colnames(df)))
})
