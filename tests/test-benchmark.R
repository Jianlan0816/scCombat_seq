library(testthat)
library(Seurat)

source(here::here("R", "core", "evaluation.R"))
source(here::here("R", "core", "benchmark_methods.R"))

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

test_that("benchmark_batch_methods returns correct metric fields for Harmony", {
  seu <- make_test_seurat()
  res <- benchmark_batch_methods(seu, returnHarmony = TRUE,
                                 returnLiger = FALSE, returnSeurat = FALSE,
                                 save_path = NULL)
  expect_named(res, "harmony")
  expect_named(res$harmony, c("seurat", "kBET_acceptance", "batch_LISI", "celltype_LISI"))
  expect_s4_class(res$harmony$seurat, "Seurat")
  expect_true(is.numeric(res$harmony$kBET_acceptance) || is.na(res$harmony$kBET_acceptance))
  expect_type(res$harmony$batch_LISI,    "double")
  expect_type(res$harmony$celltype_LISI, "double")
})

test_that("benchmark_batch_methods Seurat5 RPCA returns correct structure", {
  seu <- make_test_seurat()
  res <- benchmark_batch_methods(seu, returnHarmony = FALSE,
                                 returnLiger = FALSE, returnSeurat = TRUE,
                                 save_path = NULL)
  expect_named(res, "seurat5")
  expect_named(res$seurat5, c("seurat", "kBET_acceptance", "batch_LISI", "celltype_LISI"))
  expect_s4_class(res$seurat5$seurat, "Seurat")
})
