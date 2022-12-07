umitools_to_mtx <- function(count_file,
                            output_path = NULL,
                            ...) {
  dat <- readr::read_tsv(count_file, ...)
  
  barcodes <- unique(dat$cell)
  genes <- unique(dat$gene)
  
  dat$gene_idx <- match(dat$gene, genes)
  dat$cell_idx <- match(dat$cell, barcodes)
  
  mat <- Matrix::sparseMatrix(
    i = dat$gene_idx,
    j = dat$cell_idx,
    x = dat$count
  )
  
  if (is.null(output_path)) {
    output_path <- fs::path_dir(count_file)
  }
  
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  Matrix::writeMM(mat, file.path(output_path, "matrix.mtx"))
  R.utils::gzip(file.path(output_path, "matrix.mtx"),
                overwrite = TRUE, remove = TRUE
  )
  
  readr::write_lines(genes, file.path(output_path, "features.tsv.gz"))
  readr::write_lines(barcodes, file.path(output_path, "barcodes.tsv.gz"))
}
umitools_to_mtx("~/Work/Runs/UT_R3/output/counts.tsv.gz", output_path = "~/Work/Runs/UT_R3/output")










