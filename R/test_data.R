# library(readr)
# library(usethis)
# ec_pvals <- read_delim("test_data/ec_pvals.tsv",
#                         delim = "\t", escape_double = FALSE,
#                         trim_ws = TRUE, show_col_types = FALSE)
# usethis::use_data(ec_pvals, overwrite = T)
#
# mtb_pvals <- read_delim("test_data/mtb_pvals.tsv",
#                         delim = "\t", escape_double = FALSE,
#                         trim_ws = TRUE, show_col_types = FALSE)
# usethis::use_data(mtb_pvals, overwrite = T)
#
# edges <- read_csv("global_network/kegg_mapformula_clean/enzyme_compound_edges.csv",
#                         trim_ws = TRUE, show_col_types = FALSE)
# usethis::use_data(edges, overwrite = T)

#' Writes example input files to a folder named "test_input", within the working
#'  directory.
#'
#' @export
write_test_files <- function() {
  require(readr)
  dir.create("test_input", showWarnings = F, recursive = T)
  write_tsv(ec_pvals, file = 'test_input/ec_pvals.tsv')
  write_tsv(mtb_pvals, file = 'test_input/mtb_pvals.tsv')
  write_csv(edges, file = 'test_input/enzyme_compound_edges_kegg.csv')
}
