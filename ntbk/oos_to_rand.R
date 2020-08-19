library(data.table)

ig <- 'wash'
load_from_indic <- 's_piped'
rd <- '2020_01_01_00_00_66'
data <- fread("<<<< FILEPATH REDACTED >>>>")
row_id_col <- 'row_id'
  data <- copy(data)
  n_before <- nrow(data)
  s_ho <- readRDS("<<<< FILEPATH REDACTED >>>>")

test1 <- sapply(1:length(s_ho), function(x) {
  min(dplyr::count(s_ho[[x]], fold)$n)
  })

test2 <- sapply(1:length(s_ho), function(x) {
  length(dplyr::count(s_ho[[x]], fold)$n)
  })

resamp <- which(test1 <= 28|test2 != 5)
names(s_ho)[resamp]


for (rr in resamp) {

  s_ho[[rr]]$fold <- sample(1:5, length(s_ho[[rr]]$fold), replace = TRUE)
}


test1 <- sapply(1:length(s_ho), function(x) {
  min(dplyr::count(s_ho[[x]], fold)$n)
  })

test2 <- sapply(1:length(s_ho), function(x) {
  length(dplyr::count(s_ho[[x]], fold)$n)
  })

resamp <- which(test1 <= 28|test2 != 5)
names(s_ho)[resamp]
