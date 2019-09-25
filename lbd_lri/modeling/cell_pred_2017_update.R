#copies over 2016 to make 2017 in a cell pred object

#define a function to copy over the most recent year from a cell pred object
cell_pred_year_copy <- function (cell_p, start_year, end_year, new_end_year, out_dir, label){
  total_rows <- nrow(cell_p)
  message(paste('total cell pred rows:', total_rows))
  n_years <- end_year - start_year + 1
  year_rows <- total_rows / n_years
  message(paste('cell pred rows in one year:', year_rows))
  start_index <- total_rows - year_rows + 1
  end_index <- total_rows
  cell_pred_update <- rbind(cell_p, cell_p[start_index:end_index,])
  saveRDS(cell_pred_update, '<<<< FILEPATH REDACTED >>>>')
}

#apply the function to the cell_pred objects

#cssa
load('<<<< FILEPATH REDACTED >>>>')
cell_pred_year_copy(cell_pred, 2000, 2016, 2017, '<<<< FILEPATH REDACTED >>>>', 'cssa_')

rm(cell_pred)

#sssa
load('<<<< FILEPATH REDACTED >>>>')
cell_pred_year_copy(cell_pred, 2000, 2016, 2017, '<<<< FILEPATH REDACTED >>>>', 'sssa_')

rm(cell_pred)

#essa
load('<<<< FILEPATH REDACTED >>>>')
cell_pred_year_copy(cell_pred, 2000, 2016, 2017, '<<<< FILEPATH REDACTED >>>>', 'essa2_')

rm(cell_pred)

#wssa
load('<<<< FILEPATH REDACTED >>>>')
cell_pred_year_copy(cell_pred, 2000, 2016, 2017, '<<<< FILEPATH REDACTED >>>>', 'wssa_')

rm(cell_pred)

#noaf
load('<<<< FILEPATH REDACTED >>>>')
cell_pred_year_copy(cell_pred, 2000, 2016, 2017, '<<<< FILEPATH REDACTED >>>>', 'name_')