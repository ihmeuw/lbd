save_mbg_preds_parallel <- function(config, time_stamp, run_date, res_fit, df ,pathaddin="") {
  
  if(time_stamp==TRUE) output_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
  if(time_stamp==FALSE) output_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
  dir.create(output_dir, showWarnings = FALSE)
  
  # Save log of config file
  write.csv(<<<< FILEPATH REDACTED >>>>, row.names = FALSE)
  
  # save model
  save(res_fit, file = (paste0(<<<< FILEPATH REDACTED >>>>)))

  # save training data
  write.csv(
    df,
    file = (paste0(<<<< FILEPATH REDACTED >>>>)),
    row.names = FALSE
  )
}
