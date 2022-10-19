save_mbg_preds_parallel <- function(config, time_stamp, run_date, res_fit, df ,pathaddin="") {
  
  if(time_stamp==TRUE) output_dir <- <<<< FILEPATH REDACTED >>>>
  if(time_stamp==FALSE) output_dir <- <<<< FILEPATH REDACTED >>>>
  dir.create(output_dir, showWarnings = FALSE)
  
  # Save log of config file
  write.csv(config, <<<< FILEPATH REDACTED >>>>), row.names = FALSE)
  
  # save model
  save(res_fit, file = <<<< FILEPATH REDACTED >>>>)
  
  # save training data
  write.csv(
    df,
    file = <<<< FILEPATH REDACTED >>>>,
    row.names = FALSE
  )
  
  # Write a an empty file to indicate done with this parallel script
}
