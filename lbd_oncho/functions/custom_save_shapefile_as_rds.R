save_shapefile_as_rds <- function(shapefile_path,
                                  verbose = F,
                                  in_dir = <<<< FILEPATH REDACTED >>>>,
                                  out_dir = <<<< FILEPATH REDACTED >>>>) {
  
  # Load one function from stringr
  str_match <- stringr::str_match
  
  shape <- shapefile_path
  
  if(verbose == T) {
    message(paste0("     ", shape))
  }
  
  in_dir_no_slash <- gsub('/$','',in_dir)
  
  the_shape <- try(readOGR(dsn = in_dir_no_slash, layer = shape), silent = F)
  
  if (is.error(the_shape)) {
    # unlock_file(shape)
    return(the_shape)
  } else {
    saveRDS(the_shape, file = paste0(out_dir, shape, ".rds"))
    # unlock_file(shape)
    return("success")
  }
}
