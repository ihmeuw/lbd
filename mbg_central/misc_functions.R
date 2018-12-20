## Convenience function to pull from our Github in R code.
  pullupstream<-function(){
    system('cd <<<< FILEPATH REDACTED >>>>/upstream/mbg\ngit pull upstream develop')
  }

# miscellaneous functions
  lstype<-function(type='closure'){
      inlist<-ls(.GlobalEnv)
      if (type=='function') type <-'closure'
      typelist<-sapply(sapply(inlist,get),typeof)
      return(names(typelist[typelist==type]))
  }
  
# to listen in on parallel processes in R: http://stackoverflow.com/questions/10903787/how-can-i-print-when-using-dopar
  Log <- function(text, ...) {
    msg <- sprintf(paste0(as.character(Sys.time()), ": ", text, "\n"), ...)
    cat(msg)
    write.socket(log.socket, msg)
  }

## Master function for launch script to weight for parallel submitted models to complete.
waitformodelstofinish <- function(sleeptime=100,
                                    path =  paste0(<<<< FILEPATH REDACTED >>>>, indicator_group, '/', indicator, '/output/', run_date),
                                    rd   = run_date,
                                    lv   = loopvars,
                                    showfiles = TRUE,
                                    showcluster = FALSE){

      n_finished <- length(grep("fin_", list.files(path)))

      lv <- data.table(lv)
      names(lv) <- c("reg", "holdout")
      lv$reg <- as.character(lv$reg)

      lv[, file := paste0(path, "/", "fin__bin0_", reg, "_", holdout)]

      while(n_finished != nrow(lv)){

            n_finished <- length(grep("fin_", list.files(path)))

            message('\n====================================================================================')
            message(sprintf('=====================      Run Date: %s      ======================',rd))
            message(paste0('\nAt ',Sys.time(),' .... ',n_finished,' Models have written output.'))
            if(showfiles){
              message('\nCurrently missing models:')
              for(i in 1:nrow(lv))
                if(file.exists(lv[i, file]) == F)
                  message(paste('Region =',lv[i,1],'| Holdout =',lv[i,2]))
            }
            n_cluster <- system("qstat | grep job_ | wc | awk '{print $1}'", intern = T)
            message(paste0('\nFuthermore, there are still ', n_cluster, ' jobs running on the cluster.'))
            if(showcluster){
              system("qstat -r | grep jobname | grep -oP \"(?<=job_)(.*)\"")
            }
            message('\n====================================================================================')
            message('====================================================================================')
            message("\n")
            Sys.sleep(sleeptime)
      }

      unlink(lv$file) # Clean up by deleting extra files once done with the loop
}

## Master function for launch script to wait for post-estimation jobs to complete.
waitforpostesttofinish <- function(sleeptime = 100,
                                   indicator = indicator,
                                   indicator_group = indicator_group,
                                   run_date = run_date,
                                   strata,
                                   showfiles = TRUE) {

  path <- paste0(<<<< FILEPATH REDACTED >>>>, indicator_group, '/', indicator, '/output/', run_date, "/temp_post_est/")
  
  lv <- data.table(strata)
  names(lv) <- "stratum"
  
  lv[, filename := paste0(path, stratum, "_post_est_list.RData")]
  
  n_finished <- 0
  n_needed <- length(strata)
  
   while (n_finished != n_needed) {
     lv[, exists := file.exists(filename)]
     n_finished <- nrow(lv[exists == T])
  
     message('\n====================================================================================')
     message('\nRunning post-estimation in parallel')
     message(paste0('Run Date: ', run_date))
     message(paste0('\nAt ',Sys.time(),' .... ',n_finished,' out of ', n_needed, ' strata have written output.'))
     if (showfiles) {
      message(paste0('\nCurrently missing strata:'))
      for(i in 1:nrow(lv[exists == F])) {
        message(subset(lv, exists == F)[i, stratum])
      }
     }
     message('\n====================================================================================')
     message('====================================================================================')
     Sys.sleep(sleeptime)
 }
}

## Master function for launch script to wait for results table jobs to complete.
waitforresultstable <- function(st = strata, 
                                rd = run_date, 
                                indic = indicator,
                                ig = indicator_group,
                                baseline_year = 2000,
                                measure = 'prevalence') {

  message("Waiting for results tables to finish...")
  r_left <- length(st)
  u_left <- length(st)
  while(length(r_left) > 0 & length(u_left) > 0) {
    message(paste0("\nCurrent time: ", Sys.time()))
    str_match <- stringr::str_match
    results_dir_r <- paste0('<<<< FILEPATH REDACTED >>>>', ig, '/',indic,'/output/', rd, '/table_', baseline_year, '/')
    results_dir_u <- paste0('<<<< FILEPATH REDACTED >>>>', ig, '/',indic,'/output/', rd, '/table_', baseline_year, '_unraked/')

    r_regs_done <- list.files(results_dir_r) %>% str_match(., paste0(measure, "_(.*).csv")) %>% .[,2] %>% .[!is.na(.)]
    u_regs_done <- list.files(results_dir_u) %>% str_match(., paste0(measure, "_(.*).csv")) %>% .[,2] %>% .[!is.na(.)]

    r_left <- st[!(st %in% r_regs_done)]
    u_left <- st[!(st %in% u_regs_done)]

    message(paste0(c("  Raked tables remaining: ", paste(r_left, collapse = ", "))))
    message(paste0(c("  Unraked tables remaining: ", paste(u_left, collapse = ", "))))

    Sys.sleep(60)

  }
}

## Master function for launch script to wait for admin aggregate jobs to complete.
waitforaggregation <- function(rd = run_date, 
                               indic = indicator,
                               ig = indicator_group,
                               ages,
                               regions, 
                               holdouts, 
                               raked,
                               dir_to_search = NULL) {

  # waitformodelstofinish() analog for aggregation
  if (is.null(dir_to_search)) {
    dir_to_search <- paste0("<<<< FILEPATH REDACTED >>>>",indicator_group,"/",indicator,"/output/",run_date,"/")
  }

  lv <- expand.grid(regions, holdouts, ages, raked) %>% as.data.table
  names(lv) <- c("reg", "holdout", "age", "raked")
  lv[, file := paste0(dir_to_search, "fin_agg_", reg, "_", holdout, "_", age, "_", raked)]
  n_left <- nrow(lv)

  message("Waiting for aggregation to finish...")
  
  while(n_left > 0) {
    message(paste0("\n=====================================",
                   "\nCurrent time: ", Sys.time()))
   
    lv[, file_exists := file.exists(file)]
    n_left <- nrow(lv[file_exists == F])
    lv_left <- lv[file_exists==F]

    if (n_left > 0) {
      message(paste0("Aggregation jobs remaining: ", n_left))
      for (i in 1:nrow(lv_left)) {
        message(paste0("  ", 
                       "Region: ", lv_left[i, reg], " | ",
                       "Holdout: ", lv_left[i, holdout], " | ", 
                       "Age: ", lv_left[i, age], " | ",
                       "Raked: ", lv_left[i, raked]))
      }
    }

    Sys.sleep(60)

  }
}

## Master function for launch script to combine aggregation results.
combine_aggregation <- function(rd = run_date, 
                                indic = indicator,
                                ig = indicator_group,
                                ages,
                                regions, 
                                holdouts, 
                                raked,
                                dir_to_search = NULL,
                                delete_region_files = T,
                                merge_hierarchy_list = F) {

  # Combine aggregation objects across region

  # Args:
  #   run_date, indicator, indicator_group for this run
  #   ages: single value or vector of ages
  #   regions: vector of regions used
  #   holdouts: vector of holdouts used, e.g. just 0 or c(1,2,3,4,5,0)
  #   raked: vector of raked values, e.g. just T, just F, or c(T,F)
  #   dir_to_search: which directory to search in (defaults to share directory)
  #   delete_region_files: logical. Should we delete the region-specific intermediate files?
  #   merge_hierarchy_list: logical. Do you want to merge the sp_hierarchy_list onto your admin tables?

  # Outputs: 
  #   rdata files for each combo of age/holdout/raked
  #   each with admin_0, admin_1, admin_2 data table objects & the sp_hierarchy_list object
  #   that maps them to names of admin units

  if (is.null(dir_to_search)) {
    dir_to_search <- paste0("<<<< FILEPATH REDACTED >>>>",ig,"/",indic,"/output/",run_date,"/")
  }

  message("Combining aggregation results...")

   for(rake in raked) {
    for (holdout in holdouts) {
      for (age in ages) {
        message(paste0("\nWorking on age: ", age, " | holdout: ", holdout, " | raked: ", rake))

        # Set up lists
        ad0 <- list()
        ad1 <- list()
        ad2 <- list()
        sp_h <- list()
      
        for (reg in regions) {
          message(paste0("  Region: ", reg))
          load(paste0(dir_to_search, indic, "_", ifelse(rake, "raked", "unraked"), 
                      "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"))

          if(merge_hierarchy_list == T) {
            # Prepare hierarchy list for adm0
            ad0_list <- subset(sp_hierarchy_list, select = c("ADM0_CODE", "ADM0_NAME", "region")) %>% unique

            # Prepare hierarchy list for adm1
            ad1_list <- subset(sp_hierarchy_list, 
                               select = c("ADM0_CODE", "ADM1_CODE", "ADM0_NAME", "ADM1_NAME", "region")) %>% 
                        unique

            # Merge
            admin_0 <- merge(ad0_list, admin_0, by = "ADM0_CODE", all.y = T)
            admin_1 <- merge(ad1_list, admin_1, by = "ADM1_CODE", all.y = T)
            admin_2 <- merge(sp_hierarchy_list, admin_2, by = "ADM2_CODE", all.y = T)
            rm(ad0_list, ad1_list)
          }

          ad0[[reg]] <- admin_0
          ad1[[reg]] <- admin_1
          ad2[[reg]] <- admin_2
          sp_h[[reg]] <- sp_hierarchy_list

          rm(admin_0, admin_1, admin_2, sp_hierarchy_list)
        }

        # Get to long format & save
        message("  Combining...")
        admin_0 <- rbindlist(ad0)
        admin_1 <- rbindlist(ad1)
        admin_2 <- rbindlist(ad2)
        sp_hierarchy_list <- rbindlist(sp_h)

        message("  Saving combined file...")
        save(admin_0, admin_1, admin_2, sp_hierarchy_list,
             file = paste0(dir_to_search, indic, "_", 
                           ifelse(rake, "raked", "unraked"), 
                           "_admin_draws_eb_bin", age, "_",  
                           holdout, ".RData"))
      }
    }
  }

  if (delete_region_files == T) {
    # Make sure all full files are written
    combos <- expand.grid(ifelse(raked, "raked", "unraked"), ages, holdouts)
    files_to_check <- sapply(1:nrow(combos), function(i) {
                          paste0(dir_to_search, indic, "_", combos[i, 1], 
                          "_admin_draws_eb_bin", combos[i, 2], "_", combos[i,3], ".RData")
                        })

    if (all(file.exists(files_to_check))) {
      message("All anticipated combined files were created successfully.  Deleting intermediate files...")
      combos <- expand.grid(ifelse(raked, "raked", "unraked"), ages, regions, holdouts)
      files_to_delete <- sapply(1:nrow(combos), function(i) {
                              paste0(dir_to_search, indic, "_", combos[i, 1], 
                              "_admin_draws_eb_bin", combos[i, 2], "_", combos[i,3], "_", combos[i, 4], ".RData")
                            })
      unlink(files_to_delete)
    } else {
      warning("Did not delete intermediate files - not all output files created successfully!")
    }
  }

  # Finally, delete the "fin" files
  fin_files_to_delete <- list.files(dir_to_search, pattern = "fin_agg_", full.names=T)
  unlink(fin_files_to_delete)
}

## Master function for launch script to create parallel job submission commands.
make_qsub <- function(user = Sys.info()['user'],
                      cores = slots,
                      memory = 100,
                      proj  = 'proj_geospatial',
                      ig = indicator_group,
                      repo = repo,
                      indic = indicator,
                      reg = "test",
                      age = 0,
                      vallevel = '',
                      hard_set_sbh_wgt = TRUE,
                      pct_sbh_wgt = 100,
                      rd = run_date,
                      log_location = 'sgeoutput',
                      code,
                      saveimage=FALSE,
                      test=FALSE,
                      holdout = 0,
                      use_base_covariates=FALSE,
                      test_pv_on_child_models=TRUE,
                      constrain_children_0_inf=FALSE,
                      child_cv_folds=5,
                      fit_stack_cv=TRUE,
                      additional_job_name = indicator, 
                      shell = "r_shell.sh",
                      modeltype = 'full',
                      geo_nodes = FALSE){

   if(test) t=1 else t=0
   if(use_base_covariates) b=1 else b=0
   if(test_pv_on_child_models) pvc=1 else pvc=0
   if(constrain_children_0_inf) cc0i=1 else cc0i=0
   if(fit_stack_cv) fscv=1 else fscv=0
   if(hard_set_sbh_wgt) hssw=1 else hssw=0

   if(saveimage==TRUE)  save.image(paste0('<<<< FILEPATH REDACTED >>>>', ig, '/', indic, '/model_image_history/pre_run_tempimage_', rd, '_bin',age,'_',reg,'_',holdout,'.RData'))

   dir.create(paste0('<<<< FILEPATH REDACTED >>>>', ig, '/', indic,'/output/',rd))

   sharedir <- sprintf('<<<< FILEPATH REDACTED >>>>%s/%s',ig,indic)

   if(log_location=='sgeoutput')
     logloc = sprintf('<<<< FILEPATH REDACTED >>>>/sgeoutput/%s',user)
   if(log_location=='sharedir'){
     logloc = sprintf('%s/output/%s',sharedir,rd)
     dir.create(sprintf('%s/errors',logloc))
     dir.create(sprintf('%s/output',logloc))
   }

  ## do we want to submit to geo nodes? if so, there are a few tweaks
  if(geo_nodes == TRUE){
    shell     <- "r_shell_geos.sh"   ## for the correct path to R
    proj      <- "proj_geo_nodes"    ## correct proj for geos nodes
    node.flag <- "-l geos_node=TRUE" ## send the job to geo nodes
  }else{
    node.flag <- "" ## don't do anything special
  }

  return(sprintf("qsub -e %s/errors -o %s/output -cwd -l mem_free=%iG -pe multi_slot %i -P %s %s -N job_%s_%s_%i_%i%s %s/%s %s/%s.R %s %i %s %i %i %s %s %i %i %i %i %i %s %s %i %i",
                 logloc,logloc,memory,cores,proj,node.flag,additional_job_name, reg,age,holdout,vallevel,ig,shell,ig,code,reg,age,rd,t,holdout,indic,ig,b,pvc,cc0i,child_cv_folds,fscv,modeltype,vallevel,pct_sbh_wgt,hssw))

}

## Alternative function for launch script to create parallel job submission commands, depending on computing environment.
make_qsub_share <- function(user   = Sys.info()['user'],
                          cores  = slots,
                          memory = 100,
                          proj   = 'proj_geospatial',
                          ig     = indicator_group,
                          coderepo   = repo,
                          indic  = indicator,
                          reg    = "test",
                          age    = 0,
                          rd     = run_date,
                          log_location = 'sharedir',
                          code   = NULL,
                          saveimage = FALSE,
                          test      = FALSE,
                          additional_job_name = indicator, 
                          holdout   = 0,
                          geo_nodes = FALSE){

   # save an image
   if(saveimage==TRUE)  save.image(paste0('<<<< FILEPATH REDACTED >>>>', ig, '/', indic, '/model_image_history/pre_run_tempimage_', rd, '_bin',age,'_',reg,'_',holdout,'.RData'))

   # sort directories
   sharedir <- sprintf('<<<< FILEPATH REDACTED >>>>%s/%s',ig,indic)
   dir.create(sprintf('%s/output/%s',sharedir,rd), showWarnings = FALSE)

   if(log_location=='sgeoutput')
     logloc = sprintf('<<<< FILEPATH REDACTED >>>>/sgeoutput/%s',user)
   if(log_location=='sharedir'){
     logloc = sprintf('%s/output/%s',sharedir,rd)
     dir.create(sprintf('%s/errors',logloc), showWarnings = FALSE)
     dir.create(sprintf('%s/output',logloc), showWarnings = FALSE)
   }

  ## do we want to submit to geo nodes? if so, there are a few tweaks
  if(geo_nodes == TRUE){
    shell     <- "r_shell_geos.sh"   ## for the correct path to R
    proj      <- "proj_geo_nodes"    ## correct proj for geos nodes
    node.flag <- "-l geos_node=TRUE" ## send the job to geo nodes
    shell     <- sprintf('%s/mbg_central/share_scripts/shell_geos.sh',coderepo)
  }else{
    proj      <- "proj_geospatial"
    node.flag <- ""
    shell     <- sprintf('%s/mbg_central/share_scripts/shell_prod.sh',coderepo)
  }

  # set code to shared paralel if its NULL
  if(is.null(code)) {
    code <- sprintf('%s/mbg_central/share_scripts/parallel_model.R',coderepo)
  } else {
    code <- sprintf('%s/%s/%s.R',coderepo,ig,code)
  }

  return(sprintf("qsub -e %s/errors -o %s/output -cwd -l mem_free=%iG -pe multi_slot %i -P %s %s -N job_%s_%s_%i_%i %s %s %s %i %s %i %i %s %s fin",logloc,logloc,memory,cores,proj,node.flag,additional_job_name,reg,age,holdout,shell,code,reg,age,rd,test,holdout,indic,ig))

}

## Master function for launch script to create parallel job submission commands for post-estimation.
make_qsub_postest <- function(user = Sys.info()['user'],
                              cores = slots,
                              memory = 100,
                              proj  = 'proj_geospatial',
                              ig = indicator_group,
                              repo = repo,
                              indic = indicator,
                              stratum = "test",
                              rd = run_date,
                              log_location = 'sgeoutput',
                              code,
                              script_dir = 'mbg_central/share_scripts/',
                              keepimage=FALSE,
                              shell = "r_shell.sh",
                              geo_nodes = FALSE){

    sharedir <- paste0('<<<< FILEPATH REDACTED >>>>', ig, '/', indic, '/output/', rd, '/')
    temp_dir <- paste0(sharedir, 'temp_post_est/')
    dir.create(temp_dir, showWarnings = F)

   if(log_location=='sgeoutput')
     logloc = sprintf('<<<< FILEPATH REDACTED >>>>/sgeoutput/%s',user)
   if(log_location=='sharedir'){
     logloc = sharedir
     dir.create(sprintf('%serrors',logloc), showWarnings = F)
     dir.create(sprintf('%soutput',logloc), showWarnings = F)
   }

  ## do we want to submit to geo nodes? if so, there are a few tweaks
  node.flag <- ifelse(geo_nodes == FALSE, "", "-l geos_node=TRUE")
  if(geo_nodes == TRUE){ shell <- "r_shell_geos.sh" }

  gn <- as.character(geo_nodes)

  # last line below sets args for parallel script
  return(paste0("qsub ",
                "-e ", logloc, "errors ",
                "-o ", logloc, "output ",
                "-cwd -l mem_free=", memory, "G ",
                "-pe multi_slot ", cores, " ",
                "-P ", proj, " ", node.flag, " ",
                "-N job_pe_", indic, "_", stratum, " ",
                ig, "/", shell, " ",
                script_dir, "/", code, ".R ",
                stratum, " ", rd, " ", indic, " ", ig, " ", gn))
  
}

## Function to get non-missing cell IDs from a raster object
cellIdx <- function (x) which(!is.na(getValues(x[[1]])))

## Function to create a raster with a template and matrix of same dimensions
insertRaster <- function (raster, new_vals, idx = NULL) {

  # calculate cell index if not provided
  if (is.null(idx)) idx <- cellIdx(raster)

  # check the index makes superficial sense
  stopifnot(length(idx) == nrow(new_vals))
  stopifnot(max(idx) <= ncell(raster))

  # create results raster
  n <- ncol(new_vals)
  raster_new <- raster::brick(replicate(n,
                                        raster[[1]],
                                        simplify = FALSE))
  names(raster_new) <- colnames(new_vals)

  # update the values
  for(i in 1:n) {
    raster_new[[i]][idx] <- new_vals[, i]
  }

  return (raster_new)

}

condSim <- function (vals, weights = NULL, group = NULL, fun = NULL, ...) {
  # given a matrix of pixel-level prevalence samples `prev`
  # where each rows are pixels and columns are draws, a vector
  # of corresponding pixel populations `pop`, and an optional pixel
  # grouping factor `group`, return draws for the total deaths in each
  # group, or overall if groups are not specified

  # get dimensions of vals
  ncell <- nrow(vals)
  ndraw <- ncol(vals)

  # capture function as a string
  fun_string <- deparse(substitute(fun))

  # check fun accepts a

  # check dimensions of weights and group, set to 1 if not specified
  if (is.null(weights)) {
    weights <- rep(1, ncell)
  } else {
    if (length(weights) != ncell) {
      stop (sprintf('number of elements in weights (%i) not equal to number of cells in vals (%i)',
                    length(weights),
                    ncell))
    }
  }

  if (is.null(group)) {
    group <- rep(1, length(weights))
  } else {
    if (length(group) != ncell) {
      stop (sprintf('number of elements in group (%i) not equal to number of cells in vals (%i)',
                    length(group),
                    ncell))
    }
  }

  # otherwise, get the levels in group and create a matrix of results
  levels <- unique(na.omit(group))
  nlevel <- length(levels)

  ans <- matrix(NA,
                ncol = ndraw,
                nrow = nlevel)
  rownames(ans) <- levels

  # loop through levels in group, getting the results
  for (lvl in 1:nlevel) {

    # get an index o pixels in the level
    idx <- which(group == levels[lvl])

    # by default, calculate a weighted sum
    if (is.null(fun)) {

      # get draws and add to results
      # exception for if area has 1 cell, transpose matrix so it conforms (RB)
      if(all(dim(t(vals[idx, ]))==c(1,ndraw))){
        ans[lvl, ] <- weights[idx] %*% t(vals[idx, ])
      } else {
        ans[lvl, ] <- weights[idx] %*% vals[idx, ]
      }

    } else {

      # otherwise, apply function to each column
      ans[lvl, ] <- apply(vals[idx, ], 2, fun, weights = weights[idx], ...)

    }

  }

  # if only one level, make this a vector
  if (nlevel == 1) ans <- as.vector(ans)

  # return result
  return (ans)

}

defaultOptions <- function (resolution = 5,
                            location = loc(),
                            cores = switch(loc(), oxford = 60, seattle = 60),
                            spacetime = TRUE,
                            start = Sys.time()) {

  # built in check for if they selected more than one country

  # if these options are not already set, set them at these default values

  # list all visible objects in this environment
  object_names <- ls()

  # get them in a named list
  objects <- lapply(object_names,
                    get,
                    envir = environment())
  names(objects) <- object_names

  # empty vector of the options to keep
  keep <- c()

  # loop through identifying those that don't already exist
  for (i in 1:length(objects)) {

    # if doesn't exist yet...
    if (is.null(options(names(objects)[i])[[1]])) {

      # add it to the list
      keep <- c(keep, i)

    }
  }

  # keep the new options
  objects <- objects[keep]

  if (length(objects) > 0) {
    # notify the user
    message('\nThe following options were not defined and have been set to their defaults:\n')
    for(i in 1:length(objects)) {
      message(sprintf('    %s = %s',
                      names(objects)[i],
                      paste(objects[[i]], collapse = ', ')))
    }
    message('\n')
  }

  # add to options
  options(objects)
}


# for cleaning clutter
"rmlike" <- function(...) {
  names <- sapply(
    match.call(expand.dots = FALSE)$..., as.character)
  names = paste(names,collapse="|")
  Vars <- ls(1)
  r <- Vars[grep(paste("^(",names,").*",sep=""),Vars)]
  rm(list=r,pos=1)
}

reportTime <- function () {

  # get start time stored in options
  start <- options()$start

  # if it exists
  if (!is.null(start)) {

    # get elapsed time
    now <- Sys.time()
    diff <- difftime(now, start)

    # format nicely and report
    diff_string <- capture.output(print(diff, digits = 3))
    diff_string <- gsub('Time difference of', 'Time elapsed:', diff_string)
    message (diff_string)
  }
}

tidyInstall <- function (packages, versions = NULL, lib = .libPaths()[1]) {
  # given a vector of package names and optional vector of version numbers
  # install them from CRAN if they aren't already loaded

  npkg <- length(packages)

  # check the arguments
  stopifnot(inherits(packages, 'character'))

  if (!is.null(versions)) {
    stopifnot(length(versions) == npkg)
    stopifnot(inherits(versions, 'character'))
  }

  # run a loop if multiple required
  if (npkg > 1) {

    for (i in 1:npkg) {

      if (is.null(versions))
        versions_tmp <- NULL
      else
        versions_tmp <- versions[i]

      tidyInstall(packages[i], versions_tmp, lib)
    }

  } else {
    # otherwise just for one

    # first make sure versions and devtools are installed
    for (pkg in c('versions', 'devtools')) {
      if (!(pkg %in% installed.packages(lib.loc = lib)))
        install.packages(pkg, lib = lib)
    }

    # check whether the package is installed
    installed <- packages %in% rownames(installed.packages(lib.loc = lib))

    if (installed && is.null(versions)) {

      # if it is installed the user doesn't care about the version exit
      return (invisible())

    } else {

      # otherwise, do some installing...

      # if it's installed, check the version
      if (installed) {

        installed_version <- versions::installed.versions(packages, lib = lib)

        # if that's the required version, exit
        if (versions == installed_version) {
          return (invisible())
        } else {
          # otherwise, install the correct version
          devtools::with_lib(new = lib,
                             devtools::install_version(packages, versions))
        }

      } else {
        # if it isn't installed ...

        if (!is.null(versions)) {
          # if a version is required, install it and exit
          devtools::with_lib(new = lib,
                             devtools::install_version(packages, versions))
          return (invisible())
        } else {
          # otherwise get the most recent version and exit
          install.packages(packages, lib = lib)
          return (invisible())

        }
      }
    }
  }
}

subtractYears <- function (Date, years) {
  # given a vector of dates in Date format,
  # and a number of years (positive numeric, can be non-integer)
  # subtract the required number of years and return

  stopifnot(years >= 0)
  stopifnot(inherits(Date, 'Date'))

  n_Date <- length(Date)
  n_years <- length(years)

  if (n_years != n_Date) {
    if (n_years == 1) {
      years <- rep(years, n_Date)
    } else {
      stop('Date and years have different lengths')
    }
  }

  Date <- as.POSIXlt(Date)

  month <- round(12 * years %% 1)
  year <- floor(years)

  Date$year <- Date$year - year
  Date$mon <- Date$mon - month

  Date <- as.Date(Date)

  return (Date)

}

toYear <- function(Date) format(Date, '%Y')

targetYear <- function (period, period_end, width = 60) {
  # get the target year, based on the period, period size (in months) and end
  # date of the final period

  # convert date to cmc
  period_end <- Date2cmc(period_end)

  # get number of months preceeding
  months <- width / 2 + width * (period - 1)

  # subtract
  cmc <- period_end - months + 1

  # format as a year
  ans <- toYear(cmc2Date(cmc))

  return (ans)

}

# convert strings or numerics to a sequential numeric index
# pass any number of vectors of the same length, get a numeric ID for the
# unique ones
idx <- function (...) {
  list <- list(...)

  # get their lengths
  lengths <- lapply(list, length)
  if (!isTRUE(do.call(all.equal, lengths))) {
    stop('vectors do not appear to have the same length')
  }
  # combine them
  x <- do.call(paste, list)
  #get numeric index
  match(x, unique(x))
}

getPaths <- function(pathfile = <<<< FILEPATH REDACTED >>>>) {
  # get file paths for the key datasets
  # path points to a csv file containing named paths, the function returns a dataframe
  # of named filepaths
  paths <- read.csv(pathfile,
                    row.names = 1)
  data.frame(t(paths),
             stringsAsFactors = FALSE)
}


# given the admin level, a GAUL code and the admin1 and 2 shapefiles,
# return an SPDF with the relevant area
getPoly <- function (level, code, admin1, admin2) {

  # get admin level
  if (level == 1) {
    admin <- admin1
  } else {
    admin <- admin2
  }

  # find the reight area
  idx <- match(code, admin$GAUL_CODE)

  # if it's valid
  if (length(idx) == 1) {
    ans <- admin[idx, ]
  } else {
    # handle errors
    warning (paste0("something's wrong on row ", i))
    ans <- NULL
  }

  # return result
  return (ans)
}

# functions to determine proportion of women in the given age group
matAgeRate <- function (age_group,
                        groups = c('15-19', '20-24', '25-29', '30-34')) {
  # given a vector `age_group` reporting the age group to which mothers belong,
  # return the proportion falling in each of the age groups in `groups`.

  # check it's a character vector
  stopifnot(class(age_group) == 'character')

  # count the number in all age groups
  counts <- table(age_group)

  # add on any groups not represented a 0s
  missing_groups <- which(!(groups %in% names(counts)))

  if (length(missing_groups) > 0) {
    dummy_counts <- rep(0, length(missing_groups))
    names(dummy_counts) <- groups[missing_groups]
    counts <- c(counts, dummy_counts)
  }

  # get proportions
  props <- counts / sum(counts)

  # find the ones we want
  idx_keep <- match(groups, names(props))

  # and return these
  return (props[idx_keep])

}


# functions to determine proportion of women in the given age group
matAgeParRate <- function (ceb,
                           age_group,
                           groups = c('15-19', '20-24', '25-29', '30-34')) {
  # given vectors `ceb` and `age_group` reporting the number of children ever
  # born to mothers and the age groups to which they belong,
  # return the proportion of births falling in each of the age groups in
  # `groups`.

  # check it's a character vector
  stopifnot(class(ceb) %in% c('numeric', 'integer'))
  stopifnot(class(age_group) == 'character')

  # count the number of births in all age groups
  counts <- tapply(ceb, age_group, sum)

  # add on any groups not represented a 0s
  missing_groups <- which(!(groups %in% names(counts)))

  if (length(missing_groups) > 0) {
    dummy_counts <- rep(0, length(missing_groups))
    names(dummy_counts) <- groups[missing_groups]
    counts <- c(counts, dummy_counts)
  }

  # get proportions
  props <- counts / sum(counts)

  # find the ones we want
  idx_keep <- match(groups, names(props))

  # and return these
  return (props[idx_keep])

}

# function to tabulate maternal age rate
tabMatAgeRate <- function (age_group,
                           cluster_id,
                           groups = c('15-19', '20-24', '25-29', '30-34')){
  # given a vector `age_group` reporting the age group to which mothers belong,
  # and a vector `cluster_id` giving the cluster to which each mother belongs
  # return a matrix - with number of rows equal to the number of unique elements
  # in `cluster_id` and number of columns equal to the length of `groups` -
  # giving the proportion falling in each of the age groups in `groups`.

  # get the grouped data as a list
  ans <- tapply(age_group,
                cluster_id,
                matAgeRate,
                groups)

  # combine into a matrix
  ans <- do.call(rbind, ans)

  # and return
  return (ans)

}

# function to tabulate maternal age rate
tabMatAgeParRate <- function (ceb,
                              age_group,
                              cluster_id,
                              groups = c('15-19', '20-24', '25-29', '30-34')){
  # given vectors `ceb` and `age_group` reporting the number of children ever
  # born to mothers and the age groups to which they belong,
  # return a matrix - with number of rows equal to the number of unique elements
  # in `cluster_id` and number of columns equal to the length of `groups` -
  # giving the proportion of births falling in each of the age groups in
  # `groups`.

  # get unique cluster ids
  cluster_ids <- unique(cluster_id)

  # create dummy results matrix
  ans <- matrix(NA,
                nrow = length(cluster_ids),
                ncol = length(groups))

  rownames(ans) <- cluster_ids
  colnames(ans) <- groups

  # loop through the cluster ids calculating the rates
  for (i in 1:length(cluster_ids)) {

    # get index for the cluster
    idx_cluster <- which(cluster_id == cluster_ids[i])

    ans[i, ] <- matAgeParRate(ceb = ceb[idx_cluster],
                              age_group = age_group[idx_cluster],
                              groups = groups)

  }

  # and return
  return (ans)

}

# get predictions from predictive INLA glm models
predGLM <- function (result,
                     intercept = '(Intercept)',
                     fixed_continuous = NULL,
                     fixed_group = NULL,
                     random_group = 'cluster_id') {

  # starting means and variances
  pred_mean <- pred_var <- 0

  # add intercept terms
  if (!is.null(intercept)) {
    pred_mean <- pred_mean + result$summary.fixed[intercept, 'mean']
    pred_var <- pred_var + result$summary.fixed[intercept, 'sd'] ^ 2
  }

  # add continuous fixed effects terms
  if (!is.null(fixed_continuous)) {
    for (name in fixed_continuous) {

      # get coefficients
      coef_mean <- result$summary.fixed[name, 'mean']
      coef_var <- result$summary.fixed[name, 'sd'] ^ 2

      # get covariate
      cov <- result$model.matrix[, name]

      pred_mean <- pred_mean + cov * coef_mean
      pred_var <- pred_var + (cov ^ 2) * coef_var

    }
  }

  # add discrete fixed effects terms
  if (!is.null(fixed_group)) {
    for (name in fixed_group) {

      # find group members
      idx <- grep(sprintf('^%s*', name), rownames(result$summary.fixed))

      # get clean names
      names <- rownames(result$summary.fixed)[idx]
      names_clean <- gsub(name, '', names)

      # get coefficients
      coef_mean <- result$summary.fixed[idx, 'mean']
      coef_var <- result$summary.fixed[idx, 'sd'] ^ 2

      # get covariate
      cov <- result$model.matrix[, names]

      pred_mean <- pred_mean + as.vector(cov %*% coef_mean)
      pred_var <- pred_var + as.vector((cov ^ 2) %*% coef_var)

    }
  }

  # add discrete random effects terms
  if (!is.null(random_group)) {
    for (name in random_group) {

      # get coefficients
      coef_level <- result$summary.random[[name]][, 'ID']
      coef_mean <- result$summary.random[[name]][, 'mean']
      coef_var <- result$summary.random[[name]][, 'sd'] ^ 2

      # get covariate
      cov <- result$.args$data[, name]

      # match up the levels
      length(cov)
      length(coef_level)
      idx <- match(cov, coef_level)

      pred_mean <- pred_mean + coef_mean[idx]
      pred_var <- pred_var + coef_var[idx]

    }
  }

  # return result
  ans <- data.frame(mean = pred_mean,
                    sd = sqrt(pred_var))

  return (ans)

}

getCols <- function (df, period = 1) {
  # subset results matrix
  df[, grep(sprintf('^p%s_', period),
            colnames(df))]
}

rnormMatrix <- function (n, mean, sd) {
  # sample random normals with matrix parameters
  # returna an array as a result

  # coerce to matrix
  mean <- as.matrix(mean)
  sd <- as.matrix(sd)

  # get & check dimensions
  ncol <- ncol(mean)
  nrow <- nrow(mean)
  stopifnot(ncol(sd) == ncol)
  stopifnot(nrow(sd) == nrow)

  # convert to vector
  mean <- as.vector(mean)
  sd <- as.vector(sd)

  # sample
  draws <- rnorm(n * length(mean), mean, sd)

  # reshape
  ans <- array(draws, dim = c(nrow, ncol, n))

  return (ans)
}

accumulate <- function (mean, sd, months = c(1, 11, 24, 24), nsample = 1000) {
  # given matrices of means and standard deviations of
  # logit component death probabilities, each column giving
  # consecutive and adjacent age bins and rows giving the records,
  # calculate the logit of the cumulative mortality probability across
  # the age bins by Monte Carlo sampling

  # generate random draws from these logits
  draws <- rnormMatrix(nsample, mean, sd)

  # convert to draws of survival probabilities
  draws_p <- 1 - plogis(draws)

  # raise to power of number of months per bin
  for (i in 1:dim(draws_p)[2]) {
    draws_p[, i, ] <- draws_p[, i, ] ^ months[i]
  }

  # loop through age bins accumulating them
  for (i in 2:dim(draws_p)[2]) {
    draws_p[, i, ] <- draws_p[, i, ] * draws_p[, i - 1, ]
  }

  # convert back to logit mortality probabilities
  draws_l <- qlogis(1 - draws_p)

  # calculate the means and standard deviations of these logits
  l_mean <- apply(draws_l, c(1, 2), mean)
  l_sd <- apply(draws_l, c(1, 2), sd)

  # return as a list
  return (list(mean = l_mean,
               sd = l_sd))

}


# functions to centre and scale matrices, columnwise

getCentreScale <- function (x, exclude = NULL, na.rm = TRUE) {
  # get dataframe of centreing and scaling values to convert x
  # to the standard normal. exclude is an optional character vector
  # giving column names to exclude from scaling

  # get means and SDs for all columns
  df <- data.frame(name = colnames(x),
                   mean = colMeans(x, na.rm = na.rm),
                   sd = apply(x, 2, sd, na.rm = na.rm))
  rownames(df) <- NULL

  # replace any zero standard deviations with 1
  # to avoid divide-by-zero errors
  df$sd[df$sd == 0] <- 1

  # if any named covariates are to be excluded, set mean to 0 and sd to 1
  if (!is.null(exclude)) {
    idx <- match(exclude, df$name)
    df$mean[idx] <- 0
    df$sd[idx] <- 1
  }

  return (df)
}

centreScale <- function (x, df, inverse = FALSE) {
  # apply pre-calculated centreing/scaling to matrix x,
  # with fixed dataframe of means/sds df
  # or uncentre/unscale if inverse = TRUE

  # get the centreing/scaling dataframe if not available
  if (is.null(df))
    df <- getCentreScale(x)

  # get index to match up values with column names
  names <- colnames(x)
  idx <- match(names, df$name)

  if (any(is.na(idx))) {
    stop ('could not match up column names with the values in df')
  }

  df <- df[idx, ]

  # apply transformations
  if (!inverse) {
    # move to standard normal

    # centre
    x <- sweep(x, MARGIN = 2, STATS = df$mean, FUN = '-')
    # scale
    x <- sweep(x, MARGIN = 2, STATS = df$sd, FUN = '/')

  } else {
    # inverse case, move from standard normal to original

    # unscale
    x <- sweep(x, MARGIN = 2, STATS = df$sd, FUN = '*')
    # uncentre
    x <- sweep(x, MARGIN = 2, STATS = df$mean, FUN = '+')

  }

  return (x)

}

loc <- function () {
  # detect location based on username

  user <- system('echo $USER',
                 intern = TRUE)

  # use this to work out our location
  location <- switch(user,
                     lina1864 = 'oxford',
                     goldingn = 'seattle',
                     royburst = 'seattle')

  return (location)

}

startLog <- function(file = 'full_run.log') {
  # create a logfile and start logging
  con <- file(file)
  sink(con,
       split = TRUE)
  sink(con,
       type = 'message')

  # report session info
  message(sprintf('# starting log at %s\n', Sys.time()))
  message('# session info:\n')
  print(sessionInfo())
  message('\n# run log:\n')
}

stopLog <- function() {
  # report session info
  message('\n# session info:\n')
  print(sessionInfo())
  message(sprintf('\n# stopping log at %s', Sys.time()))
  # stop logging to the logfile
  sink()
  sink(type = 'message')
}

list2paths <- function (x) {

  # get node and branch names
  nodes <- unlist(x, use.names = FALSE)
  branches <- unlist(getBranches(x),
                     use.names = FALSE)

  # append the node name
  paths <- file.path(branches, nodes)

  # remove extra slashes
  paths <- gsub('/+', '/', paths)

  paths
}

getBranches <- function(x, parent = "") {
  # loop through levels of a tree-like list,
  # flattening names of levels into
  # slash-separated character vector

  # get element names
  names <- names(x)

  # if unnamed, repeat the parent names
  if (is.null(names)) {
    names <- rep.int(parent, length(x))
  } else{
    names <- paste0(parent, names)
  }

  # if there are more branches ahead
  if (is.list(x)) {

    # add a slash
    names <- paste0(names, "/")

    # define a function to loop through them
    getTwigs <- function(i) {
      getBranches(x[[i]], names[i])
    }

    # get the names from these
    names <- lapply(1:length(x), getTwigs)

  }

  return (names)

}

db <- function (machine = c('local'), user = 'lina1864', warnings = FALSE) {
  # check dropbox syncing status, optionally on multiple remote machines
  # requies key-pairs and firewall acces tho...
  # e.g.:
  #   synced <- db(c('local',
  #                  'zoo-mapbig02.zoo.ox.ac.uk',
  #                  'big04.well.ox.ac.uk',
  #                  'big05.well.ox.ac.uk',
  #                  'big06.well.ox.ac.uk',
  #                  'zoo-seeg02.zoo.ox.ac.uk')))

  # vectorize by recursion
  if (length(machine) > 1) {

    ans <- sapply(machine, db)

  } else {

    # for local machine
    if (machine == 'local') {

      call <- 'dropbox status'

    } else {

      call <- sprintf("ssh %s@%s 'dropbox status'",
                      user,
                      machine)

    }

    # run the system call
    result <- tryCatch({
      system(call,
             intern = TRUE,
             ignore.stderr = TRUE)
    },
    error = function(e) 'nope')

    # if the user wants warnings
    if (warnings) {

      # warn if dropbox is switched off
      if (result == "Dropbox isn't running!") {
        warning (sprintf('dropbox not running on %s',
                         machine))
      }

      # warn if the system call errored
      if (result == "nope") {
        warning (sprintf('some sort of error occurred',
                         machine))
      }

    }

    # check if it's up to date
    ans <- result == 'Up to date'
  }

  # return the result
  return (ans)

}

unpack <- function (tmp) {
  # unpack a nested list where the outer list is unnamed, but each inner
  # list contains a single named element.
  # Assign this element to the parent environment and delete the list
  # it is called on fromt he parent ennvironment

  # get name of list in calling environment
  tmp_name <- deparse(substitute(tmp))

  # get calling environment
  pf <- parent.frame()

  # unpack into a single list
  tmp <- unlist(tmp, recursive = FALSE)

  # loop through assigning to calling environment
  for (i in 1:length(tmp)) {
    assign(names(tmp)[i], tmp[[i]], envir = pf)
  }

  # remove object from calling environment
  rm(list = tmp_name, envir = pf)

  # return nothing
  return (invisible(NULL))
}

prepLines <- function (rate,
                       year,
                       ylab = '',
                       title = '',
                       xlim = c(1999, 2016),
                       ylim = c(0, max(rate)),
                       line_years = c(2000, 2005, 2010, 2015)) {
  # set up the base plot for the custom line chart

  plot(rate ~ year,
       type = 'n',
       ylab = '',
       xlab = '',
       axes = FALSE,
       xlim = xlim,
       ylim = ylim)

  for (ly in line_years) {
    lines(x = rep(ly, 2),
          y = ylim,
          lwd = 3,
          lty = 3,
          col = grey(0.8))
  }

  axis(1,
       lty = 0,
       col.axis= grey(0.4),
       line = -1)

  axis(2,
       las = 2,
       cex.axis = 0.8,
       col = grey(0.4),
       col.axis= grey(0.4))

  title(main = title,
        col.main = grey(0.35),
        cex.main = 1.2,
        line = 0.5)

  title(ylab = ylab,
        col.lab = grey(0.4),
        cex.lab = 1.2)
}

addLines <- function (rate,
                      year,
                      country,
                      countries = NULL,
                      col = grey(0.7),
                      size = 1,
                      border = grey(0.4)) {

  # given vectors: rate (y axis), year (x axis)
  # and country (grouping factor), make a nice line
  # plot with lollipop-like likes with border colour
  # 'border' and fill colour 'col' (repeated if length one)
  # If 'countries' is NULL, all countries are plotted,
  # otherwise only those named in this character vector

  # check inputs
  stopifnot(all.equal(length(rate),
                      length(year),
                      length(country)))


  # sort countries
  all_countries <- sort(unique(country))
  if (is.null(countries)) {
    countries <- all_countries
  } else {
    stopifnot(all(countries %in% all_countries))
  }
  n_ctry <- length(countries)

  # expand col and bg if needed
  if (length(col) == 1) {
    col <- rep(col, n_ctry)
  } else {
    stopifnot(length(col) == n_ctry)
  }

  # loop through countries
  for (i in 1:n_ctry) {

    ctry <- countries[i]

    idx_ctry <- which(country == ctry)

    # dark grey outline
    lines(rate[idx_ctry] ~ year[idx_ctry],
          col = border,
          lwd = 7.5 * size)
    points(rate[idx_ctry] ~ year[idx_ctry],
           col = border,
           cex = 1 * size,
           pch = 16)

    # coloured foreground
    lines(rate[idx_ctry] ~ year[idx_ctry],
          col = col[i],
          lwd = 6 * size)

    points(rate[idx_ctry] ~ year[idx_ctry],
           col = col[i],
           pch = 16,
           cex = 0.85 * size)

  }
}

addLabels <- function (rate,
                       year,
                       country,
                       countries = NULL,
                       col = grey(0.7),
                       gap = diff(range(rate)) / 60,
                       cex = 0.7,
                       adj = 0,
                       xpd = NA,
                       ...) {

  # add country names on RHS
  # arguments as before, with gap to define spacing between labels.
  # dots are passed to text
  require (plotrix)

  # check inputs
  stopifnot(all.equal(length(rate),
                      length(year),
                      length(country)))

  # sort countries
  all_countries <- sort(unique(country))
  if (is.null(countries)) {
    countries <- all_countries
  } else {
    stopifnot(all(countries %in% all_countries))
  }
  n_ctry <- length(countries)

  # expand col and bg if needed
  if (length(col) == 1) {
    col <- rep(col, n_ctry)
  } else {
    stopifnot(length(col) == n_ctry)
  }

  # keep only those for latest year
  max_year <- max(as.numeric(year))
  year_idx <- which(year == max_year)
  rate <- rate[year_idx]
  year <- year[year_idx]
  country <- country[year_idx]

  # keep only those for countries requires
  ctry_idx <- which(country %in% countries)
  rate <- rate[ctry_idx]
  year <- year[ctry_idx]
  country <- country[ctry_idx]

  # order by countries
  ctry_o <- match(countries, country)
  rate <- rate[ctry_o]
  year <- year[ctry_o]
  country <- country[ctry_o]

  # get a good gap
  y_pos <- plotrix::spreadout(rate, gap)
  text(x = max_year + 1,
       y = y_pos,
       labels = country,
       col = col,
       cex = cex,
       adj = adj,
       xpd = xpd,
       ...)

}


splitGeoNames <- function (geo) {
  # split country/years (in format 'country_year') out of rownames
  # of a geostatistical conditional simulation object and add as columns
  splits <- strsplit(rownames(geo), split = '_')
  ctry <- sapply(splits, '[', 1)
  year <- sapply(splits, '[', 2)
  geo <- data.frame(iso3 = ctry,
                    year = year,
                    geo,
                    stringsAsFactors = FALSE)
  rownames(geo) <- NULL
  return (geo)
}

getEst <- function (df, iso3, year) {
  # for each dataframe of estimates, line up the country
  # and year, add object name as prefix to the column nmes and
  # return as a dataframe

  # get object name
  prefix <- deparse(substitute(df))

  # get index
  iso3_year_target <- paste(iso3, year, sep = '_')
  iso3_year_df <- paste(df$iso3, df$year, sep = '_')
  idx <- match (iso3_year_target, iso3_year_df)

  # subset, add prefix to column names and return
  df <- df[idx, -(1:2)]
  colnames(df) <- paste(prefix, colnames(df), sep = '_')
  return (df)
}

elogit <- function (y, n) log ( (y + 0.5) / (n - y + 0.5) )

## Combine input data and covariates layers for models run by region
## Necessary for saving to Shiny tool
## Save "covs", "tv_*", and "df" to new combined snapshot in model_image_history
combine_region_image_history <- function(indicator, indicator_group, run_date, fixed_effects) {

  load(paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/model_image_history/', run_date, '_bin0_wssa_0.RData'))
  new_cov_list <- list()
  for(cov in names(cov_list)[!grepl("gaul_code", names(cov_list))]) {
    pull_raster_covs <- function(region) {
      load(paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/model_image_history/', run_date, '_bin0_', region, '_0.RData'))
      cov_raster <- cov_list[[cov]]
      return(cov_raster)
    }
    region_layers <- lapply(Regions, pull_raster_covs)
    combined_layers <- do.call(raster::merge, region_layers)
    names(combined_layers) <- gsub("layer", cov, names(combined_layers))
    if(length(names(combined_layers))==1) new_cov_list[[cov]] <- combined_layers
    if(length(names(combined_layers))!=1) assign(paste0('tv_', cov), combined_layers)
  }
  cov_list <- brick(new_cov_list)

  # Combine input data
  pull_df <- function(region) {
    load(paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/model_image_history/', run_date, '_bin0_', region, '_0.RData'))
    return(df)
  }
  df <- lapply(Regions, pull_df)
  df <- do.call(rbind.fill, df)

  covs <- cov_list

  save(list = c('df','covs',grep('^tv_*', ls(), value = TRUE)), file = paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/model_image_history/', run_date, '.RData'))

}


      cirange = function(x){
          z=quantile(x,probs=c(.025,.975),na.rm=T)
          return(z[2]-z[1])
      }
      lower = function(x) quantile(x,probs=.025,na.rm=T)
      upper = function(x) quantile(x,probs=.975,na.rm=T)

cleanup_inla_scratch <- function(run_date) {

  if(keep_inla_files==FALSE) {

    # Clean up INLA intermediate directories unless user has specified to keep them.
    inla_working_dir <- '<<<< FILEPATH REDACTED >>>>/inla_intermediate'
    inla_working_dirs <- list.dirs(inla_working_dir, recursive = FALSE)
    inla_working_dirs <- inla_working_dirs[grepl(run_date, inla_working_dirs)]
    for(inla_dir in inla_working_dirs) {
      unlink(inla_dir, recursive=TRUE)
    }

  }

  if(keep_inla_files==TRUE) {

    message('Keeping INLA intermediate files because keep_inla_files==TRUE in config.')
    message(paste0('Files stored here: <<<< FILEPATH REDACTED >>>>inla_intermediate/inla_', run_date))

  }

}

plot_raster <- function(ras,
                        indicator,
                        region = 'africa',
                        return_map = T,
                        out_dir = NULL,
                        highisbad = T,
                        min_value = NULL,
                        mid_value = NULL,
                        max_value = NULL,
                        legend_title = NULL,
                        plot_title = NULL,
                        layer_names = NULL,
                        cores = 12,
                        individual_layers = T) {

  # Inputs & Outputs ########################################################

  # Inputs:
  #     ras           = RasterBrick object (in future may extend to single layers)
  #     indicator     = your indicator name
  #     region        = region to map ('africa' only option for now)
  #     return_map    = do you want the function to return ggplot objects?
  #
  #     out_dir       = where you want the .png files to go
  #                       if null, no files returned (can use for map obj only)
  #
  #     highisbad     = are high values bad (red) in your data?
  #
  #     min_value     = minimum value for range of color bar
  #                       if not provided, minimum value in data
  #     mid_value     = inflection point for color bar
  #                       default is mean of min_value, max_value
  #     max_value     = maximum value for range of color bar
  #                       if not provided, max value in data
  #
  #     legend_title  = title for legend
  #     plot_title    = title for plot
  #     layer_names   = names for each layer of RasterBrick
  #                       most common use = year
  #
  #
  # Outputs:
  #       png files for faceted view of all layers & individual layers
  #       ggplot objects (if return_map = T) in list format
  #
  #
  #
  ##########################################################################

  # runtime configuration for cluster
  if (Sys.info()["sysname"] == "Linux") {
    j_root <- "<<<<< FILEPATH REDACTED >>>>>"
  } else {
    stop("Run on cluster for access to shapefiles, etc.")
  }

  # Load packages ------------------------------------------------------------------

  package_lib <- paste0(<<<< FILEPATH REDACTED >>>>) # Library for all MBG versioned packages.
  .libPaths(package_lib)

  pacman::p_load(data.table, tidyr, magrittr, sp, rgdal, ggplot2, dplyr, RColorBrewer, parallel)
  library(raster, lib=paste0(<<<< FILEPATH REDACTED >>>>))

  # Must load this copy of extrafont (correct fonts installed)
  library(extrafont, lib=paste0(<<<< FILEPATH REDACTED >>>>))

  # Load packages
  message("\nLoading packages...")

  # Initial checks -----------------------------------------------------------------

  # Check to see if brick or single raster
  if ("RasterBrick" %in% class(ras)) {
    is_brick = T

    if(is.null(layer_names)) stop("Must enter layer names if plotting raster brick")
    if(nlayers(ras) != length(layer_names)) stop("Length of layer_names must match number of RasterBrick layers")

  } else {
    is_brick = F
    stop("For now, must use RasterBrick")
  }

  # Load background map - by region ------------------------------------------------
  # For now, only support 'africa'

  message("Loading background map...")

  if (region == 'africa') {
    background_shp <- readOGR(dsn=<<<< FILEPATH REDACTED >>>>, layer="africa_ad0")
  } else {
    # TODO: add other regions here in else block
    stop("Region not currently supported")
  }

  background_map <- fortify(background_shp)

  # Crop and mask data
  ras <- crop(ras, extent(background_shp))
  ras <- mask(ras, background_shp)

  map_points <- mclapply(as.list(ras),
                         function(x) rasterToPoints(x) %>% as.data.table,
                        mc.cores = cores)

  names(map_points) <- as.character(layer_names)
  map_df <- mclapply(names(map_points),
                     function(n){
                       x <- map_points[[n]]
                       x <- data.table(x)
                       names(x) <- c("long", "lat", "value")
                       x[, layer := n]
                   }, mc.cores = cores)

  map_df <- rbindlist(map_df)

  # Generate a range if not specified ----------------------------------------------
  if (is.null(min_value)) min_value <- min(minValue(ras)) # get global min for layers
  if (is.null(max_value)) max_value <- max(maxValue(ras)) # get global max for layers

  if (is.null(mid_value)) {
    mid_value <- rowMeans(cbind(min_value, max_value), na.rm=TRUE)
    add_mid_value <- F
  } else {
    add_mid_value <- T
  }

  # Set up color gradient ----------------------------------------------------------
  color_gradient <- c("#4575b4", "#FFFFBF", "#d73027")
  breaks <- c(min_value, mid_value, max_value)

  # Flip if high = good
  if (highisbad == F) color_gradient <- rev(color_gradient)

  # A function to generate breaks
  my_breaks <- function(x){
    breaks <- c(min(x),mean(x),max(x))

    # Add mid_value if was originally manually specified
    if (add_mid_value == T) {
      breaks <- c(breaks, mid_value)
    }

    names(breaks) <- attr(breaks,"labels")
    breaks
  }

  # Generate a ggplot2 component for the color gradient
  color_addin <- scale_fill_gradientn(colors = color_gradient,
                                      na.value = "#a6d854",
                                      values = breaks,
                                      breaks = my_breaks,
                                      limits = c(min_value, max_value))

  # Prep other map inputs ----------------------------------------------------------

  # Create the legend
  if(is.null(legend_title)) legend_title <- "Value"

  # Prepare layer names
  if(!is.null(layer_names)) {
    inset_df <- cbind(rep(-20, length(layer_names)),
                      rep(-30, length(layer_names)),
                      as.character(layer_names),
                      as.character(layer_names)) %>% as.data.table
    names(inset_df) <- c("x", "y", "layer", "text")
    inset_df$x <- as.numeric(inset_df$x)
    inset_df$y <- as.numeric(inset_df$y)
    inset_addin <- geom_text(data = inset_df,
                             aes(x=x,
                                 y=y,
                                 label=text),
                             size = 4,
                             hjust = 0.5)
  } else {
    inset_addin <- NULL
  }

  # Generate the main map ----------------------------------------------------------

  # A custom, mostly blank theme
  theme_empty <- theme_classic() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          text = element_text(family = "Open Sans"))

  # Initialize map
  g_map <- ggplot() +
    geom_polygon(data = background_map,
                 aes(x = long,
                     y = lat,
                     group = group),
                 fill = "lightgray") +
    geom_raster(data = map_df,
                aes(x=long,
                    y = lat,
                    fill = value)) +
    geom_path(data=background_map,
              aes(x=long, y=lat, group=group),
              color='black', lwd=.1) +
    coord_equal(ratio = 1) +
    labs(fill = paste0(legend_title, "\n"),
         title = plot_title) +
    theme_empty +
    color_addin +
    inset_addin +
    facet_wrap(~layer)

  # Outputs

  message(paste0("Saving main map for ", indicator))

  if(!is.null(out_dir)) {
    png(filename = paste0(out_dir, indicator, ".png"),
        type = "cairo",
        units = "in",
        width = 11,
        height = 9,
        pointsize = 12,
        res = 450)

    print(g_map)

    dev.off()
  }

  if(return_map) {
    map_list <- list()
    map_list[["main"]] <- g_map
  }

  # Generate individual maps -------------------------------------------------------
  split_layer_list <- lapply(layer_names, function(lyr) {
      map_df_single <- map_df[layer == lyr]
      inset_df_single <- inset_df[layer == lyr]
      return_list <- list(map_df_single, inset_df_single, lyr)
      names(return_list) <- c("map_df_single", "inset_df_single", "layer")
      return(return_list)
    })

  if (individual_layers == T) {
    layer_map_list <- mclapply(split_layer_list, mc.cores = cores, FUN = function(input_list) {
      lyr <- input_list[["layer"]]
      map_df_single <- input_list[["map_df_single"]]
      inset_df_single <- input_list[["inset_df_single"]]

      inset_addin_single <- geom_text(data = inset_df_single,
                               aes(x=x,
                                   y=y,
                                   label=text),
                               size = 10,
                               hjust = 0.5)

      # Initialize map
      g_map_single <- ggplot() +
        geom_polygon(data = background_map,
                     aes(x = long,
                         y = lat,
                         group = group),
                     fill = "lightgray") +
        geom_raster(data = map_df_single,
                    aes(x=long,
                        y = lat,
                        fill = value)) +
        geom_path(data=background_map,
                  aes(x=long, y=lat, group=group),
                  color='black', lwd=.1) +
        coord_equal(ratio = 1) +
        labs(fill = paste0(legend_title, "\n"),
             title = plot_title) +
        theme_empty +
        color_addin +
        inset_addin_single

      message(paste0("Saving map for ", indicator, " layer ", lyr))

      if(!is.null(out_dir)) {
        png(filename = paste0(out_dir, indicator, "_", lyr, ".png"),
            type = "cairo",
            units = "in",
            width = 10,
            height = 9,
            pointsize = 12,
            res = 600)

        print(g_map_single)

        dev.off()
      }

      if(return_map) {
        return(g_map)
      }

    })
  }

  if(return_map) {
    map_list <- c(map_list, layer_map_list)
    return(map_list)
  }

}

# Timer Functions ---------------------------------------------------------

## Functions to time within the child stacking regions
##
## General usage:
##    require(tictoc)
##
##    tic("Step 1")
##    **your code here**
##    toc(log = T)
##
##    ticlog <- tic.log(format = F)
##    generate_time_log(ticlog)
##
##  Returns: data table with two columns
##     "step": names of events (e.g. "Step 1")
##     "time": time elapsed (as text: Xh Xm Xs)
##
##  Note: can nest tic/toc pairs

generate_time_log <- function(ticlog) {

  # Set up packages
  require(magrittr)
  require(data.table)

  # Functions in functions
  strip_time <- function(x) {
    sec <- as.numeric(x$toc - x$tic)
    time <- format_time(sec)
    name <- x$msg

    df <- c(name, time) %>%
            t %>%
            as.data.table

    names(df) <- c("step", "time")

    return(df)
  }

  format_time <- function(run_time) {
    run_time <- round(as.numeric(run_time),0)

    hours <- run_time %/% 3600
    remainder <- run_time %% 3600
    minutes <- remainder %/% 60
    seconds <- remainder %% 60

    run_time <- paste0(hours, "h ", minutes, "m ", seconds, "s")
    return(run_time)
  }

  df_out <- lapply(ticlog, strip_time) %>% rbindlist

  return(df_out)

}

graph_run_summary <- function(run_date,
                              indicator,
                              indicator_group,
                              return_graph = F) {

  # Function to graph run_summary files
  # Requires creation of a run_summary .csv file in your run_date directory

  if (Sys.info()["sysname"] == "Linux") {
    j_root <- "<<<<< FILEPATH REDACTED >>>>>"
    package_lib <- paste0(j_root,'/temp/geospatial/packages') # Library for all MBG versioned packages.
    .libPaths(package_lib)
  }

  require(data.table)
  require(magrittr)
  require(ggplot2)
  require(RColorBrewer)

  dir <- paste0("<<<< FILEPATH REDACTED >>>>", indicator_group, "/", indicator, "/output/", run_date, "/")

  file <- list.files(dir, pattern = "run_summary.*csv")
  file <- paste0(dir, file)

  df <- read.csv(file, stringsAsFactors = F) %>% as.data.table

  grab_time_hours <- function(x) {

    v_time <- unlist(strsplit(x, " "))
    hours <- v_time[1]
    hours <- substr(hours, 0, nchar(hours) - 1) %>% as.numeric

    minutes <- v_time[2]
    minutes <- substr(minutes, 0, nchar(minutes) - 1) %>% as.numeric

    seconds <- v_time[1]
    seconds <- substr(seconds, 0, nchar(seconds) - 1) %>% as.numeric

    hours <- hours + minutes/60 + seconds/(60*60)
    return(hours)

  }

  df$time <- sapply(df$time, grab_time_hours)
  df$step <- factor(df$step, levels = c("Stacking - GAM", "Stacking - GBM", "Stacking - lasso", "Stacking - ridge", "Stacking - enet",
                                        "MBG - fit model", "MBG - predict model", "Cross-validation",
                                        "Stacking - all", "MBG - all", "Entire script"))

  summary_steps <- c("Stacking - all", "MBG - all", "Entire script")

  df_graph <- subset(df, !(step %in% summary_steps))
  df_graph[, holdout := as.character(holdout != 0)]

  g_plot <- ggplot(data = df_graph, aes(x = step, y = time, color = region)) +
    stat_summary(fun.y = mean, geom = "line", aes(group = region, color = region)) +
    geom_point(aes(color=region, shape=holdout)) +
    scale_color_brewer(palette = "Set1") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Step",
         y = "Time (hours)",
         title = paste0("Run Date: ", run_date),
         color = "Region",
         shape = "Holdout?")

  png(filename = paste0(dir, "run_summary_", indicator, "_", run_date, ".png"),
      type = "cairo",
      units = "in",
      width = 8,
      height = 4.5,
      pointsize = 12,
      res = 300)

  print(g_plot)

  dev.off()

  if (return_graph == T) return(g_plot)

}


clean_model_results_table <- function(rd   = run_date,
                                      regs = Regions,
                                      ages = 0,
                                      nm   = '',
                                      indic = indicator,
                                      ig = indicator_group,
                                      stackers = stacked_fixed_effects, 
                                      coefs.sum1 = as.logical(coefs_sum1)){

  str_match <- stringr::str_match

  # Adapted by Jon Mosser from Roy Burstein's `clean_model_results()` function
  require(magrittr)

  sharedir <- paste0('<<<< FILEPATH REDACTED >>>>', ig, '/', indic, '/output/', rd, '/')

  # make loopvars
  lv <- expand.grid(regs,ages)

  # grab formatted model fit objects
  stacker_names <- strsplit(stackers, " + ", fixed=T)[[1]]
  mods <- model_fit_table(lv=lv,rd=rd,nullmodel=nm, indicator = indic, indicator_group = ig,
                          coefs.sum1 = coefs.sum1, stacker_name_vec = stacker_names)

  # add region column
  all_mods <- lapply(names(mods), function(n) {

                mod <- mods[[n]] %>% as.data.table

                r <- str_match(n,"(.*)_")[1,2]
                a <- str_match(n, "_(.*)")[1,2]
                mod[, region := r]
                mod[, age := a]
                return(mod)
               })

  all_mods <- rbindlist(all_mods)
  colorder <- c("region", "age",  
                names(all_mods)[!(names(all_mods) %in% c("region", "age"))])
  setcolorder(all_mods, colorder)

  write.csv(all_mods, paste0(sharedir, indic, "_model_results_table.csv"),
            row.names = F)

  return(mods)

}

## ~~~~~~~~~~~ make table of model results ~~~~~~~~~~~~~~~~~
# Modified from RB's code by JM on 2017-08-01

model_fit_table <- function(lv=loopvars[loopvars[,3]==0,],
                            rd=run_date,
                            nullmodel='',
                            indicator = indicator,
                            indicator_group = indicator_group,
                            holdout = 0,
                            coefs.sum1,
                            stacker_name_vec){
  # load models
  require(INLA)
  message(sprintf('Pulling together table for %s models',rd))
  tlist=list()
  sharedir <- paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', rd, '/')
  for(i in 1:nrow(lv)){
    reg <- lv[i,1]
    age <- lv[i,2]
    message(sprintf('%i of %i loopvars. %s %i',i,nrow(lv),lv[i,1],lv[i,2]))

    # Recreate the INLA data stack
    pathaddin <- paste0('_bin', age, '_', reg, '_', holdout)
    f=paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/model_image_history/', rd, pathaddin, ".RData")

    if(!file.exists(f)){
      message('FAILED TO OPEN')
    } else {
      load(f)
    }

    modnames <- child_model_names
    df = df[,paste0(child_model_names) := lapply(child_model_names,
                                                 function(x) get(paste0(x,'_cv_pred')))]

    # Get the spde
    if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()
    input_data <- build_mbg_data_stack(df = df,
                                       fixed_effects = all_fixed_effects,
                                       mesh_s = mesh_s,
                                       mesh_t = mesh_t,
                                       use_ctry_res = FALSE,
                                       use_nugget = FALSE,
                                       exclude_cs    = c(modnames),
                                       usematernnew=F)

     spde <- input_data[[2]]

     # Now get & transform model fit
     message('::::loading in INLA fit\n')
     f=paste0(sharedir, indicator, '_model_eb_bin', age, "_", reg, "_0.RData")

    if(!file.exists(f)){
      message('FAILED TO OPEN')
    } else {
      load(f)
    }

     ## now we extract what we need from the fit to get transformed spatial params
    res.field <- inla.spde2.result(res_fit, 'space', spde, do.transf=TRUE)

    ## nominal range at 0.025, 0.5, 0.975 quantiles
    range   <- inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.range.nominal[[1]])
    nom.var <- inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.variance.nominal[[1]])
    spat.hyps <- rbind(range, nom.var)
    rownames(spat.hyps) <- c('Nominal Range', 'Nominal Variance')

    ## other hyperparmas
    hyps <- summary(res_fit)$hyperpar[-(1:2), ] ## first two rows are
                                                ## theta1, theta2 which
                                                ## we have in range and
                                                ## nom.var

    colnames(spat.hyps) <- colnames(hyps)[3:5]

    ## fixed effects
    if(coefs.sum1){
      fixed <- res_fit$summary.random$covar[,1:6]
      rownames(fixed) <- stacker_name_vec
    }else{
      fixed <- summary(res_fit)$fixed[,1:6]
    }
      
    ## combine them all and just keep three quantiles

    all.res <- rbind(fixed[, (3:5 + ifelse(coefs.sum1, 1, 0))], ## indices coming from summary.random are 4:6
                     spat.hyps,
                     hyps[, 3:5])

    ## rename GPRandom rho for time
    all.res <- as.data.table(all.res, keep.rownames = T)
    setnames(all.res, "rn", "parameter")
    all.res[parameter == "GroupRho for space", parameter := "GPRandom rho for time"]

    tlist[[paste0(reg, "_", age)]] <- all.res

  }

  return(tlist)
}



check_config <-function(){
  must_haves <- c(t(read.csv(sprintf('%s/config_must_haves.csv',commondir),header=FALSE)))
  message("\nRequired covariates: ")
  for(confs in must_haves)
    if(exists(confs)) message(paste0('  ',confs,': ',get(confs))) else stop(paste0(confs,' is missing, add it to your config'))
  extras <- config$V1[!(config$V1 %in% must_haves)]  
  if(length(extras) > 0) {
    message("\nAdditional covariates: ")
    for (extra in extras) message(paste0('  ',extra,': ',get(extra)))
  }
}

## proportional_area_map ################################################

#' Create proportional area maps for count data at various admin levels
#'
#' @param data data frame or data table with at least a value and ADMX_CODE column
#' @param ad_level admin level to map (0,1,2)
#' @param value_col column in `data` holding the value to plot (proportional to size)
#' @param main_title main title for the map (default = NULL)
#' @param legend_title legend title for the map 
#' @param fill_color color to use to fill in the proportional area circles
#' @param alpha transparency (alpha) of the proportional area circles
#' @param scale_size_max max size of the proportional area circles
#'                       (passed to `scale_size_area()`)
#' @param scale_size_breaks breaks to specify for size scale legend
#'                          (passed to `scale_size_area()`)
#' @param scale_size_breaks labels to specify for size scale legend
#'                          (passed to `scale_size_area()`)
#' @param lake_river_color controls the color of the lakes/rivers
#' @param out_file output file name (png)
#' @param out_file_height height of output file (inches)
#' @param out_file_width width of ouput file (inches)
#' @param out_file_res resolution of output file (dpi) 
#' @param out_file_pointsize pointsize of text in output file (generally would first
#'                           adjust legend text/title size as below with the relevant
#'                           options, and then use this to modify line spacing, etc.)
#' @param legend_text_size adjust the size of the legend text
#' @param legend_title_size adjust the size of the legend title
#' 
#' @return ggplot object containing the map object; writes png to `out_file` if not NULL
#' @examples
#' 
#' # Not run
#' a1_df <- fread("<<<< FILEPATH REDACTED >>>>[ig]/[indic]/output/[rd]/pred_derivatives/admin_summaries/[indic]_admin_1_raked_summary.csv")
#' a1_df <- subset(a1_df, year == 2015)
#' 
#' gg_obj <- proportional_area_map(data = a1_df,
#'                                 ad_level = 1,
#'                                 value_col = "mean",
#'                                 main_title = NULL,
#'                                 legend_title = "DPT3: 2015",
#'                                 fill_col = "red",
#'                                 alpha = 0.25,
#'                                 size_scale_max = 5,
#'                                 out_file = "/my/out/file.png")


proportional_area_map <- function(data, # Needs to have ADM0_CODE, ADM1_CODE, ADM2_CODE
                                  ad_level,
                                  value_col,
                                  main_title = NULL,
                                  legend_title,
                                  fill_col,
                                  alpha = 0.25,
                                  size_scale_max = 6,
                                  size_scale_breaks = NULL,
                                  size_scale_labels = NULL,
                                  lake_river_color = "lightblue",
                                  out_file = NULL,
                                  out_file_height=12,
                                  out_file_width=12,
                                  out_file_res=300,
                                  out_file_pointsize = 16,
                                  legend_text_size = 16,
                                  legend_title_size = 18) {

  # Set up and define objects
  library("ggplot2")
  library("data.table")
  library("magrittr")
  library("raster")
  library("rgeos")
  
  # Set up admin column code
  adm_code_col <- paste0("ADM", ad_level, "_CODE")
  
  # Copy input data to avoid environment/reference issues
  df <- copy(as.data.table(data))
  
  message("Loading background map...")
  message("Note: background map only available for Africa for now")
  background_shp <- readRDS('<<<< FILEPATH REDACTED >>>>/background_map_africa.rds')
  background_map <- fortify(background_shp)
  
  message("Loading master admin shape...")
  if (ad_level == 2) {
    ad_shape <- readRDS("<<<< FILEPATH REDACTED >>>>/g_2015_2014_2_modified.rds")  
    ad_shape <- subset(ad_shape, ADM2_CODE %in% unique(df[, get(adm_code_col)]))
  } else if (ad_level == 1) {
    ad_shape <- readRDS("<<<< FILEPATH REDACTED >>>>/g2015_2014_1.rds")
    ad_shape <- subset(ad_shape, ADM1_CODE %in% unique(df[, get(adm_code_col)]))
  } else if (ad_level == 0) {
    ad_shape <- readRDS("<<<< FILEPATH REDACTED >>>>/g2015_2014_0.rds")
    ad_shape <- subset(ad_shape, ADM0_CODE %in% unique(df[, get(adm_code_col)]))
  }
  
  message("Loading masks, lakes & rivers...")
  
  ### Load & process mask for ggplot2
  mask <- raster('<<<< FILEPATH REDACTED >>>>/mask_master.tif')
  mask_df <- rasterToPoints(mask)
  mask_df <- data.frame(mask_df)
  colnames(mask_df) <- c("long", 'lat', 'mask')
  
  ### Load & process lakes & rivers for ggplot2
  lakes <- raster('<<<< FILEPATH REDACTED >>>>/lakes_all_2.tif')
  lakes_df <- rasterToPoints(lakes)
  lakes_df <- data.frame(lakes_df)
  colnames(lakes_df) <- c("long", 'lat', 'lakes')
  
  message("Getting centroids...")
  if (ad_level == 2) centroids <- gCentroid(ad_shape, byid = T, id = ad_shape$ADM2_CODE)
  if (ad_level == 1) centroids <- gCentroid(ad_shape, byid = T, id = ad_shape$ADM1_CODE)
  if (ad_level == 0) centroids <- gCentroid(ad_shape, byid = T, id = ad_shape$ADM0_CODE)
  
  centroids <- as.data.frame(centroids) %>%
    cbind(rownames(.), .) %>%
    as.data.table(.) %>%
    setnames(., names(.), c("adm_code", "x", "y"))
  centroids$adm_code <- as.numeric(levels(centroids$adm_code))[centroids$adm_code]
  setnames(centroids, "adm_code", adm_code_col)
  
  message("Merging centroids to data")
  df <- merge(df, centroids, by = eval(adm_code_col), all.x = T, all.y = F)
  n_missing <- nrow(df[is.na(x) | is.na(y),])
  
  if (n_missing > 0) {
    missing_ad_shapes <- unique(df[is.na(x) | is.na(y), get(adm_code_col)]) 
    message(paste0("The following admin shapes are missing: ", missing_ad_shapes))
  }
  
  message("Preparing to plot...")
  ad_shape_df <- fortify(ad_shape, region = adm_code_col)
  
  # A custom, mostly blank theme
  theme_empty <- theme_classic() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          text = element_text(family = "Open Sans"))
  
  # Set up an object to hold the scale size options
  if (is.null(size_scale_breaks) & is.null(size_scale_labels)) {
    
    size_scale_object <- scale_size_area(max_size = size_scale_max)
    
  } else if (!is.null(size_scale_breaks) & is.null(size_scale_labels)) {
    
    size_scale_object <- scale_size_area(max_size = size_scale_max,
                                         breaks = size_scale_breaks)
    
  } else if (is.null(size_scale_breaks) & !is.null(size_scale_labels)) {
    
    size_scale_object <- scale_size_area(max_size = size_scale_max,
                                         labels = size_scale_labels)
    
  } else if (!is.null(size_scale_breaks) & !is.null(size_scale_labels)) {
    
    size_scale_object <- scale_size_area(max_size = size_scale_max,
                                         breaks = size_scale_breaks,
                                         labels = size_scale_labels)
  
  }
  
  # For convenience change df name of plotting column
  setnames(df, value_col, "plot_me")
  
  message("Plotting...")
  g_map <- ggplot() +
    
    # Plot background
    geom_polygon(data = background_map,
                 aes(x=long,
                     y=lat,
                     group = group),
                 fill = "white")  +
    
    # Plot national country outlines
    geom_path(data = background_map,
              aes(x = long,
                  y = lat,
                  group = group),
              size = 0.5,
              color = "black")   +
    
    # Plot admin unit outlines
    geom_path(data = ad_shape_df,
              aes(x = long,
                  y = lat,
                  group = group),
              size = 0.3,
              color = "black") +
    
    # Plot lakes & rivers
    annotate(geom = 'raster', x = lakes_df$long, y = lakes_df$lat,
             fill = lake_river_color,
             alpha = 0.75) +
    
    # Plot points
    geom_point(data = df,
               aes(x = x, 
                   y = y, 
                   size = plot_me),
               color = fill_col,
               alpha = alpha) +
    
    size_scale_object +
    
    # Theme/style additions & labels
    theme_empty +
    coord_equal(ratio = 1) +
    labs(size = legend_title,
         main = title) +
    guides(size = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.text=element_text(size=legend_text_size)) +
    theme(legend.title=element_text(size=legend_title_size))
  
  if (!is.null(out_file)) {
    message(paste0("Saving to ", out_file))
    png(filename = out_file,
        type = "cairo",
        units = "in",
        width = out_file_width,
        height = out_file_height,
        pointsize = out_file_pointsize,
        res = out_file_res)
    
    print(g_map)
    dev.off()
  }
  
  return(g_map)
}

