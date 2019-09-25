#collect the arguments
print(commandArgs())

#governing variables
indicator <- as.character(commandArgs()[4])
indicator_group <- as.character(commandArgs()[5])
run_date <- as.character(commandArgs()[6])
mainrd = substr(run_date, 1, 19)
agg_slots = 5
if(indicator=='has_lri'){
  indy = 'LRI'
}else{
  indy = indicator
}

#bits to run
shiny_diagnostics = as.logical(commandArgs()[7])
make_aggregates   = as.logical(commandArgs()[8])
make_rasters      = as.logical(commandArgs()[9])
oos_statistics    = as.logical(commandArgs()[10])
coverage_report   = as.logical(commandArgs()[11])
line_plots        = as.logical(commandArgs()[12])
combine_aggs      = as.logical(commandArgs()[13])
little_csvs       = as.logical(commandArgs()[14])
arocproj_results  = as.logical(commandArgs()[15])
arocproj_graphics = as.logical(commandArgs()[16])
oos_summaries     = as.logical(commandArgs()[17])

settings = c(indicator, indicator_group, run_date, shiny_diagnostics, make_aggregates,
             make_rasters, oos_statistics, coverage_report, line_plots, combine_aggs, little_csvs,
             arocproj_results, arocproj_graphics, oos_summaries, "")
for(setty in settings){
  message(setty)
}


#Source the required files
root = '<<<< FILEPATH REDACTED >>>>'

library('rgdal')
library('seegSDM')
library('rgeos')
library('magrittr')
library('data.table')
library('sp')
library('raster')
library('ggplot2')

repo = '<<<< FILEPATH REDACTED >>>>'


#load all the damn functions
funks = list.files('<<<< FILEPATH REDACTED >>>>', pattern = 'functions',full.names = T)
for(fff in funks){
  source(fff)
}

source('<<<< FILEPATH REDACTED >>>>') #post proc functions
source('<<<< FILEPATH REDACTED >>>>') #validation functions

#set directories
shelldir = '<<<< FILEPATH REDACTED >>>>'
datadir = '<<<< FILEPATH REDACTED >>>>'
diagfolder = '<<<< FILEPATH REDACTED >>>>'
dir.create(diagfolder)

pop_measure     <- 'a0004t'
age             <- 0
holdouts         <- 0
overwrite <- F
user = Sys.info()[["user"]]
repo            <- '<<<< FILEPATH REDACTED >>>>'
log_dir = '<<<< FILEPATH REDACTED >>>>'
measure_list = c('prevalence','incidence','mortality')

temprun = '<<<< FILEPATH REDACTED >>>>'
pest_objs = '<<<< FILEPATH REDACTED >>>>'
region_list = fetch_from_rdata(temprun,item_name = 'region_list')
year_list = eval(parse(text = fetch_from_rdata(temprun,item_name = 'year_list')))
summstats = fetch_from_rdata(pest_objs, item_name = 'summstats')
config = fetch_from_rdata(temprun, 'config')
samples = as.numeric(config$V2[config$V1=='samples'])
if(grepl('1000_draws', run_date)) samples = 1000

#create link to allow make_model_diagnostics to work
if(shiny_diagnostics){
  message('shiny_diagnostics')
  for(rrr in region_list){
    input = '<<<< FILEPATH REDACTED >>>>'
    stopifnot(file.exists(input))
    output = '<<<< FILEPATH REDACTED >>>>'
    file.copy(input, output)
  }


  ##launch things to send to shiny
  source('<<<< FILEPATH REDACTED >>>>') #diagnostic functions
  make_model_diagnostics(ig = indicator_group,
                         corerepo = '<<<< FILEPATH REDACTED >>>>',
                         indic = indicator,
                         rd = run_date,
                         script_dir = '<<<< FILEPATH REDACTED >>>>')
}

#make aggregates
if(make_aggregates){
  message('make_aggregates')
  agg_jobs = list()
  iter = 1
  for(rrr in c(T,F)){
    for(mmm in measure_list){
      for(reg in region_list){
        iter = iter +1

        agg_jobs[[iter]] = system(build_qsub(job_name = paste0('agg_',mmm,'_',reg,'_',rrr),
                                             output_folder = sprintf('%s/output/', datadir),
                                             error_folder = sprintf('%s/errors', datadir),
                                             make_folders = T,
                                             shell_path = paste0(shelldir, 'shell_sing.sh'),
                                             script_path = '<<<< FILEPATH REDACTED >>>>', #aggregate_results
                                             additional_options = paste0('-v SET_MKL_THREADS=',agg_slots, ' -v sing_image=default'),
                                             project = 'proj_geo_nodes',
                                             geos_node = T,
                                             slots = agg_slots,
                                             num_tasks = 0,
                                             arguments_string = paste(indicator, indicator_group, run_date,
                                                                      rrr, 'a0004t',T, 0, 0, reg, mmm, 'blarg')),intern=T)


      }
    }
  }
  agg_jobs = unlist(agg_jobs)
}

#while those jobs are running, get started on the rasters
#make rasters
if(make_rasters){
  message('make_rasters')
  for(rak in c('raked_', '')){

    if(rak == ''){
      mlist = 'prevalence'
    }else{
      mlist = measure_list
    }

    for(measure in mlist){
      print(measure)
      ras_files = lapply(summstats, function(sstat) '<<<< FILEPATH REDACTED >>>>')
      rastys = parallel::mclapply(ras_files, function(x) build_africa_raster(x, T, rak=='raked_'), mc.cores = 3)

    }
  }
}

#make oos stuff
#check to see if stratum_ho exists. if so, use that to determine the number of folds
if(oos_statistics | oos_summaries){
  message('oos')
  strat_ho_exists = function(){
    load('<<<< FILEPATH REDACTED >>>>')
    return(exists('stratum_ho'))
  }

  if(strat_ho_exists()){
    nfolds = 0:5
  }else{
    nfolds = 0
  }

  #run or load cached statistics
  if(oos_statistics){
    oos_res <- get_is_oos_draws(ind_gp = indicator_group,
                                ind = indicator,
                                rd = run_date,
                                ind_fm = 'binomial',
                                nperiod = length(year_list),
                                yrs = year_list,
                                write.to.file = T,
                                get.oos = strat_ho_exists(),
                                year_col = 'year')
  }else{
    oos_res = readRDS('<<<< FILEPATH REDACTED >>>>')
  }

  if(oos_summaries){

    #make a ho id in case it doesn't exist
    if(!any(names(oos_res) %in% 'ho_id')){
      oos_res[, ho_id := fold]
      warning('ho_id set to fold')
    }

    #summarize the results (via pv tables and what not if we ran oos)
    source('<<<< FILEPATH REDACTED >>>>')#validation functions
    if(length(nfolds)>1) summarize_oos_results(oos_res)
  }


}

#combine the aggregations; hold if they are not done
if(make_aggregates | combine_aggs){
  #check to see if anything is still running
  if(make_aggregates){
    jobs = get_running_jobs()
    proj_name <- '<<<< PROJECT NAME REDACTED >>>>'
    #if so, wait on those
    if(any(agg_jobs %in% jobs$jobid)){
      still_to_run = agg_jobs[agg_jobs %in% jobs$jobid]
      system(paste('qsub -b y -sync y -P ', proj_name,' -hold_jid', paste(still_to_run, collapse =','), paste('-N ', indicator, '_aggs'), shQuote('echo 1')))
    }

    Sys.sleep(5)
  }

  message('combine aggregates')
  #combine the aggregations
  for(rrr in c(T,F)){
    for(measure in c('prevalence','incidence','mortality')){ #measure_list 'daly','yll','yld'
      combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                          ages = 0,
                          regions = region_list,
                          holdouts = holdouts,
                          raked = rrr,
                          delete_region_files = F,
                          measure = measure)

      #save aggregates as csvs for mapping
      if(rrr){
        for(ad in 0:2){
          adm = summarize_admin(measure, ad)
          for(vvv in c('mean', 'counts')){
            param = ifelse(vvv == 'mean', 'mean_rate', 'mean_count' )
            adm[, value := get(param)]
            outj = '<<<< FILEPATH REDACTED >>>>'
            write.csv(adm[,c(paste0('ADM',ad,'_CODE'), 'year', 'value'), with = F], file ='<<<< FILEPATH REDACTED >>>>')
          }
        } #close ad loop
      } #close logical for saving admin csvs
    } #close measure loop
  }  #close rake loop
} #close run logical


##knit the coverage report
if(coverage_report){
  message('coverage_report')
  rmarkdown::render(input = '<<<< FILEPATH REDACTED >>>>', output_format = "pdf_document",
                    output_file = '<<<< FILEPATH REDACTED >>>>',
                    params = list(run_date = run_date), intermediates_dir = diagfolder, knit_root_dir = diagfolder, output_dir = diagfolder, envir = new.env())
}

##make model line plots
if(line_plots){
  message('line_plots')
  source('<<<< FILEPATH REDACTED >>>>') #plot admin time series
}

if(little_csvs){
  message('summary csvs')
  ad1 = rbindlist(lapply(c('incidence','prevalence','mortality'), function(x) summarize_admin(x, 1, draws = T)),fill = T)
  ad0 = rbindlist(lapply(c('incidence','prevalence','mortality'), function(x) summarize_admin(x, 0)),fill = T)
  ad2 = rbindlist(lapply(c('incidence','prevalence','mortality'), function(x) summarize_admin(x, 2)),fill = T)

  save(ad1, ad0, ad2, file = '<<<< FILEPATH REDACTED >>>>')


  iter = 0
  for(aaa in list(ad0,ad1,ad2)){
    saveRDS(aaa, paste0(diagfolder,'/', indicator,'_',run_date,'_draws','_ad_',iter,'.rds'))
    write.csv(aaa[,names(aaa)[!grepl('V[0-9]',names(aaa))], with = F], file = paste0(diagfolder,'/', indicator,'_',run_date,'_summary_tables','_ad_',iter,'.csv'),
              row.names = F)
    iter = iter +1
  }

}

if(arocproj_graphics | arocproj_results){
  message('arocproj')
  source('<<<< FILEPATH REDACTED >>>>') #aroc proj prob
}

message('fin')
