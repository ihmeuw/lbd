## #############################################################################
## 
## PLOT PARAMETER FITS AND MODEL OUTPUTS
## 
## Purpose: Plot model diagnostics: RFs, model hyperparameter and covariate fits
## 
## #############################################################################


rm(list=ls())
library(ggplot2)
library(ggrepel)


## DEFINE VARIABLES
## Run date
preset_run_date <-"<<<< REDACTED >>>>"
run_date <- preset_run_date
## Number of draws
samples <- 100
## Full set of countries modeled
full_region <- 'stage1+stage2-chn-bra-mex-mys-esh-guf'
## Split into modeling regions
all_regions <- c(
  'soas',
  'ocea+seas-mys',
  'stan+mng',
  'caca-mex', # Removing Mexico from Stage 2 paper
  'ansa+trsa-bra', # Removing Brazil from Stage 2 paper
  'noaf',
  'mide+yem',
  'cssa',
  'essa-yem',
  'sssa',
  'wssa'
)
config_name <- sprintf('config_died_%s',preset_run_date)
covs_config_name <- sprintf('covs_died_%s',preset_run_date)

## Set output directory
viz_dir <- paste0('<<<< FILEPATH REDACTED >>>>',gsub('-','_',Sys.Date()))
dir.create(viz_dir, showWarnings=FALSE)

## U5M Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

doannual <- TRUE # if false will do a 5-year run (used for testing)
testing  <- FALSE # run on a subset of the data
testgauss <- FALSE # Test SBH Data As Gaussian

## Set repo location and indicator group
user            <- Sys.info()['user']                                  # user running this code
core_repo       <- sprintf('<<<< FILEPATH REDACTED >>>>',user) # main MBG repo
ig_repo         <- sprintf('<<<< FILEPATH REDACTED >>>>',user)      # u5m specific repo
core_remote     <- 'origin'                                            # gits name of lbd_core upstream repo (should stay origin)
core_branch     <- 'dev-u5m'                                           # working branch to pull lbd_core from stash
ig_remote       <- 'origin'                                            # gits name of u5m upstream repo (should stay origin)
ig_branch       <- 'develop'                                           # woring branch to pull u5m from stash  
indicator_group <- 'u5m'                                               # indicator group is u5m, this shouldn't change
indicator       <- 'died'                                              # indicator is died, this shouldn't change 

makeholdouts    <- TRUE    # set to true if this is an out of sample run
pullgit         <- FALSE   # set to true if you will pull into the repo on share before running coe
numagebins      <- 7       # this depends on data prep, right now its 5, but it will likely change to 7
custom          <- FALSE   # specific for post estimation in India as MC set up

# run some setup scripts
source(sprintf('<<<< FILEPATH REDACTED >>>>',ig_repo)) 


## 1) Make fitted parameter table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir <- sprintf('<<<< FILEPATH REDACTED >>>>',run_date)
files <- list.files(dir,'fitted_parameter')

d <- data.table()

for (f in files) {
  tmp <- fread(sprintf('%s/%s',dir,f))
  r  <- gsub('.csv','',gsub('fitted_parameter_summary_table_','',f))
  tmp$region <- r
  d   <- rbind(d,tmp)
}

# make region names nice
d[ grepl("ansa|trsa", region), region := 'South America']
d[ grepl('caca'     , region), region := 'Central America and Caribbean']
d[ grepl('cssa'     , region), region := 'Central SSA']
d[ grepl('eaas'     , region), region := 'East Asia']
d[ grepl('mide'     , region), region := 'Middle East']
d[ grepl('essa'     , region), region := 'Eastern SSA']
d[ grepl('noaf'     , region), region := 'North Africa']
d[ grepl('ocea'     , region), region := 'Southeast Asia and Oceania']
d[ grepl('soas'     , region), region := 'South Asia']
d[ grepl('stan'     , region), region := 'Central Asia']
d[ grepl('sssa'     , region), region := 'Southern SSA']
d[ grepl('wssa'     , region), region := 'Western SSA']
d$param_name[d$param_name=='access2'             ]  <- 'Travel time to Cities'
d$param_name[d$param_name=='diarrhea_prev'       ]  <- 'Diarrhea Prevalence'
d$param_name[d$param_name=='dmspntl'             ]  <- 'Light Intensity'
d$param_name[d$param_name=='dpt3_cov'            ]  <- 'DPT3 Vaccination Coverage'
d$param_name[d$param_name=='hiv_test'            ]  <- 'HIV Prevalence'
d$param_name[d$param_name=='edu_mean_stage2_gadm']  <- 'Years of Maternal Eduation'
d$param_name[d$param_name=='fertility'           ]  <- 'Fertility Ratio'
d$param_name[d$param_name=='ghslurbanicity'      ]  <- 'Urbanicity Classification'
d$param_name[d$param_name=='ihmepm25'            ]  <- 'PM 2.5'
d$param_name[d$param_name=='map_pf_incidence'    ]  <- 'Malaria (Pf) Incidence'
d$param_name[d$param_name=='stunting_mod_b'      ]  <- 'Stunting Prevalence'
d$param_name[d$param_name=='worldpop'            ]  <- 'Log Population'
d$param_name[d$param_name=='year_cov'            ]  <- 'Year'

# remove hyperparams and intercept
dd <- d[!param_name%in%c('int','age_rho','country_RE_SD','NID_RE_SD','range','nominal_variance','tau','kappa','year_rho',paste0('FE_z_level__',1:7))]

# plot
png(sprintf('<<<< FILEPATH REDACTED >>>>',viz_dir),height=800,width=1600)
ggplot(dd,aes(x=param_name,y=median,color=region,group=region)) +
    geom_hline(yintercept=1,color='black') + ylab('Mean effect size and 95% UI') + 
    xlab('Fixed Effect Coefficient for Spatial Covariates') +
    geom_pointrange(aes(ymin=lower,ymax=upper),position=position_dodge(width=.5)) + theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=14))
dev.off()

## hyper parameters
dd <- d[param_name %in% c('tau','kappa','age_rho','year_rho','NID_RE_SD','country_RE_SD')]
dd$param_name[dd$param_name=='NID_RE_SD']      <- 'SD of Source Random Effect'
dd$param_name[dd$param_name=='country_RE_SD']  <- 'SD of Country Random Effect'
dd$param_name[dd$param_name=='age_rho']        <- 'AR1 Rho for Age'
dd$param_name[dd$param_name=='year_rho']       <- 'AR1 Rho for Year'
dd$param_name[dd$param_name=='tau']            <- 'Spatial Matern Tau'
dd$param_name[dd$param_name=='kappa']          <- 'Spatial Matern Kappa'

png(sprintf('<<<< FILEPATH REDACTED >>>>',viz_dir),height=1000,width=1100)
ggplot(dd,aes(x=1,y=median,color=region,group=region)) +
  ylab('Posterior fitted hyper-parameter with 95% UI') + 
  xlab('') + facet_wrap(~param_name,scales="free") +
  geom_pointrange(aes(ymin=lower,ymax=upper),position=position_dodge(width=.5)) + theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=14),
        strip.text.x = element_text(size = 14))
dev.off()

# Copy both files to the MBG directory
file.copy(
  from=paste0(
    viz_dir,
    c('<<<< FILEPATH REDACTED >>>>','<<<< FILEPATH REDACTED >>>>')
  ),
  to=dir,
  overwrite=TRUE
)


## 2) Make national and subnational raking plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(agebin in c('under5','neonatal','infant')){
  rd           = run_date
  keep_subnats = FALSE # change this to different options (ie agg subnats to national, make seperate subnat plot)
  keep_yl      = c(2000,2005,2010,2017) 
  facet_by     = 'year' # can be year or super_region, if NULL all as one plot
  color_by     = NULL #'model_region' # could be year or super_region, if NULL all will be black
  filename     = 'gbd_comparison'
  add_iso_lab  = TRUE
  
  # pull in RF databitties
  message('Combining regional files')
  dir    <- sprintf('<<<< FILEPATH REDACTED >>>>',agebin,rd) # always pull under5 because it has all the age bins in it
  files  <- list.files(path=dir,pattern ='_rf')
  d <- data.table()
  for(f in files){
    tmp              <- fread(sprintf('%s/%s',dir,f))
    modreg           <- gsub('_rf.csv','',gsub(sprintf('died_%s_',agebin),'',f))
    tmp$model_region <- modreg
    d <- rbind(d,tmp)
  }
  
  # merge in names
  lcm <- get_location_code_mapping(shapefile_version=modeling_shapefile_version)
  d   <- merge(d,lcm,by.x='loc',by.y='loc_id',all.x=TRUE)
  # subet year list
  d <- subset(d, year %in% keep_yl)
  dorig <- copy(d)

  # drop subnationals if that argument is set
  if(keep_subnats==FALSE){
    d <- subset(d, nchar(ihme_lc_id) == 3)
  }
  d <- subset(d, !is.na(ihme_lc_id))  
  # subset to year list
  # Drop countries not in stage 2
  d <- subset(
    d, 
    ADM_CODE %in% get_adm0_codes(
      full_region,shapefile_version=modeling_shapefile_version
    )
  )
  
  print(quantile(d$raking_factor))
  
  #  plot
  png(sprintf('%s/%s_%s.png',viz_dir,filename,agebin),
      res=300,height=7,width=7,unit='in') 
  g <-
    ggplot(d, aes(x=start_point*1000,y=target*1000)) + 
      xlab('Before GBD Calibration (Mortality per 1000 live births)') +
      ylab('After GBD Calibration (Mortality per 1000 live births)') + 
      theme_minimal()
  if(!is.null(facet_by)) g <- g + facet_wrap(as.formula(paste("~", facet_by)))
  if(is.null(color_by)){
    g <- g + geom_point(size=2,alpha=.5)
  } else {
    # TODO: clean up legend name if this option is set
    g <- g + geom_point(aes(color=get(color_by)))
  }
  if(add_iso_lab == TRUE){
    # point out problems
    g <- g + geom_text_repel( 
                             aes(label=loc_nm_sh),
                             segment.size  = 0.1,
                             segment.color = "grey10",
                             size=1.5)
  }
  g <- g + geom_abline(intercept=0,slope=1,color='red') 

    
  print(g)
  dev.off()
  
  # do some subnats only
  # fix for missing ken and eth
  d <- copy(dorig)
  d <- subset(d, nchar(ihme_lc_id) > 3)
  d[,iso:=substr(ihme_lc_id,1,3)]
  color_by <- 'country'
  add_iso_lab <- FALSE
  print(quantile(d$raking_factor))
  cntrys <- data.table(iso=c("CHN", "MEX", "IDN", "BRA", "IND"),
                       country=c('China','Mexico','Indonesia','Brazil','India'))
  d <- merge(d,cntrys,by='iso')
  # plot
  png(sprintf('%s/%s_%s_subnational.png',viz_dir,filename,agebin),
      res=300,height=7,width=7,unit='in')
  g <-
    ggplot(d, aes(x=start_point*1000,y=target*1000)) + 
    xlab('Before GBD Calibration (Mortality per 1000 live births)') +
    ylab('After GBD Calibration (Mortality per 1000 live births)') + 
    theme_minimal()
  if(!is.null(facet_by)) g <- g + facet_wrap(as.formula(paste("~", facet_by)))

  g <- g + geom_point(aes(color=country))

  if(add_iso_lab == TRUE){
    # point out problems
    g <- g + geom_text_repel(data=d[year==max(keep_yl) & !raking_factor%between%c(.7,1.3)],
                             aes(label=ihme_lc_id),
                             segment.size  = 0.1,
                             segment.color = "grey10",
                             size=1.5)
  }
  g <- g + geom_abline(intercept=0,slope=1,color='red')
  print(g)
  dev.off()
  
  # plot each country
  for(c in unique(d$country)){
    tmp <- d[country==c,]
    png(sprintf('%s/%s_%s_subnational_%s.png',viz_dir,filename,agebin,c),
        res=300,height=7,width=7,unit='in') 
    g <-
      ggplot(tmp, aes(x=start_point*1000,y=target*1000)) + 
      xlab('Before GBD Calibration (Mortality per 1000 live births)') +
      ylab('After GBD Calibration (Mortality per 1000 live births)') + 
      theme_minimal()
    if(!is.null(facet_by)) g <- g + facet_wrap(as.formula(paste("~", facet_by)))
    
    g <- g + geom_point(aes(color=country))

    g <- g + geom_abline(intercept=0,slope=1,color='red') 
    
    print(g)
    dev.off()
  }
  # histograms of raking factors across subnats
  png(sprintf('%s/%s_%s_subnational_histograms.png',viz_dir,filename,agebin),
      res=300,height=7,width=7,unit='in') 
  ggplot(d,aes(x=raking_factor)) + 
    geom_histogram() + theme_minimal() + facet_grid(country~year) +
    geom_vline(xintercept = 1)
  
  dev.off()  
}
