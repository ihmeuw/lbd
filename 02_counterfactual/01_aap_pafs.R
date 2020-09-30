###################################################################################
# Generate RR and PAFs from integrated exposure response (IER) curve
###################################################################################

lri_in_dir <- '<<< FILEPATH REDACTED >>>'

#(1) Setup

# clear memory
rm(list=ls())

arg <- commandArgs()[-(1:5)]  # First args are for unix use only

if (length(arg)==0) {
  #toggle for targeted run on cluster
  arg <- c("33", #location
           2017, #year
           34, #exp grid version
           "33power2_simsd_source_priors", #ier version
           43, #output version
           1000, #draws required
           58934, #stgpr HAP run_id
           5, # cores
           "021119")  #hap map date
}

# load packages, install if missing

lib.loc <- paste0(h_root,"R/",R.Version()$platform,"/",R.Version()$major,".",R.Version()$minor)
dir.create(lib.loc,recursive=T, showWarnings = F)
.libPaths(c(lib.loc,.libPaths()))

packages <- c("data.table","magrittr","fst","parallel")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

# set working directories
home.dir <- file.path('<<< FILEPATH REDACTED >>>')
setwd(home.dir)

# set project values
location_set_id <- 22

# Set parameters from input args
this.country <- arg[1]
this.year <- arg[2]
exp.grid.version <- arg[3]
ier.version <- arg[4]
output.version <- arg[5]
draws.required <- as.numeric(arg[6])
run_id <- arg[7]
cores <- arg[8]
hap.map.date <- arg[9]


draw.cols <- paste0("draw_", 1:draws.required)
silly.cols <- paste0("draw_", 0:(draws.required-1)) #dumb formatting
hap.map.cols <- paste0("hap_map_", 1:draws.required)
hap.prop.cols <- paste0("hap_prop_", 1:draws.required)
hap.ratio.cols <- paste0("hap_ratio_",1:draws.required)
ambient.cols <- paste0("ambient_",1:draws.required)
paf.cols <- paste0("draw_", 0:(draws.required-1)) #must be saved in 0-999 format

exp_cats <- c("ambient","hap")

ier_function <- function(x,alpha,beta,gamma,tmrel) {ifelse(x > tmrel, 1 + alpha*(1-exp(-beta*(x-tmrel)^gamma)), 1)}

outcomes <- rbindlist(list(data.table(cause="lri", age=99, est_incidence=0, cause_id="322", group="children")))


# Directories -------------------------------------------------------------

#AiR PM functions#
air.function.dir <- file.path(h_root, 'code/air/_lib')

# this pulls the miscellaneous helper functions for air pollution
file.path(air.function.dir, "misc.R") %>% source

#general functions#
central.function.dir <- file.path(h_root, "code/tools/")
# this pulls the general misc helper functions
file.path(central.function.dir, "misc.R") %>% source



hap.prop.in <- '<<< FILEPATH REDACTED >>>'
hap.map.in <- '<<< FILEPATH REDACTED >>>'
ambient.grid.in <- '<<< FILEPATH REDACTED >>>'
ier.dir <- '<<< FILEPATH REDACTED >>>'

pm.sum.out <- '<<< FILEPATH REDACTED >>>'
rr.sum.out <- '<<< FILEPATH REDACTED >>>'
paf.sum.out <- '<<< FILEPATH REDACTED >>>'

dir.create(pm.sum.out, showWarnings = T)
dir.create(rr.sum.out, showWarnings = T)
dir.create(paf.sum.out, showWarnings = T)

reis <- c("air_pmhap","air_pm","air_hap")

for(rei in reis){
  
  dir.create('<<< FILEPATH REDACTED >>>', recursive = T, showWarnings = F)
  dir.create('<<< FILEPATH REDACTED >>>', recursive = T, showWarnings = F)
}

# Read in data ----------------------------------------------------------------

hap.prop <- fread(paste0(hap.prop.in,"/",this.country,".csv"))[year_id==this.year]
setnames(hap.prop,silly.cols,hap.prop.cols)
hap.prop <- hap.prop[,c(hap.prop.cols),with=F]

hap.map <- fread(paste0(hap.map.in,"/lm_pred_",hap.map.date,".csv"))[location_id==this.country & year_id==this.year]
setnames(hap.map,draw.cols,hap.map.cols)
hap.map <- hap.map[,c(hap.map.cols),with=F]

hap.ratio <- fread(paste0(hap.map.in,"/crosswalk_",hap.map.date,".csv"))
setnames(hap.ratio,draw.cols,hap.ratio.cols)
hap.ratio <- hap.ratio[,c("grouping",hap.ratio.cols),with=F]

# read in ambient grid
ambient <- read.fst(paste0(ambient.grid.in,"/",this.country,"_",this.year,".fst")) %>% as.data.table
setnames(ambient,draw.cols,ambient.cols)
ambient <- ambient[,c("latitude","longitude","weight","pop",ambient.cols),with=F]

# make data long by draw for easier manipulation
hap.prop <- melt(hap.prop,measure.vars=hap.prop.cols,value.name="hap_prop",variable.name="draw",value.factor=T)
hap.prop[,draw:=as.numeric(draw)]

hap.map <- melt(hap.map,measure.vars=hap.map.cols,value.name="hap_map",variable.name="draw",value.factor=T)
hap.map[,draw:=as.numeric(draw)]

hap.ratio <- melt(hap.ratio,id.vars="grouping",measure.vars=hap.ratio.cols,value.name="hap_ratio",variable.name="draw",value.factor=T)
hap.ratio[,draw:=as.numeric(draw)]

ambient <- melt(ambient,measure.vars=ambient.cols,value.name="ambient",variable.name="draw",value.factor=T)
ambient[,draw:=as.numeric(draw)]

# merge hap datasets together
exp <- merge(hap.prop,hap.map,by="draw")
exp <- merge(exp,hap.ratio,by="draw",allow.cartesian=T)

# calculate total exposure for those exposed to hap
exp[,hap_excess:=hap_map*hap_ratio] #hap only exposure for each group (men,women,children) (above ambient levels)

# merge on grid cell ambient
exp <- merge(exp,ambient,by="draw",allow.cartesian=T)

# calculate total exposures
exp[,hap:=hap_excess+ambient]     #total exposure for those exposed to hap, used to get value from IER

# function to calculate for each IER
calculate <- function(i){
  row <- outcomes[i]
  
  params <- fread(paste0(ier.dir,"/params_",row$cause,"_",row$age,".csv"))
  params[,draw:=1:.N]
  
  exp <- merge(exp,params,by="draw")
  
  #subset to relevant groups (LRI is the only outcome evaluated in kids)
  if(row$group=="adults"){
    exp <- exp[grouping %in% c("female","male")]
  }else if(row$group=="all"){
    exp <- exp[grouping %in% c("female","male","child")]
  }else if(row$group=="children"){
    exp <- exp[grouping %in% c("child")]
  }
  
  
  for(exp_cat in exp_cats){
    
    #calculate RR for ambient and HAP
    exp[,paste0("rr_",exp_cat):=ier_function(x=get(exp_cat),alpha=alpha,beta=beta,gamma=gamma,tmrel=tmrel)]
    
  }
  
  # incidence and mortality, for CVD outcomes scale using predefined ratios. For others, no difference, duplicate draws
  if(row$est_incidence==1){
    
    if(row$cause=="cvd_ihd"){
      ratio <- 0.141
    }else if(row$cause=="cvd_stroke"){
      ratio <- 0.553
    }
    
    yll <- copy(exp)
    yll[,type:="yll"]
    yld <- copy(exp)
    yld[,type:="yld"]
    yld[,rr_hap:= rr_hap * ratio - ratio + 1]
    yld[,rr_ambient:= rr_ambient * ratio - ratio + 1]
    
    exp <- rbind(yll,yld)
    
  }else if(row$est_incidence==0){
    
    yll <- copy(exp)
    yll[,type:="yll"]
    yld <- copy(exp)
    yld[,type:="yld"]
    
    exp <- rbind(yll,yld)
  }
  
  
  #calculate summary PM rr as weighted mean
  exp[,"rr_pm":=rr_hap*hap_prop+rr_ambient*(1-hap_prop)]
  
  #Uses equation given by Steve to recalculate hap RR. Excess risk for the combined exposure minus the excess risk due to ambient plus 1
  #This will be used in RR_max calculation for SEVs
  exp[,"rr_hap":= (rr_hap-1) - (rr_ambient-1) + 1]
  
  # population weight shift, ambient exposure, RR
  exp <- merge(exp[,lapply(.SD,weighted.mean,w=pop*weight),.SDcols=c("rr_hap","rr_ambient","rr_pm","ambient"), by=c("grouping","draw","type")], #pop_weighted PM rr to calculate PM pafs, ambient and hap rrs for SEVs, exposure for proportional split
               exp[,c("hap_prop","hap_excess","grouping","draw"),with=F] %>% unique, #HAP doesn't differ by gridcell, only by group
               by=c("grouping","draw"))
  
  exp[,pop_average_pm:=ambient+hap_prop*hap_excess]  #population average exposure (denominator of proportion for splitting pafs)
  exp[,hap_paf_ratio:=(hap_prop)*(hap_excess)/pop_average_pm] #proportion of paf attributable to hap
  exp[,ambient_paf_ratio:=ambient/pop_average_pm] #proportion of paf attributable to ambient
  
  exp[,paf_pm:=(rr_pm-1)/rr_pm] #calculate PM paf based on weighted PM RR
  exp[,paf_hap:=paf_pm*hap_paf_ratio] #proprotionally split PM paf to Hap and ambient
  exp[,paf_ambient:=paf_pm*ambient_paf_ratio]
  
  exp[,location_id:=this.country]
  exp[,year_id:=this.year]
  exp[,cause:=row$cause]
  exp[,age:=row$age]
  
  return(exp)
}

out <- lapply(1:nrow(outcomes),calculate) %>% rbindlist


#summary file of exposures

lower <- function(x){quantile(x,p=.025)}
upper <- function(x){quantile(x,p=.975)}

summary <- melt(out,id.vars=c("grouping","location_id","year_id","draw","type","cause","age"))
summary <- summary[,.(mean=mean(value),lower=lower(value),upper=upper(value)),by=c("grouping","cause","age","location_id","year_id","type","variable")]


write.csv(summary[variable %in% c("ambient","hap_prop","hap_excess","pop_average_pm","hap_paf_ratio","ambient_paf_ratio"),.(grouping,variable,location_id,year_id,mean,lower,upper)] %>% unique,
          paste0(pm.sum.out,"/",this.country,"_",this.year,".csv"),row.names=F)

# summary file of RR

write.csv(summary[variable %in% c("rr_hap","rr_ambient","rr_pm")],
          paste0(rr.sum.out,"/",this.country,"_",this.year,".csv"),row.names=F)

# summary file of PAF

write.csv(summary[variable %in% c("paf_pm","paf_hap","paf_ambient")],
          paste0(paf.sum.out,"/",this.country,"_",this.year,".csv"),row.names=F)


# Save Draws --------------------------------------------------------------

out <- out[,.(location_id,year_id,cause,age,grouping,draw,type,rr_hap,rr_pm,rr_ambient,paf_pm,paf_hap,paf_ambient)]  

# reshape long by risk

out <- melt.data.table(out, id.vars=c("location_id","year_id","cause","age","grouping","draw","type"))
out[,c("measure","risk"):=tstrsplit(variable,"_")]
out[,variable:=NULL]

# reshape wide by draw
out[,draw:=paste0("draw_",draw)]
out <- dcast.data.table(out, location_id + year_id + cause + age + grouping + type + measure + risk ~ draw, value.var="value")

setnames(out,draw.cols,paf.cols)  
out[, acause := cause]


# expand cvd_stroke to include relevant subcauses in order to prep for merge to YLDs, using your custom find/replace function
# first supply the values you want to find/replace as vectors
old.causes <- c('cvd_stroke')
replacement.causes <- c('cvd_stroke_isch',
                        "cvd_stroke_intracerebral",
                        "cvd_stroke_subarachnoid")

# then pass to your custom function
out <- findAndReplace(out,
                      old.causes,
                      replacement.causes,
                      "acause",
                      "acause",
                      TRUE) #set this option to be true so that rows can be duplicated in the table join (expanding the rows)

# now replace each cause with cause ID
out[, cause_id := acause] #create the variable
# first supply the values you want to find/replace as vectors
cause.codes <- c('cvd_ihd',
                 "cvd_stroke_isch",
                 "cvd_stroke_intracerebral",
                 "cvd_stroke_subarachnoid",
                 "lri",
                 'neo_lung',
                 'resp_copd',
                 't2_dm')

cause.ids <- c(493,
               495,
               496,
               497,
               322,
               426,
               509,
               976)

# then pass to your custom function
out <- findAndReplace(out,
                      cause.codes,
                      cause.ids,
                      "cause_id",
                      "cause_id")

out <- out[, c("risk",
               'measure',
               "type",
               "grouping",
               "location_id",
               "year_id",
               "cause",
               "acause",
               "cause_id",
               "age",
               paf.cols),
           with=F]

out[risk=="ambient",risk:="air_pm"]
out[risk=="hap",risk:="air_hap"]
out[risk=="pm",risk:="air_pmhap"]

#RR save a file for each risk (air/air_pm/air_hap), cause 
saveRiskCause <- function(this.risk,this.cause,file){
  out <- copy(file)
  
  out <- out[cause_id==this.cause 
             & risk==this.risk
             & measure=="rr"] #subset to the relevant rows
  out[,measure:=NULL]
  
  write.csv(out, paste0(home.dir,"/",this.risk,"/rr/draws/",output.version,"/",
                        this.cause, "_", this.country, "_",
                        this.year, ".csv"),
            row.names=F)
}

save <- unique(out[,.(risk, cause_id)])
save[,saveRiskCause(risk,cause_id,file=out),by=1:nrow(save)]

#expand to the appropriate age groups using a custom function
out <- lapply(unique(outcomes$cause), expandAges, input.table=out) %>% rbindlist(use.names=T)

#create necessary variables
out[type=="yld", measure_id := 3]
out[type=="yll", measure_id := 4]

out <- out[, c("risk",
               'measure',
               "measure_id",
               "age_group_id",
               "sex_id",
               "location_id",
               "year_id",
               "acause",
               "cause_id",
               paf.cols),
           with=F]



saveSex <- function(this.sex, this.cause, this.measure, this.risk, file){
  
  out <- copy(file)
  
  out <- out[sex_id==this.sex 
             & cause_id==this.cause 
             & measure_id==this.measure
             & risk==this.risk 
             & measure=="paf"] #subset to the relevant rows
  out <- out[,risk:=NULL]
  out <- out[,measure:=NULL]
  
  write.csv(out, paste0(home.dir,"/",this.risk,"/paf/draws/",output.version, "/", 
                        this.measure, "_", this.cause, "_", this.country, "_",
                        this.year, "_", this.sex, ".csv"),
            row.names=F)
  
}

#PAF save save a file for each sex, cause, yll/yld, air/air_pm/air_hap
save <- unique(out[,.(sex_id,cause_id,measure_id,risk)])
save[,saveSex(sex_id,cause_id,measure_id,risk,file=out),by=1:nrow(save)]
