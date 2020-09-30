#################################################
# Purpose: Aggregate admin-level results for an
# entire group of indicators in order to compare them
# to one another.
################################################

# Setting up the workspace
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
root <- "<<<< FILEPATH REDACTED >>>>"
n_draws<-1000 #

indicator_group<-commandArgs()[4]; message(indicator_group)
results_pull_time<-commandArgs()[5]; message(results_pull_time)

# Setting the root, loading in libraries, and functions:
setwd("<<<< FILEPATH REDACTED >>>>")
for(function_script in list.files(getwd(),pattern="*_functions.R")){message(function_script);source(function_script)};message("Central Functions Loaded.")
load_libs(c('data.table','ggplot2'))

# Defining File Paths for Results Pulling and Output
results<-fread("<<<< FILEPATH REDACTED >>>>")
results_tables_path<-"<<<< FILEPATH REDACTED >>>>"
diffs_tables_path<-"<<<< FILEPATH REDACTED >>>>"
results_time_series_path<-"<<<< FILEPATH REDACTED >>>>"

calc_diffs<-T

# Loading in the results at each admin level
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Generating empty results lists
regs<-list()
adm0<-list()
adm1<-list()
adm2<-list()
sp_h<-list() # spatial hierarchy

# Propogating results from each of the indicators into the lists
for (i in 1:nrow(results)){
  print(results$indicator_longname[i])
  load("<<<< FILEPATH REDACTED >>>>")
  regions[,indicator:=results$indicator_longname[i]]
  admin_0[,indicator:=results$indicator_longname[i]]
  admin_1[,indicator:=results$indicator_longname[i]]
  admin_2[,indicator:=results$indicator_longname[i]]
  sp_hierarchy_list[,indicator:=results$indicator_longname[i]]
  regs[[ results$indicator_longname[i]]]<-regions
  adm0[[ results$indicator_longname[i] ]]<-admin_0
  adm1[[ results$indicator_longname[i] ]]<-admin_1
  adm2[[ results$indicator_longname[i] ]]<-admin_2
  sp_h[[ results$indicator_longname[i] ]]<-sp_hierarchy_list
}

# Re-combining the results into large data.frames for the entire area
regs<-rbindlist(regs)
adm0<-rbindlist(adm0)
adm1<-rbindlist(adm1)
adm2<-rbindlist(adm2)
sp_h<-rbindlist(sp_h)
rm(sp_hierarchy_list,admin_0,admin_1,admin_2,regions,sp_h)

# Ignoring the sp_h file generated, loading in the spatial names, etc, from the level 2 DBF
admin_level<-2
sp_h<-data.table(foreign::read.dbf("<<<< FILEPATH REDACTED >>>>"))
sp_h<-sp_h[ADM0_CODE %in% adm0$ADM0_CODE,]

# Generating summary tables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
draw_colnames<-paste0("V",seq(1:n_draws))

# Defining functions to get upper and lower confidence intervals:
lower_confint<-function(x){z=quantile(x,probs=c(.025),na.rm=T);return(z)}
upper_confint<-function(x){z=quantile(x,probs=c(.975),na.rm=T);return(z)}

# Making difference files across time (in this case, 2000 and 2015):
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# I'm sure this could be refactored into a more elegant solution, but this should work for now:
  if(calc_diffs){
    # Making sure things are ordered correctly such that overlapping the confidence intervals works
    adm0<-adm0[order(year,indicator,ADM0_CODE)]
    adm1<-adm1[order(year,indicator,ADM1_CODE)]
    adm2<-adm2[order(year,indicator,ADM2_CODE)]

    # Getting subsets of the data.tables that will serve as the IDs
    adm0_ids<-adm0[year==2000,][,list(indicator,ADM0_CODE)]
    adm1_ids<-adm1[year==2000,][,list(indicator,ADM1_CODE)]
    adm2_ids<-adm2[year==2000,][,list(indicator,ADM2_CODE)]

    adm0_diffs_list<-list()
    adm1_diffs_list<-list()
    adm2_diffs_list<-list()

    # Getting all permutations of the landmark years (2000,2005,2010,2015, etc...)
      year_combos<-CJ(year_1=seq(2000,2015,5),year_2=seq(2000,2015,5))
      year_combos<-year_combos[year_1<year_2,] # Only keep the combinations in which the second year is later than the first year

    for (y in 1:nrow(year_combos)){ # For each of the  "landmark" years (2000,2005,2010,2015..)
      year_1<-year_combos[y]$year_1
      year_2<-year_combos[y]$year_2
      message(paste0("Making diffs for ",year_1," to ",year_2))

      # Admin 0 differences
      message("   admin 0")
        est_yr1<-adm0[year==year_1,][,draw_colnames,with=F] # Get just the draw columns, ordered the same way
        est_yr2<-adm0[year==year_2,][,draw_colnames,with=F]
        adm0_diffs<-100*(est_yr2-est_yr1)/est_yr1 # basically overlapping the data.tables as matrices to calculate across rows and columns simultaneously
        adm0_diffs[,year:=paste0("pct_",year_1,"_",year_2)]
        adm0_diffs<-cbind(adm0_ids,adm0_diffs) # Adding back on the ID columns

        adm0_diffs<-adm0_diffs[,list(lower=apply(.SD,FUN="lower_confint",MARGIN=1),
                                     mean=apply(.SD,FUN="mean",MARGIN=1),
                                     upper=apply(.SD,FUN="upper_confint",MARGIN=1),
                                     ADM0_CODE,year,indicator),.SDcols=draw_colnames]
        adm0_diffs<-merge(adm0_diffs,unique(sp_h[,list(ADM0_CODE,ADM0_NAME)]),by="ADM0_CODE",all.x=T)
        adm0_diffs_list[[paste0("pct_",year_1,"_",year_2)]]<-adm0_diffs
        rm(adm0_diffs,est_yr1,est_yr2)


      # Admin 1 differences
        message("   admin 1")
        est_yr1<-adm1[year==year_1,][,draw_colnames,with=F] # Get just the draw columns, ordered the same way
        est_yr2<-adm1[year==year_2,][,draw_colnames,with=F]
        adm1_diffs<-100*(est_yr2-est_yr1)/est_yr1 # basically overlapping the data.tables as matrices to calculate across rows and columns simultaneously
        adm1_diffs[,year:=paste0("pct_",year_1,"_",year_2)]

        adm1_diffs<-cbind(adm1_ids,adm1_diffs) # Adding back on the ID columns
        adm1_diffs<-adm1_diffs[,list(lower=apply(.SD,FUN="lower_confint",MARGIN=1),
                                     mean=apply(.SD,FUN="mean",MARGIN=1),
                                     upper=apply(.SD,FUN="upper_confint",MARGIN=1),
                                     ADM1_CODE,year,indicator),.SDcols=draw_colnames]
        adm1_diffs<-merge(adm1_diffs,unique(sp_h[,list(ADM1_CODE,ADM1_NAME,ADM0_NAME)]),by="ADM1_CODE",all.x=T)
        adm1_diffs_list[[paste0("pct_",year_1,"_",year_2)]]<-adm1_diffs
        rm(adm1_diffs,est_yr1,est_yr2)


      # Admin 2 differences
        message("   admin 2")
        est_yr1<-adm2[year==year_1,][,draw_colnames,with=F]
        est_yr2<-adm2[year==year_2,][,draw_colnames,with=F]
        adm2_diffs<-100*(est_yr2-est_yr1)/est_yr1
        adm2_diffs[,year:=paste0("pct_",year_1,"_",year_2)]
        adm2_diffs<-cbind(adm2_ids,adm2_diffs) # Adding back on the ID columns

        adm2_diffs<-adm2_diffs[,list(lower=apply(.SD,FUN="lower_confint",MARGIN=1),
                                     mean=apply(.SD,FUN="mean",MARGIN=1),
                                     upper=apply(.SD,FUN="upper_confint",MARGIN=1),
                                     ADM2_CODE,year,indicator),.SDcols=draw_colnames]
        adm2_diffs<-merge(adm2_diffs,unique(sp_h[,list(ADM2_CODE,ADM2_NAME,ADM1_NAME,ADM0_NAME)]),by="ADM2_CODE",all.x=T)
        adm2_diffs_list[[paste0("pct_",year_1,"_",year_2)]]<-adm2_diffs
        rm(adm2_diffs,est_yr1,est_yr2)
    }

    adm0_diffs<-rbindlist(adm0_diffs_list);rm(adm0_diffs_list)
    adm1_diffs<-rbindlist(adm1_diffs_list);rm(adm1_diffs_list)
    adm2_diffs<-rbindlist(adm2_diffs_list);rm(adm2_diffs_list)
    message("Done creating admin-based % difference files.")

    save(adm0_diffs,adm1_diffs,adm2_diffs,file=diffs_tables_path)
    #save(adm0_diffs,adm1_diffs,file=diffs_tables_path)
    rm(adm0_diffs,adm1_diffs,adm2_diffs)

  } # closing if-calc-diffs

# Getting confidence intervals of the mean estimates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generating a weighted mean for the entire spatial extent (Africa, for now):
  all_regions <- regs[,lapply(.SD,weighted.mean,w=pop, na.rm=T),by=c("year","indicator"), .SDcols=draw_colnames]

  all_regions<-all_regions[,list(lower=apply(.SD,FUN="lower_confint",MARGIN=1),
                                 mean=apply(.SD,FUN="mean",MARGIN=1),
                                 upper=apply(.SD,FUN="upper_confint",MARGIN=1),
                                 year,indicator),.SDcols=draw_colnames]

  adm0<-adm0[,list(lower=apply(.SD,FUN="lower_confint",MARGIN=1),
                   mean=apply(.SD,FUN="mean",MARGIN=1),
                   upper=apply(.SD,FUN="upper_confint",MARGIN=1),
                   ADM0_CODE,year,indicator,pop),.SDcols=draw_colnames]
  adm0<-merge(adm0,unique(sp_h[,list(ADM0_CODE,ADM0_NAME)]),by="ADM0_CODE",all.x=T)

  adm1<-adm1[,list(lower=apply(.SD,FUN="lower_confint",MARGIN=1),
                   mean=apply(.SD,FUN="mean",MARGIN=1),
                   upper=apply(.SD,FUN="upper_confint",MARGIN=1),
                   ADM1_CODE,year,indicator,pop),.SDcols=draw_colnames]
  adm1<-merge(adm1,unique(sp_h[,list(ADM1_CODE,ADM1_NAME,ADM0_NAME)]),by="ADM1_CODE",all.x=T)

  adm2<-adm2[,list(lower=apply(.SD,FUN="lower_confint",MARGIN=1),
                   mean=apply(.SD,FUN="mean",MARGIN=1),
                   upper=apply(.SD,FUN="upper_confint",MARGIN=1),
                   ADM2_CODE,year,indicator,pop),.SDcols=draw_colnames]
  adm2<-merge(adm2,unique(sp_h[,list(ADM2_CODE,ADM2_NAME,ADM1_NAME,ADM0_NAME)]),by="ADM2_CODE",all.x=T)


  save(results,adm0,adm1,adm2,all_regions, file=results_tables_path)


  ## Plotting results for visual inspection:
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pdf(file=results_time_series_path,width=12,height=7) ## Starting a PDF of the results

  ## Subset only the indicators that you really want to show:
  plot_indicators<-results[final_set==1]$indicator_longname ## Only taking the ones in the "final" list

  ## A plot of all final indicators, together
  ggplot(all_regions[indicator %in% plot_indicators], aes(year, mean,fill=indicator,color=indicator))+
    geom_point()+
    geom_line(data=all_regions)+
    labs(title="All Regions, All Indicators")+
    geom_ribbon(data=all_regions,
                aes(ymin=lower,ymax=upper),
                alpha=0.3)

  ## Plot for each indicator, by itself
  for (i in plot_indicators){
    by_indicator<-ggplot(all_regions[indicator==i], aes(year, mean))+
      geom_point()+
      geom_line()+
      geom_ribbon(aes(ymin=lower,ymax=upper),
                  alpha=0.3)+
      labs(title=i)
    plot(by_indicator)
  }

  ## Plot for each country, all indicators.
  ## Making a table with the unique name/gaul code combinations to iterate through
  name_table<-unique(adm0[,list(ADM0_CODE,ADM0_NAME)])[order(ADM0_NAME)]

  for(code in unique(name_table$ADM0_CODE)){
    country_name<-name_table[ADM0_CODE==code,]$ADM0_NAME
    print(country_name)
    p<-ggplot(adm0[ADM0_CODE==code & indicator %in% plot_indicators,], aes(year, mean,fill=indicator,color=indicator))+
      geom_point()+
      geom_line(data=adm0[ADM0_CODE==code & indicator %in% plot_indicators,])+
      labs(title=country_name)+
      geom_ribbon(data=adm0[ADM0_CODE==code & indicator %in% plot_indicators,],
                  aes(ymin=lower,ymax=upper),
                  alpha=0.3)
    plot(p)
  }

  dev.off() ## Closing out the PDF.

