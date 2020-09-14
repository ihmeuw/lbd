#set regions and directory to search
reg_list <- c('cssa', 'wssa', 'essa', 'name', 'sssa')
run_date <- '2019_03_27_17_08_00'
in_dir <- '<<<< FILEPATH REDACTED >>>>'
setwd(in_dir)

#load in the child models
for (reg in reg_list){
  load(paste0('child_model_list_', reg, '_0.RData'))
  enet <- child_models$enet
  message(paste0("enet coefficients for ", reg))
  print(coef(enet, s = enet$cv_1se_lambda))
}
