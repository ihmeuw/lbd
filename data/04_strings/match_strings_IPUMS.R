library(data.table)
library(feather)
library(parallel)
library(doParallel)
files <- list.files('<<<< FILEPATH REDACTED >>>>', 
                    pattern = '.feather', 
                    full.names = TRUE)

dt <- data.table()

message("Make cluster")
cl <- makeCluster(10)
clusterEvalQ(cl, .libPaths('<<<< FILEPATH REDACTED >>>>'))
message("Register cluster")
registerDoParallel(cl)
message("Start foreach")
#Read in each .dta file in parallel - returns a list of data frames
top <- foreach(i=1:length(files), .packages = c('feather')) %dopar% {
  dta <- read_feather(files[i])
  dta <- unique(dta[,c('nid', 'iso3', 't_type', 'sewage')])
  return(dta)
}
message("Foreach finished")
message("Closing cluster")
stopCluster(cl)

topics <- rbindlist(top, fill=T, use.names=T)
colnames(topics)[colnames(topics)=="t_type"] <- "toilet"
write.csv(topics, '<<<< FILEPATH REDACTED >>>>', row.names = FALSE)


defs <- read.csv('<<<< FILEPATH REDACTED >>>>')

test <- as.data.table(topics)
test <- test[!(nid %in% defs$nid)]


defs <- rbind(defs, test, fill = TRUE)
write.csv(defs, '<<<< FILEPATH REDACTED >>>>', row.names = FALSE)
