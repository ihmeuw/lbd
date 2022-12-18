####ANC EXTRACTION CODE######################################################################

# It's purpose is to extract the ANC data from the read_epp_data() function by country and year, clean the data, and output csvs for futher analysis.
# The base set up was copied from the HIV collapse code. 
  
###############################################################################################

## LBD Base Set Up
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))
## Set repo
commondir      <- sprintf("<<<< FILEPATH REDACTED >>>>")
package_list <- c(t(read.csv(sprintf("<<<< FILEPATH REDACTED >>>>",commondir), header = FALSE)))
user <- Sys.info()["user"]
core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>")
if (!dir.exists(core_repo)) core_repo <- "<<<< FILEPATH REDACTED >>>>"
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>")

setwd(core_repo)

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, "<<<< FILEPATH REDACTED >>>>"))
mbg_setup(package_list = package_list, repos = core_repo)


library(reshape2)
library(dplyr)

source("<<<< FILEPATH REDACTED >>>>")

##################################################################################################
                         
countries <- as.list(c("AGO","BDI","CMR","DJI","COD","ETH","LSO","MLI","MOZ","NGA", "RWA","SOM","UGA","SSD", "GAB","TZA",
"GHA","ERI","LBR","CAF","GMB","GNB","SWZ","MDG","GIN","NER", "SDN","BWA","NAM","TGO","MAR",
"TCD","MWI","BFA","COG","BEN","CIV","ZMB","CPV","SLE","MRT","ZWE", "COM", "GNQ", "KEN", "SEN", "SLE", "ZAF", "HTI", "DOM", "KHM"))

yrs <- c("2019")
output <- list()
for (c in countries){
  output[[c]] <- as.list(c(c))
  for (y in yrs){
    #create directory for saving
    dir <- paste0("<<<< FILEPATH REDACTED >>>>")
    files <- list.files(dir, full.names = T)
    files <- as.list(grep("PJNZ", files, value = T))
    
    for(f in files) {
    if (file.exists(f)){
      data <- read_epp_data(f,c)
      for (pop in names(data)){
        
        anc_prev <- as.data.frame(data[[pop]]$anc.prev)
        anc_prev$site <- rownames(anc_prev)
        anc_In <- as.data.frame(data[[pop]]$anc.used)
        
        ##created anc_In for site identification 
        anc_In$site <- rownames(anc_prev)
        names(anc_In)[1]<-paste("In")
        anc_In$In <- as.numeric(anc_In$In)
        
        #continuing prevalence tidying
        anc_prev <- melt(anc_prev, id.vars=c("site"))
        names(anc_prev) <- c("site", "year", "prev")
        
        #merge"anc_In" lookup table to long format anc_prev
        anc_prev <- (merge(anc_In, anc_prev, by = 'site'))
        
        anc_n <- as.data.frame(data[[pop]]$anc.n)
        anc_n$site <- rownames(anc_n)
        anc_n <- melt(anc_n, id.vars=c("site"))
        names(anc_n) <- c("site", "year", "n")
        
        anc <- merge(anc_prev, anc_n, by = c("site", "year"))
        anc <- anc[which(!is.na(anc$prev)), ]
        anc$prev <- anc$prev * 100
        anc$Group <- pop
        anc
        
        data[[pop]] <- anc
        message(paste0("processed ", pop, "_", c, "_", y))
       
      }
      df <- ldply(data, data.frame)
      df$UNAIDS_file <- substr(f, 40, (nchar(f)-5))
      files[[f]] <- df
      
    }
    
    final <- ldply(files, data.frame)
    final[1:2] <- NULL
    final <- na.omit(final)
    colnames(final)[which(names(final) == "site")] <- "Site"
    colnames(final)[which(names(final) == "year")] <- "Year"
    colnames(final)[which(names(final) == "prev")] <- "Prev"
    colnames(final)[which(names(final) == "n")] <- "N"
    colnames(final)[which(names(final) == "in")] <- "In"
    col_order <- c("Group", "Site", "Year",
                   "Prev", "N", "In", "UNAIDS_file")
    final <- final[, col_order]
    output[[c]][[y]] <- final
                
                
    write.csv(output[[c]][[y]], 
              file = paste0("<<<< FILEPATH REDACTED >>>>") ,
              row.names = F)
    }
  }
}

