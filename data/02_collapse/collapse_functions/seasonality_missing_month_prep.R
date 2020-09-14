############################################
# Grab month information from ghdx and
# use the midpoint month if seasonality
#     microdata is unavailable
############################################

#input data
library(data.table)

source('<<<< FILEPATH REDACTED >>>>')
# Define the query
q <- pub_meta_template <- (
  "
  SELECT
    node.nid   as nid,
    fdft.field_time_value as start_date,
    fdft.field_time_value2 as end_date
  FROM
    ghdx.node node
  LEFT JOIN
    ghdx.field_data_field_time fdft
    ON node.nid = fdft.entity_id
  "
)

cooper <- execute_sql(statement=q, conn_def='ghdx') %>% setDT()

cooper$start_date <- as.character(cooper$start_date)
cooper$end_date <- as.character(cooper$end_date)

cooper <- cooper[!is.na(start_date)]
cooper <- cooper[nid %in% lri_data$nid,]

cooper[,start_year := substr(start_date, 1, 4)]
cooper[,end_year := substr(end_date, 1, 4)]
cooper[,start_month := substr(start_date, 6, 7)]
cooper[,end_month := substr(end_date, 6, 7)]

cooper[,time_period := paste0(start_month,'/',start_year,'-',end_month,'/',end_year)]
cooper <- cooper[,c('nid','time_period')]

mimo <- merge(lri_data, cooper, by = 'nid', all.x = T)
mimo <- unique(subset(mimo, is.na(scalar) & is.na(month))[,c('nid','country','survey_series','file_path','month','scalar','region_name','time_period')])
setnames(mimo, 'time_period','ghdx_date')

mimo = mimo[, start_mo := as.integer(substr(ghdx_date, 1, 2))]
mimo = mimo[, end_mo := as.integer(substr(ghdx_date, 9, 10))]
mimo = mimo[, num_months := ((end_mo - start_mo) %% 12) + 1]

mimo <- subset(mimo, !is.na(num_months))
mimo = mimo[,-c("month", "scalar")]

write.csv(mimo, '<<<< FILEPATH REDACTED >>>>', row.names = F)

#find midpoint month
mimo[, yr_span := 0]
mimo[end_mo < start_mo, yr_span := 1]
mimo[yr_span != 1, mid_mo := (start_mo + end_mo)/2]
mimo[yr_span == 1, mid_mo := (((12 - start_mo) + end_mo)/2) + start_mo - 12]
mimo[mid_mo <= 0, mid_mo := 12 - mid_mo]
mimo[mid_mo > 12, mid_mo := mid_mo - 12]

#read in scalars
season <- fread('<<<< FILEPATH REDACTED >>>>')
#merge ALL monthtly scalars on each survey row
seas <- season[, c("month", "scalar", "region_name")]
seas <- melt(seas, c("region_name", "month"), "scalar")
seas <- dcast(seas, region_name ~ month, value.var = "value")
mimo = merge(mimo, seas, by = "region_name")

#midpoint scalars
mimo$mid_mo = as.numeric(mimo$mid_mo)
mimo$num_months = as.numeric(as.character(mimo$num_months))
mimo = mimo[, mid_scalar := 0]
mimo = mimo[, avg_scalar := 0]



#maybe just run a for loop over each row? (with if statements)
for(i in 1:nrow(mimo)) {
  ## assign midpoint scalars
  #odd # month durations
  if(mimo[i, num_months] %% 2 != 0) { 
    mo = mimo[i, mid_mo]
    mo_col = as.character(mo)
    mo_val = mimo[i, get(mo_col)]
    mimo[i, mid_scalar := mo_val]
  }
  #even # month durations
  if(mimo[i, num_months] %% 2 == 0) {
    mo = mimo[i, mid_mo]
    mo_1 = as.character(mo - 0.5)
    mo_1 = ifelse(mo_1 == 0, "12", mo_1)
    mo_2 = as.character(mo + 0.5)
    mo_avg = (mimo[i, get(mo_1)] + mimo[i, get(mo_2)]) / 2
    mimo[i, mid_scalar := mo_avg]
  }
  
  ## assign average month scalars
  #non-year span surveys
  if(mimo[i, yr_span] != 1) {
    mos = c(mimo[i, start_mo]:mimo[i, end_mo])
    sclrs = vector(mode = "numeric", length = length(mos))
    for(m in 1:length(mos)) {
      mo_col = as.character(mos[m])
      sclrs[m] = mimo[i, get(mo_col)]
    }
    mimo = mimo[i, avg_scalar := mean(sclrs)]
  }
  #year span surveys
  if(mimo[i, yr_span] == 1) {
    mos_1 = c(mimo[i, start_mo]:12)
    mos_2 = c(1:mimo[i, end_mo])
    mos = c(mos_1, mos_2)
    sclrs = vector(mode = "numeric", length = length(mos))
    for(m in 1:length(mos)) {
      mo_col = as.character(mos[m])
      sclrs[m] = mimo[i, get(mo_col)]
    }
    mimo = mimo[i, avg_scalar := mean(sclrs)]
  }
}

#save data
write.csv(mimo, '<<<< FILEPATH REDACTED >>>>', row.names = F)
