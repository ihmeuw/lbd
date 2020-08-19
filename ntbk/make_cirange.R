setwd("<<<< FILEPATH REDACTED >>>>")

library(raster)

s_imp <- brick('sani_imp_total_uci_ras.tif') - brick('sani_imp_total_lci_ras.tif')
s_piped <- brick('sani_piped_uci_ras.tif') - brick('sani_piped_lci_ras.tif')
s_unimp <- brick('sani_unimp_cr_uci_ras.tif') - brick('sani_unimp_cr_lci_ras.tif')
s_od <- brick('sani_no_facility_uci_ras.tif') - brick('sani_no_facility_lci_ras.tif')

w_imp <- brick('water_imp_total_uci_ras.tif') - brick('water_imp_total_lci_ras.tif')
w_piped <- brick('water_piped_uci_ras.tif') - brick('water_piped_lci_ras.tif')
w_unimp <- brick('water_unimp_cr_uci_ras.tif') - brick('water_unimp_cr_lci_ras.tif')
w_od <- brick('water_no_facility_uci_ras.tif') - brick('water_no_facility_lci_ras.tif')


setwd("<<<< FILEPATH REDACTED >>>>")
writeRaster(s_imp,
  filename = 's_imp_cirange_unraked_2000_2017.tif',
  format = 'TGiff', overwrite = TRUE)

writeRaster(s_piped,
  filename = 's_piped_cirange_unraked_2000_2017.tif',
  format = 'GTiff', overwrite = TRUE)

writeRaster(s_unimp,
  filename = 's_unimp_cirange_unraked_2000_2017.tif',
  format = 'GTiff', overwrite = TRUE)

writeRaster(s_od,
  filename = 's_od_cirange_unraked_2000_2017.tif',
  format = 'GTiff', overwrite = TRUE)

writeRaster(w_imp,
  filename = 'w_imp_cirange_unraked_2000_2017.tif',
  format = 'GTiff', overwrite = TRUE)

writeRaster(w_piped,
  filename = 'w_piped_cirange_unraked_2000_2017.tif',
  format = 'GTiff', overwrite = TRUE)

writeRaster(w_unimp,
  filename = 'w_unimp_cirange_unraked_2000_2017.tif',
  format = 'GTiff', overwrite = TRUE)

writeRaster(w_od,
  filename = 'w_surface_cirange_unraked_2000_2017.tif',
  format = 'GTiff', overwrite = TRUE)
