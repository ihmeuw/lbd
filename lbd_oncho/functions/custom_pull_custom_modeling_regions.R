#' @title Pull custom modeling regions
#'
#' @description Define modeling regions that are not simple combinations of the
#'   default MBG regions (in other words, regions that are not combinations of
#'   four-letter MBG regions such as "wssa" or "seas+ocea" or ISO-3 codes such
#'   as 'ZAF' or 'CHN').
#'
#' @param custom_region character vector of custom modeling regions
#'
#' @return Returns a named list of custom modeling regions with associated
#'    "standard" (non-custom) modeling regions that can be directly interpreted
#'    by get_adm0_codes().
#'
pull_custom_modeling_regions <- function(custom_regions){
  custom_regions <- tolower(custom_regions)
  
  # FULL LIST OF ALL REFERENCE REGIONS
  # If you need to add a new custom region, add it to this list
  ref_reg_list <- list(
    'africa'          = 'noaf+essa+wssa+cssa+sssa-yem',
    'middle_east'     = 'mide+stan-pak',
    'eastern_europe'  = "blr+est+lva+ltu+mda+ukr",
    'latin_america'   = 'caca+trsa+ansa',
    'south_asia'      = 'soas+chn_d2+pak',
    'central_america' = 'caca',
    'south_america'   = 'ansa+trsa',
    # se_asia was historically inclusive of East Asia + SE Asia
    'se_asia'         = 'eaas+seas+ocea+png',
    
    #data coverage regions
    'africa_dcp' = 'noaf+essa+wssa+cssa+sssa+yem',
    'middle_east_dcp' = 'mide+stan-yem-pak',
    'latin_america_dcp' = 'latin_america+cub',
    'south_asia_dcp' = 'south_asia-mdv-syc',
    'se_asia_dcp' = 'eaas+seas+png+idn+phl+tls+mys+twn',
    
    'stage1' = 'noaf+essa+wssa+cssa+sssa-yem',
    # ONLY stage 2 countries (not inclusive of Stage 1)
    'stage2' = 'ansa+caca+stan+eaas+mide+ocea+soas+seas+trsa+yem',
    'stage3' = 'all-stage1-stage2',
    
    # India + all disputed territories it claims: Jammu & Kashmir, Arunachal
    #   Pradesh, Aksay Chin, Demjok, Tirpani, Bara Hotii, & Samdu
    'ind_all_territories' = 'ind+ind_d1+ind_d2+ind_d3+chn_d2',
    
    # Getting into modeler-defined regions
    'cssa_diarrhea'  = 'cssa-ago',
    'sssa_diarrhea'  = 'sssa+ago-swz',
    'essa_diarrhea'  = 'essa+swz-ken_d1-yem',
    'cssa_diarrhea2' = 'cssa+cmr+tcd-ago-cod',
    'essa_diarrhea2' = 'essa+sdn+swz-ken_d1-yem',
    'name_diarrhea2' = 'noaf-esh-egy-egy_d1-sdn-sdn_d2',
    'sssa_diarrhea2' = 'sssa+ago-swz',
    'wssa_diarrhea2' = 'wssa+cmr+tcd-cmr-tcd',
    
    'essa_edu' = 'dji+eri+eth+som+ssd+sdn',
    'cssa_edu' = 'cssa+essa+nga+swz+cmr-dji-eri-eth-som-ssd-sdn-ken_d1-yem',
    'name_edu' = 'noaf-egy-egy_d1-sdn-sdn_d2-esh',
    'sssa_edu' = 'sssa-swz',
    'wssa_edu' = 'wssa-cmr-nga',
    
    'essa_sdn' = 'essa+sdn-ken_d1-yem',
    
    'cessa' = 'cssa+essa-ken_d1-yem',
    'cwssa' = 'cssa+wssa',
    
    'namelite' = 'noaf-dza-lby-tun-sdn_d1-sdn_d2-egy_d1-esh',
    
    # Region 'cessa2' originally contained American Samoa--this has been dropped
    'cessa2' = 'cssa+essa-ago-zmb-mwi-moz-ken_d1-yem',
    'sssa2'  = 'sssa+ago+zmb+mwi+moz',
    
    'cssa_cam'   = 'cssa+cmr',
    'wssa_nocam' = 'wssa-cmr',
    
    'name_hi' = 'noaf-sdn-sdn_d2-egy_d1-esh',
    'essa_hi' = 'zmb+ken',
    'essa_lo' = 'essa+sdn-zmb-ken-ken_d1-yem',
    'cssa_hi' = 'gnq+cog+gab',
    'cssa_lo' = 'cssa-gnq-cog-gab',
    'wssa_hi' = 'gha',
    'wssa_lo' = 'wssa-gha-mrt',
    'sssa_hi' = 'bwa+nam+zaf',
    'sssa_lo' = 'sssa-bwa-nam-zaf',
    
    'essa_hilo' = 'essa+lso+sdn+swz+zwe-ken_d1-yem',
    
    # Vaccines
    'vax_soas' = 'soas+pak-syc',
    'vax_seas' = 'seas+ocea-asm-fji-kir-wsm-ton',
    'vax_eaas' = 'eaas',
    
    # Need to move these back to vax_seas when 0long issue fixed
    # Not using this for now
    'vax_seas_0long' = 'asm+fji+kir+wsm+ton',
    
    'vax_caeu' = 'arm+aze+geo+kgz+mda+tjk+tkm+ukr+uzb',
    
    'vax_crbn' = 'cub+dma+dom+grd+hti+jam+lca+vct+vir',
    'vax_ctam' = 'blz+col+cri+slv+gtm+hnd+mex+nic+pan+ven',
    'vax_ansa' = 'ansa-col-ven',
    'vax_trsa' = 'trsa',
    
    'vax_name' = 'noaf+mide+afg+omn',
    'vax_cssa' = 'cssa',
    'vax_essa' = 'essa+syc',
    'vax_sssa' = 'sssa',
    'vax_wssa' = 'wssa',
    
    # For modelers who are still using the defunct NAME region
    'name_historic' = 'noaf-esh',
    
    # Custom regions for Diarrhea/ORT/WASH/LRI/HAP
    'dia_afr_horn' = 'dji+eri+eth+sdn+som+ssd+yem',
    'dia_cssa' = 'ago+caf+cod+cog+gab+gnq+stp',
    'dia_wssa' = 'ben+bfa+civ+cmr+cpv+gha+gin+gmb+gnb+lbr+mli+mrt+ner+nga+sen+sle+tcd+tgo',
    'dia_name' = 'dza+egy+esh+lby+mar+tun',
    'dia_sssa' = 'bwa+nam+zaf',
    'dia_mcaca' = 'blz+cri+cub+dma+dom+grd+gtm+hnd+hti+jam+lca+mex+nic+pan+slv+vct',
    'dia_s_america' = 'bol+bra+col+ecu+guf+guy+per+pry+sur+tto+ven',
    'dia_central_asia' = 'kgz+tjk+tkm+uzb',
    'dia_chn_mng' = 'chn+mng',
    'dia_se_asia' = 'khm+lao+mmr+mys+tha+vnm',
    'dia_malay' = 'idn+phl+png+tls',
    'dia_south_asia' = 'bgd+btn+ind+lka+npl+pak',
    'dia_mid_east' = 'afg+irn+irq+jor+pse+syr',
    'dia_essa' = 'bdi+com+ken+lso+mdg+moz+mwi+rwa+swz+syc+tza+uga+zmb+zwe',
    'dia_oceania' = 'asm+fji+fsm+kir+mhl+slb+ton+vut+wsm',
    
    # Custom regions for Focal 3
    'lf_s_asia' = 'bgd+ind+lka+npl',
    'lf_se_asia' = 'brn+idn+khm+lao+mmr+mys+phl+png+tha+tls+vnm',
    'lf_s_se_asia' = 'soas+seas+ocea+idn+phl+png+lka-btn-mdv-syc+brn',
    'lf_se_asia_pacific' = 'brn+idn+khm+lao+mmr+mys+phl+png+tha+tls+vnm+asm+cok+fsm+fji+pyf+kir+mhl+ncl+niu+plw+wsm+ton+tuv+vut+wlf',
    'lf_endem_afr' = 'essa+wssa+cssa+sdn+zwe+egy+yem-cpv',
    'lf_hispaniola' = 'hti+dom',
    'lf_non_mbg' = 'asm+bra+cok+fsm+fji+pyf+guy+kir+mdv+mhl+ncl+niu+plw+wsm+ton+tuv+vut+wlf',
    'oncho_endem_afr' = 'essa+wssa+cssa+sdn+yem-cpv-zmb-mdg-mrt-com-dji-eri'
    
  )
  # Warn if there are any custom regions not in the reference list
  missing_regions <- custom_regions[ !(custom_regions %in% names(ref_reg_list)) ]
  if( length(missing_regions) > 0 ){
    message(paste0('WARNING: The following custom regions are not defined: ',
                   paste(missing_regions, collapse=','))
    )
  }
  # Return a named list of all custom regions
  custom_regions_list <- ref_reg_list[ names(ref_reg_list) %in% custom_regions ]
  return(custom_regions_list)
}
