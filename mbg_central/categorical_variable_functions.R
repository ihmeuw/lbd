#ras: a rasterlayer or brick
#the return is one layer or brick per non-NA value in the input raster. Don't use with floats.
create_categorical_raster = function(ras){
    
    #find the unique values in the raster
    uniq_vals = unique(as.vector(unique(ras)))
    uniq_vals = uniq_vals[!is.na(uniq_vals)]
    
    #create the new indicator rasters
    new_rasters = lapply(uniq_vals, function(x) ras == x)
    
    layer_name = names(ras)    
    
    #get overall object name(s)
    layer_name = unique(gsub('\\.[0-9]+$',"", layer_name))
        
    #set names both of the raster and its sublayers
    names(new_rasters) = paste0(layer_name,'_',uniq_vals)
    
    #iteratively set the object names
    for(rrr in 1:length(new_rasters)){
        names(new_rasters[[rrr]]) = rep(names(new_rasters)[rrr], dim(new_rasters[[rrr]])[3])
    }
        
    return(new_rasters)
        
}

#function to update the fixed effects variable
adjust_equation_categorical = function(effects_equation = fixed_effects, cat_var ='', cat_ras ){
    #format the raster names into an equation format
    eq_add = paste0(names(cat_ras), collapse = ' + ')
    effects_equation = sub(cat_var, eq_add, effects_equation)
    
    return(effects_equation)
    
}

#function to update the cov layers option
update_cov_layers = function(original_layers, new_layers){
    return(append(original_layers, new_layers))
}


