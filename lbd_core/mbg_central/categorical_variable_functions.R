#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ras a rasterlayer or brick
#' @return one layer or brick per non-NA value in the input raster. Don't use with floats.
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname create_categorical_raster
#' @export
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


#' @title FUNCTION_TITLE
#' @description function to update the fixed effects variable
#' @param effects_equation PARAM_DESCRIPTION, Default: fixed_effects
#' @param cat_var PARAM_DESCRIPTION, Default: ''
#' @param cat_ras PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname adjust_equation_categorical
#' @export
adjust_equation_categorical = function(effects_equation = fixed_effects, cat_var ='', cat_ras ){
    #format the raster names into an equation format
    eq_add = paste0(names(cat_ras), collapse = ' + ')
    effects_equation = sub(cat_var, eq_add, effects_equation)
    
    return(effects_equation)
    
}


#' @title FUNCTION_TITLE
#' @description unction to update the cov layers option
#' @param original_layers PARAM_DESCRIPTION
#' @param new_layers PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname update_cov_layers
#' @export
update_cov_layers = function(original_layers, new_layers){
    return(append(original_layers, new_layers))
}


