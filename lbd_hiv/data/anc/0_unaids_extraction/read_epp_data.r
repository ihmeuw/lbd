############################################################################
####  Function to read prevalence data used in EPP fitting (from .xml)  ####
############################################################################

# Description of the inputs:
# pjnz <- string file path to the pjnz file from UNAIDS that is being extracted
# loc <- string ISO3 code for the country being extracted.  This parameter was added in this instance and is used to govern some processing exceptions.


read_epp_data <- function(pjnz, loc = loc){
  xmlfile <- grep(".xml", unzip(pjnz, list = TRUE)$Name, value = TRUE)
  con <- unz(pjnz, xmlfile)
  epp.xml <- scan(con, "character", sep = "\n")
  close(con)
  
  if (!require("XML", quietly = TRUE))
    stop("read_epp_data() requires the package 'XML'. Please install it.", call. = FALSE)
  
  obj <- xmlTreeParse(epp.xml)
  
  r <- xmlRoot(obj)[[1]]
  eppSetChildren.idx <- which(xmlSApply(r, xmlAttrs) == "eppSetChildren")
  country <- xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetCountry")]][[1]])
  
  epp.data <- list() # declare list to store output
  attr(epp.data, loc) <- country
  
  for (eppSet.idx in 1:xmlSize(r[[eppSetChildren.idx]])) {
    print(eppSet.idx)
    tmp.eppSet <- r[[eppSetChildren.idx]][[eppSet.idx]][[1]]
    n.iter <- 1
    if (xmlSize(xmlSApply(tmp.eppSet, xmlAttrs)) > 0 & xmlSize(which(xmlSApply(tmp.eppSet, xmlAttrs) == "siteNames")) == 0) {
      newSetChildren.idx <- which(xmlSApply(tmp.eppSet, xmlAttrs) == "eppSetChildren")
      n.iter <- xmlSize(tmp.eppSet[[newSetChildren.idx]])
    }
    used.indices <- c()
    if (length(which(xmlSApply(tmp.eppSet, xmlAttrs) == "groupToAssignTo")) > 0) {
      hold <- tmp.eppSet[[which(xmlSApply(tmp.eppSet, xmlAttrs) == "groupToAssignTo")]]
      if (any(xmlSApply(hold, xmlAttrs) == "epp2011.core.sets.ProjectionSet")) {
        if (xmlSize(hold[[which(xmlSApply(hold, xmlAttrs) == "epp2011.core.sets.ProjectionSet")]]) > 1) {
          n.iter <- n.iter + length(which(xmlSApply(tmp.eppSet, xmlAttrs) == "groupToAssignTo"))
        }        
      }
    }
    for (eppSet.idx.2 in 1:n.iter) {
      print(eppSet.idx.2)
      eppSet <- tmp.eppSet
      if (xmlSize(xmlSApply(tmp.eppSet, xmlAttrs)) > 0 & xmlSize(which(xmlSApply(tmp.eppSet, xmlAttrs) == "siteNames")) == 0)
        eppSet <- tmp.eppSet[[newSetChildren.idx]][[eppSet.idx.2]][[1]]
      
      if (xmlSize(eppSet) == 0) {
        
        tmp_len <- 0
        tmp_index <- 0
        obj_size <- 0
        while (obj_size == 0) {
          tmp_index <- tmp_index + 1
          if (tmp_index != eppSet.idx & !(tmp_index %in% used.indices)) {
            eppSet <- r[[eppSetChildren.idx]][[tmp_index]][[1]]
            tmp_len <- length(which(xmlSApply(eppSet, xmlAttrs) == "groupToAssignTo"))
            if (tmp_len > 0) {
              tmp_eppSet <- eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "groupToAssignTo")]]
              set_location <- ifelse(xmlSize(tmp_eppSet) == 1, 1, which(xmlSApply(tmp_eppSet, xmlAttrs) == "epp2011.core.sets.ProjectionSet"))
              tmp_eppSet <- tmp_eppSet[[set_location]]
              obj_size <- xmlSize(tmp_eppSet)
            }
          }
        }
        used.indices <- c(used.indices, tmp_index)
        
        eppSet <- eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "groupToAssignTo")]]
        set_location <- ifelse(xmlSize(eppSet) == 1, 1, which(xmlSApply(eppSet, xmlAttrs) == "epp2011.core.sets.ProjectionSet"))
        eppSet <- eppSet[[set_location]]
      }
      if (eppSet.idx.2 > 1 & loc != "PNG") {
        eppSet <- eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "groupToAssignTo")]]
        set_location <- ifelse(xmlSize(eppSet) == 1, 1, which(xmlSApply(eppSet, xmlAttrs) == "epp2011.core.sets.ProjectionSet"))
        eppSet <- eppSet[[set_location]]
      }
      eppName <- xmlToList(eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "name")]][["string"]])
      print(eppName)
      
      # ANC data
      siteNames.idx <- which(xmlSApply(eppSet, xmlAttrs) == "siteNames")
      siteSelected.idx <- which(xmlSApply(eppSet, xmlAttrs) == "siteSelected")
      survData.idx <- which(xmlSApply(eppSet, xmlAttrs) == "survData")
      survSampleSizes.idx <- which(xmlSApply(eppSet, xmlAttrs) == "survSampleSizes")
      
      siteNames <- xmlSApply(eppSet[[siteNames.idx]][[1]], xmlSApply, xmlToList, FALSE)
      siteIdx <- as.numeric(xmlSApply(eppSet[[siteNames.idx]][[1]], xmlAttrs)) ## 0 based
      
      nsites <- length(siteNames)
      nANCyears <- max(as.integer(xmlSApply(eppSet[[survData.idx]][["array"]][[1]][[1]], xmlAttrs))) + 1
      
      # ANC site used
      anc.used <- rep(FALSE, nsites)
      anc.used[as.integer(xmlSApply(eppSet[[siteSelected.idx]][[1]], xmlAttrs)) + 1] <- as.logical(xmlSApply(eppSet[[siteSelected.idx]][[1]], xmlSApply, xmlToList, FALSE))
      
      # ANC prevalence
      anc.prev <- matrix(NA, nsites, nANCyears)
      rownames(anc.prev) <- siteNames
      colnames(anc.prev) <- 1985 + 0:(nANCyears - 1)
      for (clinic.idx in 1:nsites) {
        clinic <- eppSet[[survData.idx]][["array"]][[clinic.idx]][[1]]
        prev <- as.numeric(xmlSApply(clinic, xmlSApply, xmlToList, FALSE))
        idx <- as.integer(xmlSApply(clinic, xmlAttrs)) + 1
        anc.prev[clinic.idx, idx] <- prev
      }
      anc.prev[is.na(anc.prev)] <- 0.0 
      anc.prev[anc.prev == -1] <- NA
      anc.prev <- anc.prev/100
      
      # ANC sample sizes
      anc.n <- matrix(NA, nsites, nANCyears)
      rownames(anc.n) <- siteNames
      colnames(anc.n) <- 1985 + 0:(nANCyears - 1)
      for (clinic.idx in 1:nsites) {
        clinic <- eppSet[[survSampleSizes.idx]][["array"]][[clinic.idx]][[1]]
        n <- as.numeric(xmlSApply(clinic, xmlSApply, xmlToList, FALSE))
        idx <- as.integer(xmlSApply(clinic, xmlAttrs)) + 1
        anc.n[clinic.idx, idx] <- n
      }
      anc.n[anc.n == -1] <- NA
      

      hhsUsed.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyIsUsed")
      hhsHIV.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyHIV")
      hhsSampleSize.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveySampleSize")
      hhsSE.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyStandardError")
      hhsYear.idx <- which(xmlSApply(eppSet, xmlAttrs) == "surveyYears")

      nhhs <- max(try(xmlSize(eppSet[[hhsYear.idx]][[1]])),
                  try(xmlSize(eppSet[[hhsHIV.idx]][[1]])),
                  try(xmlSize(eppSet[[hhsSE.idx]][[1]])),
                  try(xmlSize(eppSet[[hhsUsed.idx]][[1]])))

      hhs <- data.frame(year = rep(NA, nhhs), prev = rep(NA, nhhs), se = rep(NA, nhhs), n = rep(NA, nhhs), used = rep(NA, nhhs))

      hhs$year[as.integer(xmlSApply(eppSet[[hhsYear.idx]][[1]], xmlAttrs)) + 1] <- as.numeric(xmlSApply(eppSet[[hhsYear.idx]][[1]], xmlSApply, xmlToList, FALSE))
      hhs$prev[as.integer(xmlSApply(eppSet[[hhsHIV.idx]][[1]], xmlAttrs)) + 1] <- as.numeric(xmlSApply(eppSet[[hhsHIV.idx]][[1]], xmlSApply, xmlToList, FALSE))/100
      hhs$se[as.integer(xmlSApply(eppSet[[hhsSE.idx]][[1]], xmlAttrs)) + 1] <- as.numeric(xmlSApply(eppSet[[hhsSE.idx]][[1]], xmlSApply, xmlToList, FALSE))/100
      hhs$n <- 20
      hhs$used[as.integer(xmlSApply(eppSet[[hhsUsed.idx]][[1]], xmlAttrs)) + 1] <- as.logical(xmlSApply(eppSet[[hhsUsed.idx]][[1]], xmlSApply, xmlToList, FALSE))

      hhs <- subset(hhs, !is.na(prev))
    
      epp.data[[eppName]] <- list(country = country,
                                  region = eppName,
                                  anc.used = anc.used,
                                  anc.prev = anc.prev,
                                  anc.n = anc.n,
                                  hhs = hhs)
   }
 }
  

  class(epp.data) <- "eppd"
  
  return(epp.data)
}

