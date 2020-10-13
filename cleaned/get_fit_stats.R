
get_fit_stats <- function(dat_all, pred_raster, orig_dat, mo = "00") {
  monthly <- FALSE
  if (mo != "00") monthly <- TRUE
  points <- orig_dat[which(orig_dat$shape_type == 'point'), c('lat', 'long', 'shape_type')]
  
  dat0.pts <- dat_all[dat_all$PA==0, c('long', 'lat')]
  dat0.preds <- extract(preds_mean, dat0.pts)
  
  dat1.pts <- dat_all[dat_all$PA==1, c('long', 'lat')]
  dat1.preds <- extract(preds_mean, dat1.pts)
  
  pos_mean_preds <- na.omit(dat1.preds)
  neg_mean_preds <- na.omit(dat0.preds)
  
  sensitivity <- c()
  fpr <- c()
  dif = 0.01
  for (thresh in seq(0.00, 1, dif)){
    tp <- length(pos_mean_preds[which(pos_mean_preds >= thresh)])
    fn <- length(pos_mean_preds[which(pos_mean_preds < thresh)])
    tn <- length(neg_mean_preds[which(neg_mean_preds <= thresh)])
    fp <- length(neg_mean_preds[which(neg_mean_preds > thresh)])
    sensitivity <- append(sensitivity, tp / (tp + fn))
    fpr <- append(fpr, 1 - (tn / (fp + tn)))
  }
  
  distance <- c()
  for (i in 1:length(fpr)){
    distance <- append(distance, dist(rbind(c(fpr[i], sensitivity[i]), c(0, 1))))
  }
  opt_thresh <- max(seq(0.00, 1, dif)[which(distance == min(distance))])
  auc <- integrate.xy(fpr, sensitivity, 0, 1)
  thresh <- opt_thresh
  tp <- length(pos_mean_preds[which(pos_mean_preds >= thresh)])
  fn <- length(pos_mean_preds[which(pos_mean_preds < thresh)])
  tn <- length(neg_mean_preds[which(neg_mean_preds <= thresh)])
  fp <- length(neg_mean_preds[which(neg_mean_preds > thresh)])
  rmse <- sqrt(mean((1-pos_mean_preds)^2))
  error_rate <- (fp+fn)/(tp+tn+fp+fn)
  
  fit_stats <- data.frame(opt_thresh, auc, rmse, error_rate)
  return(fit_stats)
}
