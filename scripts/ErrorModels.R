library(ShortRead)
library(magrittr)
library(dplyr)

# function to subsample fastq files 
subsample_fastq <- function(input_files, n_reads = 100000, overwrite = TRUE) {
  
  subsampled_files <- gsub("\\.fastq\\.gz$", "_sub.fastq.gz", input_files)
  
  for (i in seq_along(input_files)) {
    reads <- readFastq(input_files[i])
    
    # check if file exists for filt_subsample_custom 
    if (file.exists(subsampled_files[i])) {
      if (overwrite) {
        file.remove(subsampled_files[i])
        cat("Removed existing file:", subsampled_files[i], "\n")
      } else {
        cat("File exists, skipping:", subsampled_files[i], "\n")
        next
      }
    }
    reads <- readFastq(input_files[i])
    
    if (length(reads) > n_reads) {
      set.seed(123)
      sample_idx <- sample(length(reads), n_reads)
      reads <- reads[sample_idx]
    }
    
    writeFastq(reads, subsampled_files[i], compress = TRUE)
  }
  
  return(subsampled_files)
}

# subsample filtered files 
filt_subsample <- subsample_fastq(filtF, n_reads = 50000)

# default loess error model  
err_loess <- learnErrors(filt_subsample, multithread = TRUE)
plotErrors(err_loess, nominalQ=TRUE)

# binned quality score error model
err_binned <- learnErrors(filt_subsample, 
                          errorEstimationFunction = makeBinnedQualErrfun(c(scores)), 
                          multithread=TRUE)
plotErrors(err_binned, nominalQ=TRUE)

# subsample for the modified loess (needs more than the previous two)
filt_subsample_custom <- subsample_fastq(filtF, n_reads = 100000, overwrite = TRUE)

# modified loess error function from Github
loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  
  # check if we have sufficient data
  total_observations <- sum(trans)
  if(total_observations < 1000) {
    warning("Insufficient data for error estimation: ", total_observations, " observations")
    return(NULL)
  }
  
  for (nti in c("A","C","G","T")) {
    for (ntj in c("A","C","G","T")) {
      if (nti != ntj) {
        
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        
        # check for sufficient data for this transition
        if(sum(tot) < 10) {
          # if insufficient data, use a simple error rate
          simple_rate <- sum(errs) / sum(tot)
          if(is.na(simple_rate) || simple_rate == 0) simple_rate <- 1e-7
          est <- rbind(est, rep(simple_rate, length(qq)))
          next
        }
        
        rlogp <- log10((errs+1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # remove rows with NA or infinite values
        df <- df[is.finite(df$rlogp) & df$tot > 0,]
        
        if(nrow(df) < 3) {
          # not enough data points for loess, use simple approach
          simple_rate <- sum(errs) / sum(tot)
          if(is.na(simple_rate) || simple_rate == 0) simple_rate <- 1e-7
          est <- rbind(est, rep(simple_rate, length(qq)))
          next
        }
        
        # try loess fitting with error handling
        tryCatch({
          # use smaller span if we have limited data
          span_val <- max(0.5, min(2, 3/nrow(df)))
          mod.lo <- loess(rlogp ~ q, df, weights = log10(df$tot + 1), span = span_val)
          pred <- predict(mod.lo, qq)
          
          # handle NA predictions
          if(all(is.na(pred))) {
            simple_rate <- sum(errs) / sum(tot)
            if(is.na(simple_rate) || simple_rate == 0) simple_rate <- 1e-7
            pred <- rep(log10(simple_rate), length(qq))
          } else {
            # fill in NA values
            maxrli <- max(which(!is.na(pred)))
            minrli <- min(which(!is.na(pred)))
            pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
            pred[seq_along(pred) < minrli] <- pred[[minrli]]
          }
          
          est <- rbind(est, 10^pred)
          
        }, error = function(e) {
          # fallback to simple error rate
          simple_rate <- sum(errs) / sum(tot)
          if(is.na(simple_rate) || simple_rate == 0) simple_rate <- 1e-7
          est <- rbind(est, rep(simple_rate, length(qq)))
          return(est)
        })
        
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # check if we have a valid matrix
  if(nrow(est) != 12 || ncol(est) != length(qq)) {
    warning("Error matrix has wrong dimensions: ", nrow(est), "x", ncol(est))
    return(NULL)
  }
  
  # set error rate bounds
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity without referencing undefined X40
  estorig <- est
  
  # ensure error rates don't increase with quality
  for(i in 1:nrow(est)) {
    for(j in 2:ncol(est)) {
      if(est[i,j] > est[i,j-1]) {
        est[i,j] <- est[i,j-1]
      }
    }
  }
  
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  
  # final validation
  if(any(is.na(err)) || any(is.infinite(err))) {
    warning("Error matrix contains NA or infinite values")
    return(NULL)
  }
  
  return(err)
}

# test custom error model to run 
tryCatch({
  err_custom <- learnErrors(filt_subsample_custom, 
                            errorEstimationFunction = loessErrfun_mod, 
                            multithread = TRUE)
  
  if(!is.null(err_custom)) {
    plotErrors(err_custom, nominalQ=TRUE)
  } else {
    print("Custom Loess did not work!")
  }
  
}, error = function(e) {
  err_custom <- learnErrors(filt_subsample_custom, multithread = TRUE)
  plotErrors(err_custom, nominalQ=TRUE)
})