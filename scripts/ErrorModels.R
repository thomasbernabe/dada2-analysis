library(ShortRead)
library(magrittr)
library(dplyr)

# function to subsample fastq files 
subsample_fastq <- function(input_files, n_reads = 100000) {
  
  subsampled_files <- gsub("\\.fastq\\.gz$", "_sub.fastq.gz", input_files)
  
  for (i in seq_along(input_files)) {
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
err_default <- learnErrors(filt_subsample, multithread = TRUE)
plotErrors(err_default, nominalQ=TRUE)

# pacbio error model 
err_pacbio <- learnErrors(filt_subsample, 
                           errorEstimationFunction = PacBioErrfun, 
                           multithread = TRUE)
plotErrors(err_pacbio, nominalQ=TRUE)

# modified loess error function from Github
loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  
  for (nti in c("A","C","G","T")) {
    for (ntj in c("A","C","G","T")) {
      if (nti != ntj) {
        
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # modified loess with custom weights and span (Gulliem Salazar's solution)
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot), span = 2)
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # set error rate bounds
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # Enforce monotonicity
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  
  # return
  return(err)
}

# custom error model 
err_custom <- learnErrors(filt_subsample, 
                          errorEstimationFunction = loessErrfun_mod, 
                          multithread = TRUE)
plotErrors(err_custom, nominalQ=TRUE)

# binned quality score error model
err_binned <- learnErrors(filt_subsample, 
                          errorEstimationFunction = makeBinnedQualErrfun(c(scores)), 
                          multithread=TRUE)
plotErrors(err_binned, nominalQ=TRUE)