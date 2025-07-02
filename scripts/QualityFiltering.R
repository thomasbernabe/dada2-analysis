# filtered file names (single-end)
filtF <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))

# quality filtering for single end reads 
out <- filterAndTrim(fn, filtF, 
                     truncLen = 0,      
                     maxEE = 2,       
                     rm.phix = TRUE, 
                     compress = TRUE, 
                     multithread = TRUE)
