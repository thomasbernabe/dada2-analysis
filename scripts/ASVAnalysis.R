library(dplyr)
library(tidyr)

# run denoising 
dada_loess <- dada(filt_subsample, err=err_loess, multithread = TRUE)
dada_custom <- dada(filt_subsample_custom, err=err_custom, multithread = TRUE)
dada_binned <- dada(filt_subsample, err=err_binned, multithread = TRUE)

# create sequence tables 
seqtab_loess <- makeSequenceTable(dada_loess)
seqtab_custom <- makeSequenceTable(dada_custom)
seqtab_binned <- makeSequenceTable(dada_binned)

general_stats <- data.frame(
  Error_Model = c("Default Loess", "Custom Loess", "Binned Quality"),
  Total_ASVs = c(ncol(seqtab_loess), ncol(seqtab_custom), ncol(seqtab_binned)),
  Total_Reads = c(sum(seqtab_loess), sum(seqtab_custom), sum(seqtab_binned)),
  Mean_ASV_Length = c(
    round(mean(nchar(colnames(seqtab_loess)))),
    round(mean(nchar(colnames(seqtab_custom)))),
    round(mean(nchar(colnames(seqtab_binned))))
  )
)

# pick first sample 
sample_index <- 1
sample_name <- rownames(seqtab_loess)[sample_index]

# extract ASV data for this sample 
sample_loess <- seqtab_loess[sample_index, ]
sample_custom <- seqtab_custom[sample_index, ]
sample_binned <- seqtab_binned[sample_index, ]

# remove any zero abundance ASVs in this sample
sample_loess <- sample_loess[sample_loess > 0]
sample_custom <- sample_custom[sample_custom > 0]
sample_binned <- sample_binned[sample_binned > 0]

asv_count <- data.frame(
  Error_Model = c("Default Loess", "Custom Loess", "Binned Quality"),
  ASVs_in_Sample = c(length(sample_loess), length(sample_custom), length(sample_binned)),
  Total_Reads_in_Sample = c(sum(sample_loess), sum(sample_custom), sum(sample_binned))
)

# detailed ASV breakdown for sample 
analyze_asv <- function(sample_data, model_name) {
  asv_sort <- sort(sample_data, decreasing = TRUE) # sort by descending
  
  asv_info <- data.frame(
    ASV_Rank = 1:length(asv_sort),
    ASV_Sequence = names(asv_sort),
    ASV_Length = nchar(names(asv_sort)),
    Read_Count = as.numeric(asv_sort),
    Relative_Abundance = round(as.numeric(asv_sort) / sum(asv_sort) * 100, 2)
  )
  
  # model identifier ? test first 
  asv_info$Error_Model <- model_name
  
  return(asv_info)
}

# analyze each model 
asv_loess <- analyze_asv(sample_loess, "Default Loess")
asv_custom <- analyze_asv(sample_custom, "Custom Loess")
asv_binned <- analyze_asv(sample_binned, "Binned Quality")

# comparing selected samples across the 3 models 
count_asvs_per_sample <- function(seqtab, model_name) {
  asv_count_sample <- apply(seqtab, 1, function(x) sum(x > 0))
  read_count_sample <- rowSums(seqtab)
  
  result <- data.frame(
    Sample = names(asv_count_sample),
    Error_Model = model_name,
    ASVs = asv_count_sample,
    Total_Reads = read_count_sample
  )
  
  return(result)
}

# get the ASV counts for all samples 
all_samples_loess <- count_asvs_per_sample(seqtab_loess, "Default Loess")
all_samples_custom <- count_asvs_per_sample(seqtab_custom, "Custom Loess")
all_samples_binned <- count_asvs_per_sample(seqtab_binned, "Binned Quality")

combined_samples <- rbind(all_samples_loess, all_samples_custom, all_samples_binned)

asv_comparison <- combined_samples %>%
  select(Sample, Error_Model, ASVs) %>%
  pivot_wider(names_from = Error_Model, values_from = ASVs)

print("ASV Counts per Sample:")
print(asv_comparison)

# summary statistics 
summary_stats <- combined_samples %>%
  group_by(Error_Model) %>%
  summarise(
    Mean_ASVs_per_Sample = round(mean(ASVs), 1),
    Median_ASVs_per_Sample = median(ASVs),
    Min_ASVs = min(ASVs),
    Max_ASVs = max(ASVs),
    SD_ASVs = round(sd(ASVs), 1),
    .groups = 'drop'
  )

print(summary_stats)

if(!dir.exists(path.out)) dir.create(path.out)

# create csv files with stats

# overall comparison
write.csv(general_stats, file.path(path.out, "overall_asv_comparison.csv"), row.names = FALSE)
# per sample comparison
write.csv(asv_comparison, file.path(path.out, "asv_counts_per_sample.csv"), row.names = FALSE)
# save summary stats
write.csv(summary_stats, file.path(path.out, "asv_summary_statistics.csv"), row.names = FALSE)

sample_details <- rbind(
  asv_loess[, c("Error_Model", "ASV_Rank", "ASV_Length", "Read_Count", "Relative_Abundance")],
  asv_custom[, c("Error_Model", "ASV_Rank", "ASV_Length", "Read_Count", "Relative_Abundance")],
  asv_binned[, c("Error_Model", "ASV_Rank", "ASV_Length", "Read_Count", "Relative_Abundance")]
)

write.csv(sample_details, 
          file.path(path.out, paste0("detailed_asvs_", sample_name, ".csv")), 
          row.names = FALSE)

print("\n Summary")
print("Expected for ATCC mock community: ~20 species (ASVs)")
print(paste("Default Loess found:", general_stats$Total_ASVs[1], "total ASVs"))
print(paste("Custom Loess found:", general_stats$Total_ASVs[2], "total ASVs")) 
print(paste("Binned Quality found:", general_stats$Total_ASVs[3], "total ASVs"))