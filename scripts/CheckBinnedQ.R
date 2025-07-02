fn <- file.path(path, "m84036_230702_205216_s2.MAS16S_Fwd_01--MAS16S_Rev_13.hifi_reads.fastq.gz")

foo <- plotQualityProfile(fn)
scores <- unique(foo$data$Score) |> sort()

print(scores)