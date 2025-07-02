# test with onyl first file
first_file <- fn[1]
first_output <- file.path(nop_path, basename(first_file), "\n")

file.info(first_file)

# read a few lines of first file 
con <- gzfile(first_file, "r")
first_lines <- readLines(con, n = 20)
close(con)
print(first_lines)

prim_test <- removePrimers(first_file, first_output,
                           primer.fwd = F27_kinnex,
                           primer.rev = dada2:::rc(R1492_kinnex),
                           orient = TRUE,
                           verbose = TRUE)



# 
# Count total lines in the file
count_lines <- function(file) {
  con <- gzfile(file, "rt")
  line_count <- 0
  chunk_size <- 10000
  
  repeat {
    chunk <- readLines(con, n = chunk_size, warn = FALSE)
    if (length(chunk) == 0) break
    line_count <- line_count + length(chunk)
  }
  close(con)
  return(line_count)
}

total_lines <- count_lines(first_file)
cat("Total lines:", total_lines, "\n")
cat("Divisible by 4:", total_lines %% 4, "\n")
cat("Complete records:", floor(total_lines / 4), "\n")

if (total_lines %% 4 != 0) {
  cat("WARNING: File has incomplete FASTQ record at the end!\n")
  cat("Missing lines:", 4 - (total_lines %% 4), "\n")
}




# Calculate the number of lines for complete records only
complete_records <- floor(total_lines / 4)  # 1,213,447 complete records
complete_lines <- complete_records * 4      # 4,853,788 lines (instead of 4,853,790)

cat("Creating corrected file with", complete_records, "complete records\n")
cat("Using", complete_lines, "lines instead of", total_lines, "\n")

# Read only the complete records
con <- gzfile(first_file, "rt")
corrected_lines <- readLines(con, n = complete_lines, warn = FALSE)
close(con)

# Create corrected filename
corrected_file <- file.path(dirname(first_file), 
                            paste0(tools::file_path_sans_ext(basename(first_file)), 
                                   "_corrected.fastq.gz"))

# Write corrected file
con_out <- gzfile(corrected_file, "wt")
writeLines(corrected_lines, con_out)
close(con_out)

cat("Corrected file saved as:", corrected_file, "\n")

# Verify the corrected file
corrected_lines_check <- count_lines(corrected_file)
cat("Corrected file lines:", corrected_lines_check, "\n")
cat("Divisible by 4:", corrected_lines_check %% 4, "\n")