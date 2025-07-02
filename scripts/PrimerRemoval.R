F27_kinnex <- "AGRGTTYGATYMTGGCTCAG"
R1492_kinnex <- "RGYTACCTTGTTACGACTT"

nop_path <- file.path(path, "noprimers")
if (!dir.exists(nop_path)) dir.create(nop_path)

# remove the primers
nop <- file.path(nop_path, basename(fn))
prim <- removePrimers(fn, nop, 
                      primer.fwd = F27_kinnex,
                      primer.rev = dada2:::rc(R1492_kinnex),
                      orient = TRUE,
                      verbose = TRUE)

# debugging first 
prim <- removePrimers(fn[1], nop[1], 
                      primer.fwd = F27_kinnex,
                      primer.rev = dada2:::rc(R1492_kinnex),
                      orient = TRUE,
                      verbose = TRUE)


# create output path for corrected file
corrected_output <- file.path(nop_path, 
                              paste0(tools::file_path_sans_ext(basename(first_file)), 
                                     "_corrected.fastq.gz"))