# Magrittr Pipe Operator
`%>%` <- magrittr::`%>%`

# Create temporal bisulphite genome
bisGenome <- function(fastaGenome){
    out_name <- file.path(tempdir(), "bisGenome.fa")
    gen <- Biostrings::readDNAStringSet(genome)
    gen_r <- lapply(gen, function(x){
        tmp <- stringr::str_replace_all(x, pattern = "C", replacement = "T")
    })
    seqinr::write.fasta(sequences = gen_r, names = names(gen_r), file.out = out_name, as.string = T)
}

