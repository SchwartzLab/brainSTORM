# Magrittr Pipe Operator
`%>%` <- magrittr::`%>%`

# Make Temporary directory for STORMseq processing
mkTmpDir <- function(){
    if(!dir.exists("./STORMtmp_dir")){dir.create("./STORMtmp_dir")}
}
rmTmpDir <- function(){
    if(dir.exists("./STORMtmp_dir")){unlink("./STORMtmp_dir", recursive = TRUE)}
}

# Create temporal bisulphite genome
bisGenome <- function(fastaGenome){
    mkTmpDir()
    tmpDir <- "./STORMtmp_dir"
    out_name <- file.path(tmpDir, "bisGenome.fa")
    gen <- Biostrings::readDNAStringSet(fastaGenome)
    gen_r <- lapply(gen, function(x){
        tmp <- stringr::str_replace_all(x, pattern = "C", replacement = "T")
    })
    seqinr::write.fasta(sequences = gen_r, names = names(gen_r), file.out = out_name, as.string = T)
}

# Make GTF from BED
mkGTF <- function(bedAnnotationm, outName = NULL){
    mkTmpDir()
    tmpDir <- file.path(getwd(), "STORMtmp_dir")
    if(is.null(outName)){outName <- file.path(tmpDir, "geneAnnot.gtf")}
    com <- paste("/apps/RH7U2/general/kentUtils/v377/bin/bedToGenePred",
                  bedAnnotation, "/dev/stdout | genePredToGtf file /dev/stdin", outName)
    system(com)
}

# bedToGenePred in.bed /dev/stdout | genePredToGtf file /dev/stdin out.gtf
