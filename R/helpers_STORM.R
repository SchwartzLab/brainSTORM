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
    outName <- paste0(strsplit(fastaGenome, split = "/") %>% unlist %>% tail(1),
                      ".bisulp")
    out_name <- file.path(getwd(), "STORMtmp_dir", outName)
    gen <- Biostrings::readDNAStringSet(fastaGenome)
    gen_r <- lapply(gen, function(x){
        tmp <- stringr::str_replace_all(x, pattern = "C", replacement = "T")
    })
    seqinr::write.fasta(sequences = gen_r, names = names(gen_r), file.out = out_name, as.string = T)
    out_name
}

# Make GTF from BED
mkGTF <- function(bedAnnotation, outName = NULL, source = "user"){
    mkTmpDir()
    tmpDir <- file.path(getwd(), "STORMtmp_dir")
    if(is.null(outName)){outName <- file.path(tmpDir, "tmp_geneAnnot.gtf")}
    com <- paste0("/apps/RH7U2/general/kentUtils/v377/bin/bedToGenePred ",
                  bedAnnotation, " /dev/stdout | /apps/RH7U2/general/",
                  "kentUtils/v377/bin/genePredToGtf file /dev/stdin ", outName)
    system(com)
    tmp <- read.delim(outName, header = FALSE)
    tmp[,2] <- source
    write.table(tmp, file = outName, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}

# STAR genome generator
mkSTARgenome <- function(fastaGenome, bedAnnotation = NULL, outDir = NULL,
                         nCores = 2, maxReadLength = 101){
    outName <- paste0(strsplit(fastaGenome, split = "/") %>% unlist %>% tail(1),
                      ".STAR")
    if(is.null(outDir)){
        mkTmpDir()
        outDir <- file.path(getwd(), "STORMtmp_dir", outName)
    }else if(!is.null(outDir)){
        outDir <- outDir
    }
    genome <- txtools::tx_load_genome(fastaGenome)
    lGen <- sum(Biostrings::width(genome))
    nRef <- length(Biostrings::width(genome))
    if(nRef > 1000){
        genomChrBinNbits <- floor(min(18,log2(max(lGen/nRef, maxReadLength))))
    }else{
        genomChrBinNbits <- 18
    }
    genomeindexNb <- floor(min(14, log2(lGen)/2 - 1))
    if(is.null(bedAnnotation)){
        com <- paste("/apps/RH7U2/general/STAR/2.7.5c/bin/Linux_x86_64/STAR",
                     "--runMode genomeGenerate",
                     "--runThreadN", nCores,
                     "--genomeDir", outDir,
                     "--genomeFastaFiles", fastaGenome,
                     "--genomeChrBinNbits", genomChrBinNbits,
                     "--genomeSAindexNbases", genomeindexNb)
        system(com)
    }else{
        tmpF <- tempfile() %>% strsplit(split = "/") %>% unlist %>% tail(1)
        tmpGTF <- file.path(getwd(), "STORMtmp_dir", tmpF)
        mkGTF(bedAnnotation, tmpGTF)
        com <- paste("/apps/RH7U2/general/STAR/2.7.5c/bin/Linux_x86_64/STAR",
                     "--runMode genomeGenerate",
                     "--runThreadN", nCores,
                     "--genomeDir", outDir,
                     "--genomeFastaFiles", fastaGenome,
                     "--sjdbOverhang", maxReadLength - 1,
                     "--sjdbGTFfile", tmpGTF,
                     "--genomeChrBinNbits", genomChrBinNbits,
                     "--genomeSAindexNbases", genomeindexNb)
        system(com)
    }

    invisible(file.remove(tmpGTF))
    outDir
}

# Library complexity reports lce.txt
libComplexReport <- function(META, maxExtrapolation = 2.01e6, steps = 1e5, verbose = FALSE){
    if(all(file.exists(META$BAM))){
        for(file in META$BAM){
            com <- paste0("/apps/RH7U2/gnu/preseq/2.0.1/preseq lc_extrap -P -B ",
                          "-e ", maxExtrapolation, " -s ", steps, " ", file, " -o ",
                          gsub(pattern = ".bam", replacement = ".lce.txt", x = file),
                          " &")
            system(com)
        }
        Sys.sleep(time = 5)
        lastBam <- gsub(pattern = ".bam", replacement = ".lce.txt",x = META$BAM)
        while(min(difftime(Sys.time(), file.info(lastBam)$mtime, units = "secs")) < 60){
            Sys.sleep(time = 2)
        }
    }else{
        stop("Files ", paste(META$BAM[!file.exists(META$BAM)], collapse = " "),
             " do not exist")
    }
    if(verbose){cat("DONE: Library complexity reports.")}
}

