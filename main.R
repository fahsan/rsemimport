#main.R
#Fasih Ahsan, Teitell Lab UCLA Pathology and Laboratory Medicine
#9.29.17
#Takes as input several RSEM gene.results files that are specially edited, and creates a list of gene count/TPM/FPKM...
#... matrices for use in DESeq2, PCA, GSEA, and other methods that require normalized or count data for RNA-seq experiments. 
#@param files type character vector input describing location and associated metadata with the RNA-seq RSEM files of interest.
#@param importer function used to import in files. 
#@return R type list of matrices, containing gene count, TPM abudance, and FPKM abudance matrices. 
#@references (This code is adapted wholly from):
#' Charlotte Soneson, Michael I. Love, Mark D. Robinson (2015):
#' Differential analyses for RNA-seq: transcript-level estimates
#' improve gene-level inferences. F1000Research.
#' \url{http://dx.doi.org/10.12688/f1000research.7563.1}
#@dependencies utils read.delim (if readr not installed) capture.output (R CRAN)

#Declare function and inputs
rsemimport <- function(files, importer = NULL) {
  
  #Set read.delimiter function using importer. 
  if(is.null(importer)) {
    # If the readr utility is not installed, feed importer function read.delim to enter package data. 
    if(!requireNamespace("readr", quietly = TRUE)) {
      message("reading in files with read.delim function (install 'readr' package to speed this process up)")
      importer <- read.delim
     }
     #If the readr utility is installed, using readr function read_tsv to import the column data from each tab file. 
       else {
       message("reading in files with read_tsv")
       readrStatus <- TRUE
       importer <- function(x) readr::read_tsv(x, progress = FALSE, 
                                               #Add column labels to each file for downstream extraction. 
                                              col_names=c("Isoform", "Chromosome", "Strand", "txStart", "txEnd", "cdsStart", 
                                                                               "cdsEnd", "exonCount", "exonStarts", "exonEnds", "Hugo", "alignID", "Canonical",
                                                                               "length", "effective_length", "tx_expected_count", "txTPM", "txFPKM", "IsoPct",
                                                                               "gene_expected_count", "geneTPM", "geneFPKM"), 
                                              col_types=readr::cols())
     }
   }
  
  #For loop through index of samples in the files supplied. 
  for (sample in seq_along(files)) {
    
    #Output sample
    message(sample," ", appendLF=FALSE)
    
    #Capture data from files provided in
    out <- capture.output({
      files.to.import <- as.data.frame(importer(files[sample]))
    }, type = "message")
    
    #Extract gene exp counts, TPM, and FPKM abundance values. 
    geneIdCol <- "Hugo"
    expcountsCol <- "gene_expected_count"
    TPMabundanceCol <- "geneTPM"
    FPKMabundanceCol <- "geneFPKM"
    effectivelengthsCol <- "effective_length"
    unique(files.to.import)
    
    #Kill command if all params are not set. 
    stopifnot(all(c(geneIdCol, TPMabundanceCol, FPKMabundanceCol, expcountsCol, effectivelengthsCol) %in% names(files.to.import)))
    
    #Iterate through each sample and prepare matrix columns. 
    if (sample == 1) {
      gene.matrix <- matrix(nrow=nrow(files.to.import), ncol=length(files))
      rownames(gene.matrix) <- files.to.import[[geneIdCol]]
      colnames(gene.matrix) <- names(files)
      
      tpmabundance.matrix <- gene.matrix
      expcounts.matrix <- gene.matrix
      fpkmabundance.matrix <- gene.matrix
      effective_length.matrix <- gene.matrix
     }
     
    #Combine all matrix columns together.  
    tpmabundance.matrix[,sample] <- files.to.import[[TPMabundanceCol]]
    expcounts.matrix[,sample] <- files.to.import[[expcountsCol]]
    fpkmabundance.matrix[,sample] <- files.to.import[[FPKMabundanceCol]]
    effective_length.matrix[,sample] <- files.to.import[[effectivelengthsCol]]
    
    }
  
  #Remove duplicate rows from the matrices (since gene-level counts were duplicated based on transcript ID). 
  tpmabundance.matrix <- tpmabundance.matrix[!duplicated(rownames(tpmabundance.matrix)),]
  fpkmabundance.matrix <- fpkmabundance.matrix[!duplicated(rownames(fpkmabundance.matrix)),]
  expcounts.matrix <- expcounts.matrix[!duplicated(rownames(expcounts.matrix)),]
  effective_length.matrix <- effective_length.matrix[!duplicated(rownames(effective_length.matrix)),]
  
  #Convert all matrix counts with 0 to 1 to prevent DESeq2 error. 
  effective_length.matrix[effective_length.matrix == 0] <- 1
  
  #@return the list. 
   message("")
   return(list(TPM = tpmabundance.matrix, FPKM = fpkmabundance.matrix, counts = expcounts.matrix, lengths = effective_length.matrix, countsFromAbundance = "no"))
}
#done
