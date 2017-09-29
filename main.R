#main.R
#Takes as input
rsemimport <- function(files, importer = NULL) {
  
  for (sample in seq_along(files)) {
    message(sample," ", appendLF=FALSE)
    
    out <- capture.output({
      files.to.import <- as.data.frame(importer(files[sample]))
    }, type = "message")
    
    colnames(files.to.import) <- c("Isoform", "Chromosome", "Strand", "txStart", "txEnd", "cdsStart", 
                                   "cdsEnd", "exonCount", "exonStarts", "exonEnds", "Hugo", "alignID", "Canonical",
                                   "length", "effective_length", "tx_expected_count", "txTPM", "txFPKM", "IsoPct",
                                   "gene_expected_count", "geneTPM", "geneFPKM")
    geneIdCol <- "Hugo"
    expcountsCol <- "gene_expected_count"
    TPMabundanceCol <- "geneTPM"
    FPKMabundanceCol <- "geneFPKM"
    unique(files.to.import)
                            
    stopifnot(all(c(geneIdCol, TPMabundanceCol, FPKMabundanceCol, expcountsCol) %in% names(files.to.import)))
    
    if (sample == 1) {
      gene.matrix <- matrix(nrow=nrow(files.to.import), ncol=length(files))
      rownames(gene.matrix) <- files.to.import[[geneIdCol]]
      colnames(gene.matrix) <- names(files)
      
      tpmabundance.matrix <- gene.matrix
      expcounts.matrix <- gene.matrix
      fpkmabundance.matrix <- gene.matrix
     }
     
     tpmabundance.matrix[,sample] <- files.to.import[[TPMabundanceCol]]
     expcounts.matrix[,sample] <- files.to.import[[expcountsCol]]
     fpkmabundance.matrix[,sample] <- files.to.import[[FPKMabundanceCol]]
    }
   message("")
   return(list(TPM = tpmabundance.matrix, FPKM = FPKMabundance.matrix, expcounts = expcounts.matrix, countsFromAbundance = "no")
}
