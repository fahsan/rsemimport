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


#Declare function and inputs
rsemimport <- function(files, importer = NULL) {
  
  #For loop through index of samples in the files supplied. 
  for (sample in seq_along(files)) {
    
    #Output sample
    message(sample," ", appendLF=FALSE)
    
    #Capture data from files provided in
    out <- capture.output({
      files.to.import <- as.data.frame(importer(files[sample]))
    }, type = "message")
    
    #Add titles to all data in files to import
    colnames(files.to.import) <- c("Isoform", "Chromosome", "Strand", "txStart", "txEnd", "cdsStart", 
                                   "cdsEnd", "exonCount", "exonStarts", "exonEnds", "Hugo", "alignID", "Canonical",
                                   "length", "effective_length", "tx_expected_count", "txTPM", "txFPKM", "IsoPct",
                                   "gene_expected_count", "geneTPM", "geneFPKM")
    #Extract gene exp counts, TPM, and FPKM abundance values. 
    geneIdCol <- "Hugo"
    expcountsCol <- "gene_expected_count"
    TPMabundanceCol <- "geneTPM"
    FPKMabundanceCol <- "geneFPKM"
    unique(files.to.import)
    
    #Kill command if all params are not set. 
    stopifnot(all(c(geneIdCol, TPMabundanceCol, FPKMabundanceCol, expcountsCol) %in% names(files.to.import)))
    
    #Iterate through each sample and prepare matrix columns. 
    if (sample == 1) {
      gene.matrix <- matrix(nrow=nrow(files.to.import), ncol=length(files))
      rownames(gene.matrix) <- files.to.import[[geneIdCol]]
      colnames(gene.matrix) <- names(files)
      
      tpmabundance.matrix <- gene.matrix
      expcounts.matrix <- gene.matrix
      fpkmabundance.matrix <- gene.matrix
     }
     
    #Combine all matrix columns together.  
    tpmabundance.matrix[,sample] <- files.to.import[[TPMabundanceCol]]
     expcounts.matrix[,sample] <- files.to.import[[expcountsCol]]
     fpkmabundance.matrix[,sample] <- files.to.import[[FPKMabundanceCol]]
    }
   
  #@return the list. 
   message("")
   return(list(TPM = tpmabundance.matrix, FPKM = FPKMabundance.matrix, expcounts = expcounts.matrix, countsFromAbundance = "no")
}
