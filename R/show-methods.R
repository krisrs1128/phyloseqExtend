# show-phyloseqExtend -----------------------------------------------------
#' @title Print a phyloseqExtend object
#'
#' @export
#' @docType methods
setMethod("show", "phyloseqExtend", function(object){
  cat("phyloseqExtend-class experiment-level object", fill=TRUE)
  # print otu_table (always there).
  cat(paste("otu_table()   OTU Table:         [ ", ntaxa(otu_table(object)), " taxa and ",
            nsamples(otu_table(object)), " samples ]", sep = ""), fill = TRUE)
  # print Sample Data if there
  if(!is.null(sample_data(object, FALSE))){
    cat(paste("sample_data() Sample Data:       [ ", dim(sample_data(object))[1], " samples by ",
              dim(sample_data(object))[2],
              " sample variables ]", sep = ""), fill = TRUE)
  }

  # print tax Tab if there
  if(!is.null(tax_table(object, FALSE))){
    cat(paste("tax_table()   Taxonomy Table:    [ ", dim(tax_table(object))[1], " taxa by ",
              dim(tax_table(object))[2],
              " taxonomic ranks ]", sep = ""), fill = TRUE)
  }

  # print tree if there
  if(!is.null(phy_tree(object, FALSE))){
    cat(paste("phy_tree()    Phylogenetic Tree: [ ", ntaxa(phy_tree(object)), " tips and ",
              phy_tree(object)$Nnode,
              " internal nodes ]", sep = ""),
        fill = TRUE
    )
  }

  # print refseq summary if there
  if(!is.null(refseq(object, FALSE))){
    cat(paste("refseq()      ", class(refseq(object))[1], ":      [ ", ntaxa(refseq(object)), " reference sequences ]", sep = ""), fill=TRUE)
  }

  # print spectra summary if there
  spec_obj <- spectra(object)
  if(!is.null(spec_obj)) {
    cat(paste("spectra()     ", class(spec_obj)[1], ":           [ ", nrow(spec_obj), " samples by ", ncol(spec_obj), " intensities ]", sep = ""), fill = TRUE)
  }

})
