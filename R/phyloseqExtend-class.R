
# phyloseqExtend-class ----------------------------------------------------
#'  @title Assemble several data objects into a phyloseqExtend object.
#'
#'  @description This creates a phyloseqExtend object from the multiple sources
#'    of data common in microbiome experiments.
#'
#'  @param otu_table_object An object of class otu_table, to incorporate into
#'    the phyloseqExtend object. Defaults to NULL.
#'  @param sample_data_object An object of class sam_data, to incorporate into
#'    the phyloseqExtend object. Defaults to NULL.
#'  @param phy_tree_object An object of class phylo, to incorporate into the
#'    phyloseqExtend object. Defaults to NULL.
#'  @param spectra_object An object of class spectra, to incorporate into the
#'    phyloseqExtend object. Defaults to NULL.
#'  @param refseq_object An object of class refseq, to incorporate into the
#'    phyloseqExtend object. Defaults to NULL.
#'
#'  @return phyloseqExtend_object An object of class phyloseqExtend, with slots
#'    containing each of the input data objects, and NULL for potential but
#'    unoccupied slots.
phyloseqExtend <- function(otu_table_object = NULL,
                           sample_data_object = NULL,
                           phy_tree_object = NULL,
                           refseq_object = NULL,
                           spectra_object = NULL) {
  phyloseqExtend_object <- new("phyloseqExtend",
                               otu_table = otu_table_object,
                               sam_data = sample_data_object,
                               phy_tree = phy_tree_object,
                               refseq = refseq_object,
                               spectra = spectra_object)
  return (phyloseqExtend_object)
}
