
# get-gene ----------------------------------------------------------------
#' @title Get gene data associated with a UniProt XML
#' @param uni_xml An xml_document associated with a uniprot data base,
#' for example, read_html(getURL("http://www.uniprot.org/uniprot/A0A015QDR1.xml"))
#' @importFrom rvest html_nodes html_attr html_attrs html_text
#' @importFrom magrittr %>%
#' @export
get_gene <- function(uni_xml) {
  gene_xml <- uni_xml %>% html_nodes("gene > name")
  data.frame(type = html_attr(gene_xml, "type"),
             evidence = html_attr(gene_xml, "evidence"),
             gene_id = html_text(gene_xml))
}

# get-protein -------------------------------------------------------------
#' @title Get protein data associated with a UniProt XML
#' @param uni_xml An xml_document associated with a uniprot data base,
#' for example, read_html(getURL("http://www.uniprot.org/uniprot/A0A015QDR1.xml"))
#' @importFrom rvest html_nodes html_attr html_attrs html_text
#' @importFrom magrittr %>%
#' @export
get_protein <- function(uni_xml) {
  protein_xml <- uni_xml %>%
    html_nodes("protein > submittedname > fullname")
  data.frame(evidence = html_attr(protein_xml, "evidence"),
             protein_name = html_text(protein_xml))
}

# get-organism ------------------------------------------------------------
#' @title Get organism data associated with a UniProt XML
#' @importFrom rvest html_nodes html_attr html_attrs html_text
#' @importFrom magrittr %>%
#' @param uni_xml An xml_document associated with a uniprot data base,
#' for example, read_html(getURL("http://www.uniprot.org/uniprot/A0A015QDR1.xml"))
#' @export
get_organism <- function(uni_xml) {
  evidence <- uni_xml %>%
    html_nodes("organism") %>%
    html_attr("evidence")
  scientific_name <- uni_xml %>%
    html_nodes("organism > name") %>%
    html_text
  ncbi_taxonomy <- uni_xml %>%
    html_nodes("organism > dbreference") %>%
    html_attr("id")
  data.frame(evidence, scientific_name, ncbi_taxonomy)
}

# get-authors -------------------------------------------------------------
#' @title Get authorship data associated with a UniProt XML
#' @importFrom rvest html_nodes html_attr html_attrs html_text
#' @importFrom magrittr %>%
#' @param uni_xml An xml_document associated with a uniprot data base,
#' for example, read_html(getURL("http://www.uniprot.org/uniprot/A0A015QDR1.xml"))
#' @export
get_authors <- function(uni_xml) {
  uni_xml %>%
    html_nodes("authorlist > person") %>%
    html_attrs %>% unlist %>% unname
}

# get-citation ------------------------------------------------------------
#' @title Get citation data associated with a UniProt XML
#' @param uni_xml An xml_document associated with a uniprot data base,
#' for example, read_html(getURL("http://www.uniprot.org/uniprot/A0A015QDR1.xml"))
#' @importFrom magrittr %>%
#' @importFrom rvest html_nodes html_attr html_attrs html_text
#' @export
get_citation <- function(uni_xml) {
  uni_xml %>%
    html_nodes("citation") %>%
    html_attrs %>% unlist
}

# get-comment -------------------------------------------------------------
#' @title Get comment data associated with a UniProt XML
#' @importFrom rvest html_nodes html_attr html_attrs html_text
#' @param uni_xml An xml_document associated with a uniprot data base,
#' for example, read_html(getURL("http://www.uniprot.org/uniprot/A0A015QDR1.xml"))
#' @importFrom magrittr %>%
#' @export
get_comment <- function(uni_xml) {
  uni_xml %>%
    html_nodes("comment > text") %>%
    html_text
}

# get-source --------------------------------------------------------------
#' @title Get source data associated with a UniProt XML
#' @param uni_xml An xml_document associated with a uniprot data base,
#' for example, read_html(getURL("http://www.uniprot.org/uniprot/A0A015QDR1.xml"))
#' @importFrom magrittr %>%
#' @importFrom rvest html_nodes html_attr html_attrs html_text
#' @export
get_source <- function(uni_xml) {
  uni_xml %>%
    html_nodes("source > dbreference") %>%
    html_attrs %>% unlist
}

# get-sequence ------------------------------------------------------------
#' @title Get sequence data associated with a UniProt XML
#' @param uni_xml An xml_document associated with a uniprot data base,
#' for example, read_html(getURL("http://www.uniprot.org/uniprot/A0A015QDR1.xml"))
#' @importFrom rvest html_nodes html_attr html_attrs html_text
#' @importFrom magrittr %>%
#' @export
get_sequence <- function(uni_xml) {
  sequence <- uni_xml %>%
    html_nodes("sequence") %>%
    html_text
  gsub("\n", "\\\\n", sequence) # escape the \n's
}
