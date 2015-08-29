
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


# get-db ------------------------------------------------------------------
#' @title Get database data associated with a UniProt XML
#' @param uni_xml An xml_document associated with a uniprot data base,
#' for example, read_html(getURL("http://www.uniprot.org/uniprot/A0A015QDR1.xml")))
#' @importFrom rvest html_nodes html_attr html_attrs html_text
#' @importFrom magrittr %>%
#' @importFrom plyr rbind.fill
#' @export
get_db_ref <- function(uni_xml) {
  db_ref <- uni_xml %>%
    html_nodes("entry > dbreference")
  dbs <- html_attr(db_ref, "type")

  n_dbs <- length(unique(dbs))
  dbs_list <- vector(length = n_dbs, mode = "list")
  names(dbs_list) <- unique(dbs)

  # Loop over available databases
  for(i in seq_len(n_dbs)) {
    cur_db_xml <- db_ref[dbs == unique(dbs)[i]]
    cur_db_list <- list()

    # extract every element from the current data base
    for(db_elem in seq_len(length(cur_db_xml))) {
      db_elem_data <- cur_db_xml[db_elem] %>% html_nodes("property")
      db_elem_vec <- db_elem_data %>% html_attr("value")
      names(db_elem_vec) <- db_elem_data %>% html_attr("type")
      db_elem_vec["id"] <- html_attr(cur_db_xml[[db_elem]], "id")
      cur_db_list[[db_elem]] <- data.frame(t(db_elem_vec))
    }
    dbs_list[[i]] <- rbind.fill(cur_db_list)
  }
  return (dbs_list)
}

# get-uniprot -------------------------------------------------------------
#' @title Extract information from UniProt Accession ID
#' @param id The UniProt ID to extract information from.
#' @return A list containing information about the specified UniProt ID.
#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes html_attr html_attrs html_text
#' @importFrom RCurl getURL
#' @importFrom magrittr %>%
#' @examples
#' # compare this with http://www.uniprot.org/uniprot/A0A015QDR1
#' get_uniprot_info("A0A015QDR1")
#' @export
get_uniprot_info <- function(id) {
  uni <- "http://www.uniprot.org/uniprot"
  uni_xml <- read_html(getURL(paste0(uni, "/", id, ".xml")))

  result <- list()
  result$uniprot_id <- id
  result$lineage <- uni_xml %>% html_nodes("taxon") %>% html_text
  result$gene <- get_gene(uni_xml)
  result$protein <- get_protein(uni_xml)
  result$organism <- get_organism(uni_xml)
  result$authors <- get_authors(uni_xml)
  result$citation <- get_citation(uni_xml)
  result$comment <- get_comment(uni_xml)
  result$db_source <- get_source(uni_xml)
  result$ref_dbs <- get_db_ref(uni_xml)
  return (result)
}

