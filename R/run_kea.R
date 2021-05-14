#' Run KEA3 API based on a set of gene symbols
#'
#' This function takes in HGNC gene symbols and connect to KEA3 API and returns results
#'
#' @param gene_set vector, HGNC gene symbols based on the differentially phosphorylated peptides
#' @param lib searched kea libraries (default is "kinases" which will return only kinase libraries like ChengKSIN, PTMsigDB, PhosDAll)
#'
#' @return list, tables from each KEA3 library
#'
#' @export
#'
#'
#'

run_kea <- function(gene_set, lib = c("kinases")) {

  rbind(stk_pamchip_87102_mapping, ptk_pamchip_86402_mapping) %>%
    dplyr::filter(ID %in% gene_set)

  url <- "https://amp.pharm.mssm.edu/kea3/api/enrich/"

  response <-  httr::POST(url = url, body = list(gene_set = gene_set), encode = "json")

  if (response$status_code == 200) {

    json <-  httr::content(response, "text")
    results <-  jsonlite::fromJSON(json)

    if(lib == "kinases") {
      results[names(results) %in% c("PTMsigDB", "ChengKSIN", "PhosDAll")]
    }

    else {
      results
    }

    # body_res$data[[1]] %>% enframe() %>% unnest(cols = c("value")) %>%
    #   pivot_wider(names_from = name, values_from = value)

  }
  else {message("Couldn't connect to KEA3 API successfully")}


}
