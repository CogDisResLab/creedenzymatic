#' Run KEA3 based on  Rank UKA table
#'
#' reads UKA table and checks for correct format
#'
#' This function takes in UKA table and rank and quartile kinases based on the absolute Score values
#'
#' @param df dataframe, UKA table output (requires at least Kinase and Z columns)
#' @param ..., arguments passed to rank_kinases function
#'
#' @return dataframe, Ranked and quartiled UKA table
#'
#' @export
#'
#'
#'

# get_kea3_data <- function(filename) {
#   prefix <- "kea3_list"
#   input_path <- file.path(prefix, filename)
#   proteins <- read.table(input_path, stringsAsFactors = F)[,1]
#   url <- "https://amp.pharm.mssm.edu/kea3/api/enrich/"
#   encode <- "json"
#   payload <-  list(gene_set = proteins)
#
#   response <-  POST(url = url, body = payload, encode = encode)
#   json <-  content(response, "text")
#
#   results <-  fromJSON(json)
#   result_names <- names(results)
#   write.csv(results[[4]], file.path(prefix, "results", paste(filename, result_names[4], "csv", sep = ".")))
#   write.csv(results[[9]], file.path(prefix, "results",  paste(filename, result_names[9], "csv", sep = ".")))
#   write.csv(results[[11]], file.path(prefix, "results",  paste(filename, result_names[11], "csv", sep = ".")))
# }
#
# lapply(files, get_kea3_data)
