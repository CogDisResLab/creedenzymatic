#' Run PTM-SEA API using a gct file as input
#'
#' This function takes in a gct file (created by the read_prmsea function) and run PTM-SEA API and returns results
#'
#' @param gene_set vector, HGNC gene symbols based on the differentially phosphorylated peptides
#' @param lib searched kea libraries "iptmnet" or "ptm-sea" or "all (default is "iptmnet" which uses the iptmnet mapping)
#' @param nperm number of permutations
#' @param min.overlap minimum overlap of target peptides with referernce peptides sets
#' @param ... additional arguments passed to the ssGSEA_ce function
#'
#' @return list
#'
#' @export
#'
#'
#'

run_ptmsea <- function(gct_object, lib = "iptmnet", nperm = 1000, min.overlap = 1, ...) {


  message("running ptm-sea ...")
  res <- ssGSEA2_ce(input.ds = gct_object,
                    output.prefix = "test_ce",
                    gene.set.databases =
                      if(lib == "iptmnet") {c(list(ptm_sea_iptmnet_mapping_stk), list(ptm_sea_iptmnet_mapping_ptk))}
                      else if(lib == "ptm-sea") {
                      c(list(ptk_pamchip_86402_onlyChipPeps_dbs), list(stk_pamchip_87102_onlyChipPeps_dbs))}
                    else {

                      c(list(ptmsea_all_dbs))
                    }
                    ,
                    nperm = nperm,
                    export.signat.gct=F,
                    sample.norm.type = 	"rank",
                    weight = 	0.75,
                    statistic =	"area.under.RES",
                    output.score.type = 	"NES",
                    min.overlap = min.overlap,
                    correl.type =	"z.score", ...
  )

  res


}
