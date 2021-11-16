library(cmapR)
library(dplyr)
library(readr)

#source("data-raw/PTM_SEA_DB/gct-io.R")
source("data-raw/PTM_SEA_DB/ssGSEA2.0.R")
source("data-raw/PTM_SEA_DB/ssGSEA2.0_ce_version.R")

res <- ssGSEA2(
  input.ds=opt$input.ds,
  output.prefix=opt$output.prefix,
  gene.set.databases=opt$gene.set.databases,
  sample.norm.type=opt$sample.norm.type,
  weight=opt$weight,
  statistic=opt$statistic,
  output.score.type=opt$output.score.type,
  nperm=opt$nperm,
  min.overlap=opt$min.overlap,
  correl.type=opt$correl.type,
  export.signat.gct=opt$export.signat.gct,
  extended.output=opt$extended.output,
  #	fdr.pvalue=opt$fdr.pvalue,
  global.fdr=opt$global.fdr,
  par=opt$par,
  spare.cores=spare.cores,
  log.file=log.file
)

ex <- cmapR::parse_gctx("data-raw/PTM_SEA_DB/PI3K_pert_logP_n2x23936.gct")
ptmsea_all <- cmapR::parse_gmt("data-raw/PTM_SEA_DB/ptm.sig.db.all.flanking.human.v1.8.1.gmt.txt")

res <- ssGSEA2(input.ds = "data-raw/PTM_SEA_DB/PI3K_pert_logP_n2x23936.gct",
        output.prefix = "test",
        gene.set.databases = "data-raw/PTM_SEA_DB/ptm.sig.db.all.flanking.human.v1.8.1.gmt.txt",
        nperm = 5,
        export.signat.gct=F
        )

input_ds <- read_delim("data-raw/PTM_SEA_DB/PDCL15vsWT_input_table.txt", delim = " ")

input_ds %>% select(Peptide, totalGeoMeanLFC) %>% distinct() -> input_ds_mod

input_ds_mod %>% left_join(ptk_pamchip_86402_array_layout_ptmsea, by = c("Peptide" = "ID")) %>%
  select(PTM_SEA_ID, totalGeoMeanLFC) %>% distinct(PTM_SEA_ID, .keep_all = T) -> input_ds_mod_filtered

input_ds_mod_filtered %>%
  column_to_rownames("PTM_SEA_ID") %>%
  as.matrix() -> mat_ds

(my_new_ds <- new("GCT", mat=mat_ds))

cmapR::write_gct(my_new_ds,"data-raw/PTM_SEA_DB/PDCL15vsWT_input_PTM-SEA.gct")


res <- ssGSEA2(input.ds = "data-raw/PTM_SEA_DB/PDCL15vsWT_input_PTM-SEA.gct_n1x167.gct",
               output.prefix = "test",
               gene.set.databases = "data-raw/PTM_SEA_DB/ptk_pamchip_86402_onlyChipPeps_dbs.gmt",
               nperm = 1000,
               export.signat.gct=F,
               sample.norm.type = 	"rank",
               #global.fdr = T,
               weight = 	0.75,
               statistic =	"area.under.RES",
               output.score.type = 	"NES",
               min.overlap = 1,
               correl.type =	"z.score"
)

res <- ssGSEA2_ce(input.ds = my_new_ds,
               output.prefix = "test_ce",
               gene.set.databases = c(#"data-raw/PTM_SEA_DB/stk_pamchip_87102_onlyChipPeps_dbs.gmt",
                                      "data-raw/PTM_SEA_DB/ptk_pamchip_86402_onlyChipPeps_dbs.gmt"),
               nperm = 1000,
               export.signat.gct=F,
               sample.norm.type = 	"rank",
               #global.fdr = T,
               weight = 	0.75,
               statistic =	"area.under.RES",
               output.score.type = 	"NES",
               min.overlap = 1,
               correl.type =	"z.score"
)




res %>%  as.data.frame() -> output_res

res %>% lapply(`[[`, 1) %>%
  lapply(`[[`, 1) %>% as.data.frame() %>%
  pivot_longer(cols = 1:ncol(.), names_to = "Set", values_to = "ES") %>%
  mutate(Kinase = str_extract(Set, "_.+"),
         Kinase = str_remove(Kinase, "_"),
         Kinase = str_replace(Kinase, "\\.", "/") %>% toupper()
  ) %>% dplyr::select(Kinase, Score = ES)

data_frame_res <- as.data.frame(do.call(rbind, res)) %>%
  rownames_to_column("Set") %>%
  unnest(V1) %>% unnest(V1)

data.table::rbindlist(res, fill=F, use.names = F) -> data_frame_res

head(output_res)


output_res %>% mutate(Kinase = str_extract(id, "_.+"),
                      Kinase = str_remove(Kinase, "_")
) -> output_res


output_res %>% select(id, fdr.pvalue.totalGeoMeanLFC) %>%
  ggplot(aes(reorder(id, -log10(fdr.pvalue.totalGeoMeanLFC)), -log10(fdr.pvalue.totalGeoMeanLFC))) +
  geom_col() + coord_flip()

output_res %>% select(Kinase, fdr.pvalue.totalGeoMeanLFC) %>%
  ggplot(aes(reorder(Kinase, -log10(fdr.pvalue.totalGeoMeanLFC)), -log10(fdr.pvalue.totalGeoMeanLFC))) +
  geom_col() + coord_flip() +
  theme_bw() +
  labs(x = "", y = "-log10 adj P Value")


ex <- cmapR::parse_gctx("test_ce-combined.gct")
ex <- cmapR::parse_gctx("test_ce-scores.gct")

ptmsea_all <- cmapR::parse_gmt("data-raw/PTM_SEA_DB/stk_pamchip_87102_onlyChipPeps_dbs.gmt")
readLines(ptmsea_all)

rl_ex <- readLines("data-raw/PTM_SEA_DB/stk_pamchip_87102_onlyChipPeps_dbs.gmt")
rl_ex

for (x in c(list(ptk_pamchip_86402_onlyChipPeps_dbs), list(stk_pamchip_87102_onlyChipPeps_dbs))) {
  print(x)
}
