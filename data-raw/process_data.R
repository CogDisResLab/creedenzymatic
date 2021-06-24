library(readr)
library(dplyr)

load("data-raw/runDb.RData")
DB_ptk <- DB %>% ungroup()

load("data-raw/runDb_stk.RData")
DB_stk <- DB %>% ungroup()

uka_db_full <- rbind(DB_ptk,DB_stk)
saveRDS(uka_db_full, "data/uka_map_db.RDS")


# kinome mapping file
kinome_mp_file_v1 <- read_delim("data-raw/kinase_mapping.txt", delim = "\t",
                      na = c("N/A", "Not found by Ethan", "Not found by Jake", "Not found",
                             "Not Found By Ethan", "Not Found By Jake"
                      )
)

kinome_mp_file_v2 <- read_delim("data-raw/2021_05_20-creedenzymatic_map.txt", delim = "\t",
                                na = c("N/A","Not found by Ethan", "Not found by Jake", "Not found",
                                       "Not Found By Ethan", "Not Found By Jake")) %>%
  filter(!is.na(hgnc_symbol)) %>%
  mutate_at(c("class", "group", "family", "subfamily", "krsa_id", "uka_id", "kea3_id", "ptmsea_id"), toupper) %>%
  select(1:12)

kinome_mp_file_v3 <- read_delim("data-raw/2021_06_11-creedenzymatic_map.txt", delim = "\t") %>%
  mutate_at(c("group", "family", "subfamily", "krsa_id", "uka_id", "kea3_id", "ptmsea_id"), toupper) %>%
  filter(!is.na(hgnc_symbol)) %>%
  mutate(subfamily = ifelse(subfamily == "N/A", family, subfamily)) %>%
  select(1:26)


# peptide to HGNC

stk_pamchip_87102_mapping <- read_delim("data-raw/2021_10_13-JFC_complete-stk_peptides_mapping.txt", delim = "\t") %>%
  select(ID, HGNC)
ptk_pamchip_86402_mapping <- read_delim("data-raw/2021_10_13-JFC_complete-ptk_peptides_mapping.txt", delim = "\t") %>%
  select(ID, HGNC)

KRSA::KRSA_coverage_STK_PamChip_87102_v2$Substrates %>% unique() -> stk_peps
KRSA::KRSA_coverage_PTK_PamChip_86402_v1$Substrates %>% unique() %>% as.character() -> ptk_peps

setdiff(stk_peps, stk_pamchip_87102_mapping$ID)

add_stk <- tibble::tibble(
  ID = c("H2B1B_ 27_40","E1A_ADE05_212_224"),
  HGNC = c("H2BC3", NA)
)

setdiff(ptk_peps, ptk_pamchip_86402_mapping$ID)

add_ptk <- tibble::tibble(
  ID = c("CD3Z_147_159"),
  HGNC = c("CD247")
)


rbind(stk_pamchip_87102_mapping, add_stk) -> stk_pamchip_87102_mapping
rbind(ptk_pamchip_86402_mapping, add_ptk) -> ptk_pamchip_86402_mapping

usethis::use_data(uka_db_full,
                  kinome_mp_file_v1,
                  kinome_mp_file_v2,
                  kinome_mp_file_v3,
                  stk_pamchip_87102_mapping,
                  ptk_pamchip_86402_mapping,
                  overwrite = T)

