library(tidyverse)

load("data-raw/runDb.RData")
DB_ptk <- DB %>% ungroup()

load("data-raw/runDb_stk.RData")
DB_stk <- DB %>% ungroup()

uka_db_full <- rbind(DB_ptk,DB_stk)
saveRDS(uka_db_full, "data/uka_map_db.RDS")


# kinome mapping file
kinome_mp_file_v1 <- read_delim("data-raw/KinomeMapping/kinase_mapping.txt", delim = "\t",
                      na = c("N/A", "Not found by Ethan", "Not found by Jake", "Not found",
                             "Not Found By Ethan", "Not Found By Jake"
                      )
)

kinome_mp_file_v2 <- read_delim("data-raw/KinomeMapping/2021_05_20-creedenzymatic_map.txt", delim = "\t",
                                na = c("N/A","Not found by Ethan", "Not found by Jake", "Not found",
                                       "Not Found By Ethan", "Not Found By Jake")) %>%
  filter(!is.na(hgnc_symbol)) %>%
  mutate_at(c("class", "group", "family", "subfamily", "krsa_id", "uka_id", "kea3_id", "ptmsea_id"), toupper) %>%
  select(1:12)

kinome_mp_file_v3 <- read_delim("data-raw/KinomeMapping/2021_06_11-creedenzymatic_map.txt", delim = "\t") %>%
  mutate_at(c("group", "family", "subfamily", "krsa_id", "uka_id", "kea3_id", "ptmsea_id"), toupper) %>%
  filter(!is.na(hgnc_symbol)) %>%
  mutate(subfamily = ifelse(subfamily == "N/A", family, subfamily)) %>%
  select(1:26)


kinome_mp_file_v4 <- read_delim("data-raw/KinomeMapping/2021_07_09-creedenzymatic_map.txt", delim = "\t",
                                na = c("N/A")
                                ) %>%
  mutate_at(c("group", "family", "subfamily", "krsa_id", "uka_id", "kea3_id", "ptmsea_id"), toupper) %>%
  filter(!is.na(hgnc_symbol)) %>%
  mutate(subfamily = ifelse(subfamily == "N/A", family, subfamily)) %>%
  select(1:26)


kinome_mp_file <- kinome_mp_file_v4


# peptide to HGNC

stk_pamchip_87102_mapping <- read_delim("data-raw/Pamchips_Layout/2021_10_13-JFC_complete-stk_peptides_mapping.txt", delim = "\t") %>%
  select(ID, HGNC)
ptk_pamchip_86402_mapping <- read_delim("data-raw/Pamchips_Layout/2021_10_13-JFC_complete-ptk_peptides_mapping.txt", delim = "\t") %>%
  select(ID, HGNC)

stk_pamchip_87102_array_layout <- readxl::read_xlsx("data-raw/Pamchips_Layout/87102_ArrayLayout.xlsx") %>%
  filter(ID != "#REF", Ser != "NA", Thr != "NA")  %>% select(ID, Ser, Thr, UniprotAccession) %>%
  mutate(Ser = gsub("\\[", "", Ser), Ser = gsub("\\]", "", Ser),Ser = gsub("\\s+", "", Ser),
         Thr = gsub("\\[", "", Thr), Thr = gsub("\\]", "", Thr),Thr = gsub("\\s+", "", Thr),
         UniprotAccession = gsub("†", "", UniprotAccession),
         ID = gsub("\\s+", "", ID)
         )

ptk_pamchip_86402_array_layout <- readxl::read_xlsx("data-raw/Pamchips_Layout/86402_ArrayLayout.xlsx") %>%
  filter(ID != "#REF", Tyr != "NA") %>%
  select(ID, Tyr, UniprotAccession) %>%
  mutate(Tyr = gsub("\\[", "", Tyr), Tyr = gsub("\\]", "", Tyr), Tyr = gsub("\\s+", "", Tyr)) %>%
  mutate(UniprotAccession = gsub("†", "", UniprotAccession),
         ID = gsub("\\s+", "", ID))



# format O60934;S343-p
ptk_pamchip_86402_array_layout %>%
  separate_rows(Tyr, sep = ",") %>%
  mutate(Tyr = paste0("Y", Tyr),
         PTM_SEA_ID = paste0(UniprotAccession, ";", Tyr, "-p")) -> ptk_pamchip_86402_array_layout_ptmsea


stk_pamchip_87102_array_layout %>%
  separate_rows(Ser, sep = ",") %>%
  separate_rows(Thr, sep = ",") %>%
  mutate(Ser = paste0("S", Ser),
         Thr = paste0("T", Thr)) %>%
  pivot_longer(2:3,names_to = "AA", values_to = "Position") %>%
  filter(nchar(Position) != 1) %>%
  select(-AA) %>%
  mutate(PTM_SEA_ID = paste0(UniprotAccession, ";", Position, "-p")) -> stk_pamchip_87102_array_layout_ptmsea

# PTM DB
ptmsea_all <- cmapR::parse_gmt("data-raw/PTM_SEA_DB/ptm.sig.db.all.uniprot.human.v1.9.0.gmt")

covertList <- function(x) {
  tibble(entry = x$entry,
         head = x$head
  )

}

processedDB <- map_df(ptmsea_all, covertList)

processedDB %>% filter(grepl("Kinase",head, ignore.case = T)) %>%
  filter(grepl(";Y", entry)) %>%
  group_by(head) %>% mutate(len = n()) %>%
  ungroup() -> filterDBS_ptk


filterDBS_ptk %>% mutate(ids = entry, ids = str_remove(ids, ";u")) %>%
  filter(ids %in% ptk_pamchip_86402_array_layout_ptmsea$PTM_SEA_ID) %>%
  select(head, entry) %>%
  group_by(head) %>% mutate(len = n()) %>%
  ungroup() %>%
  group_by(head) %>%
  mutate(entry = list(unique(entry))) %>% ungroup() %>%
  distinct() %>%
  select(head, entry, len) -> onlyChipPeps_dbs

onlyChipPeps_dbs <- setNames(as.list(as.data.frame(t(onlyChipPeps_dbs))), onlyChipPeps_dbs$head)

cmapR::write_gmt(onlyChipPeps_dbs, "data-raw/PTM_SEA_DB/ptk_pamchip_86402_onlyChipPeps_dbs.gmt")


processedDB %>% filter(grepl("Kinase",head, ignore.case = T)) %>%
  filter(grepl(";Y|S", entry)) %>%
  group_by(head) %>% mutate(len = n()) %>%
  ungroup() -> filterDBS_stk


filterDBS_stk %>% mutate(ids = entry, ids = str_remove(ids, ";u")) %>%
  filter(ids %in% stk_pamchip_87102_array_layout_ptmsea$PTM_SEA_ID) %>%
  select(head, entry) %>%
  group_by(head) %>% mutate(len = n()) %>%
  ungroup() %>%
  group_by(head) %>%
  mutate(entry = list(unique(entry))) %>% ungroup() %>%
  distinct() %>%
  select(head, entry, len) -> onlyChipPeps_dbs_stk

onlyChipPeps_dbs_stk <- setNames(as.list(as.data.frame(t(onlyChipPeps_dbs_stk))), onlyChipPeps_dbs_stk$head)

cmapR::write_gmt(onlyChipPeps_dbs_stk, "data-raw/PTM_SEA_DB/stk_pamchip_87102_onlyChipPeps_dbs.gmt")



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
                  kinome_mp_file,
                  kinome_mp_file_v1,
                  kinome_mp_file_v2,
                  kinome_mp_file_v3,
                  kinome_mp_file_v4,
                  stk_pamchip_87102_mapping,
                  stk_pamchip_87102_array_layout_ptmsea,
                  ptk_pamchip_86402_mapping,
                  ptk_pamchip_86402_array_layout_ptmsea,
                  overwrite = T)


