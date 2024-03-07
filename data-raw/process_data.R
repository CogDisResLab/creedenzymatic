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

kinome_mp_file_v5 <- read_delim("data-raw/KinomeMapping/2021_07_09-creedenzymatic_map_KA_edits.txt", delim = "\t",
                                na = c("N/A")) %>%
  mutate_at(c("group", "family", "subfamily", "krsa_id", "uka_id", "kea3_id", "ptmsea_id"), toupper) %>%
  filter(!is.na(hgnc_symbol)) %>%
  mutate(subfamily = ifelse(subfamily == "N/A", family, subfamily)) %>%
  mutate(kea3_id = ifelse(kea3_id == "NOT FOUND" | is.na(kea3_id), hgnc_symbol, kea3_id))

kinome_mp_file_v6 <- read_delim("data-raw/KinomeMapping/2024_03_07-creedenzymatic_map_AH-ASI_edits.txt", delim = "\t",
                                na = c("N/A")) %>%
  mutate_at(c("group", "family", "subfamily", "krsa_id", "uka_id", "kea3_id", "ptmsea_id"), toupper) %>%
  filter(!is.na(hgnc_symbol)) %>%
  mutate(subfamily = ifelse(is.na(subfamily), family, subfamily)) %>%
  mutate(kea3_id = ifelse(kea3_id == "NOT FOUND" | is.na(kea3_id), hgnc_symbol, kea3_id)) %>%
  filter(!is.na(hgnc_id))


kinome_mp_file <- kinome_mp_file_v6


# peptide to HGNC

stk_pamchip_87102_mapping <- read_delim("data-raw/Pamchips_Layout/2021_10_13-JFC_complete-stk_peptides_mapping-KA_edit.txt", delim = "\t") %>%
  select(ID, HGNC)

ptk_pamchip_86402_mapping <- read_delim("data-raw/Pamchips_Layout/2021_10_13-JFC_complete-ptk_peptides_mapping.txt", delim = "\t") %>%
  select(ID, HGNC)

#setdiff(stk_peps, stk_pamchip_87102_mapping$ID)

add_stk <- tibble::tibble(
  ID = c("H2B1B_ 27_40","E1A_ADE05_212_224"),
  HGNC = c("H2BC3", NA)
)

#setdiff(ptk_peps, ptk_pamchip_86402_mapping$ID)

add_ptk <- tibble::tibble(
  ID = c("CD3Z_147_159"),
  HGNC = c("CD247")
)


rbind(stk_pamchip_87102_mapping, add_stk) -> stk_pamchip_87102_mapping
rbind(ptk_pamchip_86402_mapping, add_ptk) -> ptk_pamchip_86402_mapping


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

ptk_pamchip_86402_array_layout <- readxl::read_xlsx("data-raw/Pamchips_Layout/86402_ArrayLayout.xlsx") %>%
  filter(ID != "#REF", Tyr != "NA") %>%
  select(ID, Tyr, UniprotAccession) %>%
  mutate(
    UniprotAccession = gsub("†", "", UniprotAccession),
    Tyr = str_replace(Tyr, "\\[", ""),
    Tyr = str_replace(Tyr, "\\]", ""),
    Tyr = gsub("\\s+", "", Tyr),
    ID = gsub("\\s+", "", ID)) %>%
  separate_rows(Tyr, sep = ",") %>%
  mutate(site_residue = "Y") %>%
  select(substrate_ac = UniprotAccession,
         site_residue,
         site_position = Tyr
         ) %>% filter(site_position != "") -> iptm_format_ptk

utils::write.table(iptm_format_ptk %>% filter(site_position != ""),
                   "data-raw/iptm_net_input_ptk_chip.csv", col.names = F, sep = "\t", quote = F, row.names = F)

#iptm_format_list <- split(iptm_format, seq(nrow(iptm_format))) %>% map(as.list)
#iptmnetr::get_ptm_enzymes_from_list(sites) -> iptm_res

iptmnetr::get_ptm_enzymes_from_file("data-raw/iptm_net_input_ptk_chip.csv") -> iptm_res_ptk

iptm_res_ptk %>%
  filter(!enz_name %in% c("conventional protein kinase C", "IkappaB kinase complex (human)",
                          "AKT kinase", "aurora kinase", "IKBKG", "Abl fusion"
  )) %>%
  mutate(enz_name = case_when(
    enz_name == "FAK" ~ "PTK2",
    T ~ enz_name
  )) -> iptm_res_ptk

setdiff(iptm_res_ptk$enz_name %>% unique(), kinome_mp_file$hgnc_symbol) -> iptmnet_hgnc_missing
iptm_res_ptk %>% filter(!enz_name %in% iptmnet_hgnc_missing) -> iptm_res_ptk


iptm_res_ptk %>%
  select(head = enz_name, sub_id, site) %>%
  mutate(entry = paste0(sub_id, ";", site,"-p;u")) %>%
  select(head, entry) %>%
  group_by(head) %>%
  mutate(entry = list(unique(entry)), len = n()) %>% ungroup() %>%
  distinct() -> iptm_res_ptk

left_join(iptm_res_ptk, kinome_mp_file %>% select(hgnc_symbol, ptmsea_id), by = c("head" = "hgnc_symbol")) %>% View()

iptm_res_ptk <- setNames(as.list(as.data.frame(t(iptm_res_ptk))), iptm_res_ptk$head)

cmapR::write_gmt(iptm_res_ptk, "data-raw/PTM_SEA_DB/ptm_sea_iptmnet_mapping_ptk.gmt")

ptm_sea_iptmnet_mapping_ptk <- readLines("data-raw/PTM_SEA_DB/ptm_sea_iptmnet_mapping_ptk.gmt")


# iptm to gmt
"O43561-2;Y226-p;u"

iptm_res_ptk %>% select(enz_name, sub_id, site) %>%
  mutate(entry = paste0(sub_id, ";", site, "-p;u")) %>%
  select(head = enz_name, entry) %>% View()

setdiff(iptm_res_ptk$enz_name %>% unique(), kinome_mp_file$hgnc_symbol)


## STK iptment
stk_pamchip_87102_array_layout_ptmsea <- readxl::read_xlsx("data-raw/Pamchips_Layout/87102_ArrayLayout.xlsx") %>%
  filter(ID != "#REF") %>%
  select(Ser,Thr, substrate_ac = UniprotAccession) %>%
  mutate(Ser = ifelse(Ser == "[]", NA, Ser),
         Thr = ifelse(Thr == "[]", NA, Thr),
         Ser = str_replace(Ser, "\\[", ""),
         Ser = str_replace(Ser, "\\]", ""),
         Thr = str_replace(Thr, "\\[", ""),
         Thr = str_replace(Thr, "\\]", ""),
         ) %>%
    separate_rows(Ser, sep = ", ") %>%
    separate_rows(Thr, sep = ", ") %>%
  rename(S = Ser, `T` = Thr) %>%
  pivot_longer(cols = 1:2, names_to = "site_residue", values_to = "site_position") %>%
  filter(!is.na(site_position), substrate_ac != "NA") -> iptm_format_stk


utils::write.table(iptm_format_stk,
                   "data-raw/iptm_net_input_stk_chip.csv", col.names = F, sep = "\t", quote = F, row.names = F)


iptmnetr::get_ptm_enzymes_from_file("data-raw/iptm_net_input_stk_chip.csv") -> iptm_res_stk

iptm_res_stk %>%
  filter(!enz_name %in% c("conventional protein kinase C", "IkappaB kinase complex (human)",
                          "AKT kinase", "aurora kinase", "IKBKG"
                          )) %>%
  mutate(enz_name = case_when(
    enz_name == "hMAP3K7/Phos:1" ~ "MAP3K7",
    enz_name == "hAURKA" ~ "AURKA",
    enz_name == "hNUAK1" ~ "NUAK1",
    enz_name == "hCDK9" ~ "CDK9",
    enz_name == "hATM" ~ "ATM",
    enz_name == "hIKBKE" ~ "IKBKE",
    enz_name == "hAURKB" ~ "AURKB",
    enz_name == "Pdk1" ~ "PDK1",
    enz_name == "CHK2" ~ "CHEK2",
    T ~ enz_name
  )) -> iptm_res_stk

setdiff(iptm_res_stk$enz_name %>% unique(), kinome_mp_file$hgnc_symbol) -> iptmnet_hgnc_missing
iptm_res_stk %>% filter(!enz_name %in% iptmnet_hgnc_missing) -> iptm_res_stk


iptm_res_stk %>%
  select(head = enz_name, sub_id, site) %>%
  mutate(entry = paste0(sub_id, ";", site,"-p;u")) %>%
  select(head, entry) %>%
  group_by(head) %>%
  mutate(entry = list(unique(entry)), len = n()) %>% ungroup() %>%
  distinct() -> iptm_res_stk

iptm_res_stk <- setNames(as.list(as.data.frame(t(iptm_res_stk))), iptm_res_stk$head)

cmapR::write_gmt(iptm_res_stk, "data-raw/PTM_SEA_DB/ptm_sea_iptmnet_mapping_stk.gmt")

ptm_sea_iptmnet_mapping_stk <- readLines("data-raw/PTM_SEA_DB/ptm_sea_iptmnet_mapping_stk.gmt")



# format O60934;S343-p
ptk_pamchip_86402_array_layout %>%
  mutate(
         PTM_SEA_ID = paste0(substrate_ac, ";", site_residue, "-p")) %>%
  dplyr::select(substrate_ac, PTM_SEA_ID) -> ptk_pamchip_86402_array_layout_ptmsea


stk_pamchip_87102_array_layout %>%
  separate_rows(Ser, sep = ",") %>%
  separate_rows(Thr, sep = ",") %>%
  mutate(Ser = paste0("S", Ser),
         Thr = paste0("T", Thr)) %>%
  pivot_longer(2:3,names_to = "AA", values_to = "Position") %>%
  filter(nchar(Position) != 1) %>%
  select(-AA) %>%
  mutate(PTM_SEA_ID = paste0(UniprotAccession, ";", Position, "-p")) %>%
  dplyr::select(ID, PTM_SEA_ID) -> stk_pamchip_87102_array_layout_ptmsea

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


# filterDBS_ptk %>% mutate(ids = entry, ids = str_remove(ids, ";u"), ids = str_remove(ids, "-2")) %>%
#   filter(ids %in% ptk_pamchip_86402_array_layout_ptmsea$PTM_SEA_ID) %>%
#   select(head, entry) %>%
#   group_by(head) %>% mutate(len = n()) %>%
#   ungroup() %>%
#   group_by(head) %>%
#   mutate(entry = list(unique(entry))) %>% ungroup() %>%
#   distinct() %>%
#   select(head, entry, len) -> onlyChipPeps_dbs

# onlyChipPeps_dbs <- setNames(as.list(as.data.frame(t(onlyChipPeps_dbs))), onlyChipPeps_dbs$head)

# cmapR::write_gmt(onlyChipPeps_dbs, "data-raw/PTM_SEA_DB/ptk_pamchip_86402_onlyChipPeps_dbs.gmt")

ptk_pamchip_86402_onlyChipPeps_dbs <- readLines("data-raw/PTM_SEA_DB/ptk_pamchip_86402_onlyChipPeps_dbs.gmt")

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

stk_pamchip_87102_onlyChipPeps_dbs <- readLines("data-raw/PTM_SEA_DB/stk_pamchip_87102_onlyChipPeps_dbs.gmt")


ptmsea_all_dbs <- readLines("data-raw/PTM_SEA_DB/ptm.sig.db.all.uniprot.human.v1.9.0.gmt")

KRSA::KRSA_coverage_STK_PamChip_87102_v2$Substrates %>% unique() -> stk_peps
KRSA::KRSA_coverage_PTK_PamChip_86402_v1$Substrates %>% unique() %>% as.character() -> ptk_peps

setdiff(stk_peps, stk_pamchip_87102_mapping$ID)



usethis::use_data(uka_db_full,
                  kinome_mp_file,
                  kinome_mp_file_v1,
                  kinome_mp_file_v2,
                  kinome_mp_file_v3,
                  kinome_mp_file_v4,
                  kinome_mp_file_v5,
                  kinome_mp_file_v6,
                  stk_pamchip_87102_mapping,
                  stk_pamchip_87102_array_layout_ptmsea,
                  ptk_pamchip_86402_mapping,
                  ptk_pamchip_86402_array_layout_ptmsea,
                  ptk_pamchip_86402_onlyChipPeps_dbs,
                  stk_pamchip_87102_onlyChipPeps_dbs,
                  ptm_sea_iptmnet_mapping_ptk,
                  ptm_sea_iptmnet_mapping_stk,
                  ptmsea_all_dbs,
                  overwrite = T)
