library(tidyverse)

krsa_table <- read_delim("inst/extdata/KRSA_example.txt", delim = "\t")
uka_table <- read_delim("inst/extdata/UKA_example.txt", delim = "\t") %>% select(Kinase, Score)

mp_file <- read_delim("inst/extdata/kinase_mapping.txt", delim = "\t",
                      na = c("N/A", "Not found by Ethan", "Not found by Jake", "Not found",
                             "Not Found By Ethan", "Not Found By Jake"
                      )
)


read_krsa(krsa_table) -> krsa_table_ranked

left_join(krsa_table_ranked, select(mp_file, KRSA, Uniprot_Gene), by = c("Kinase" = "KRSA")) %>%
  mutate(KinaseFamily = Kinase) -> full_krsa


read_uka(uka_table) -> uka_table_ranked
left_join(uka_table_ranked, select(mp_file, UKA, Uniprot_Gene, KRSA), by = c("Kinase" = "UKA")) %>%
  rename(KinaseFamily = KRSA) -> full_UKA


rbind(full_krsa, full_UKA) -> dd1

dd1 %>% filter(Qrt <= 1) %>% pull(KinaseFamily) %>% unique() -> sigkins

quartile_figure(dd1, sigkins)


