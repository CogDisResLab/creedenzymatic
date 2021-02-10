library(tidyverse)

krsa_table <- read_delim("inst/extdata/KRSA_example.txt", delim = "\t")

uka_table <- read_delim("inst/extdata/UKA_example.txt", delim = "\t") %>% select(Kinase, Score)


perc <- ecdf(1:length(krsa_table$Score))

krsa_table %>%
  mutate(Rank = dense_rank( abs(Score)) ) %>% View()
  mutate(Method = "KRSA") %>%
  mutate(Qnt = perc(Rank),
         Qrt = case_when(
           Qnt <= 0.25 ~ 1,
           Qnt <= 0.5 ~ 2,
           Qnt <= 0.75 ~ 3,
           Qnt <= 1 ~ 4,
         )

           ) %>% View()




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

dd1 %>%
  select(Uniprot_Gene, KinaseFamily, Qrt, Method) %>%
  pivot_wider(names_from = Method, values_from = Qrt) %>%
  pivot_longer(3:4, names_to = "Method", values_to = "Qrt") %>%
  mutate(
         present = ifelse(is.na(Qrt), "1", "2"), Qrt = ifelse(present == "1", 2, Qrt)
  ) %>%
  filter(KinaseFamily %in% sigkins) %>%
  ggplot(aes(reorder(Uniprot_Gene, KinaseFamily), Method)) + geom_point(aes(size = Qrt, shape = present)) +
  #coord_flip() +
  scale_radius(trans = "reverse") +
  #scale_size(trans = "reverse", range = c(1,5), breaks = c(1,2,3,4)) +
  #geom_text(position = position_dodge(width = 1), aes(x=`Kinase Family`, y=5, label = `Kinase Family`)) +
  theme_bw() + facet_grid(. ~ KinaseFamily, scales = "free", space = "free") +
  scale_shape_manual(values=c(1, 19)) +
  theme(axis.text.x = element_text(angle = 30, size = 7.5, vjust = 0.7)) +
  labs(x = "", y = "")


