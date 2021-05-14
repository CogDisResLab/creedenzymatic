# Header Start ---------------------------
# User: Justin Creeden / justincreeden@gmail.com
# Date created: 2021_05_14
# Script: run creedenzymatic analysis
# Purpose: Testing
# Notes (if any): This script is part of test run
#
# Important: This header is automatically generated when user creates new document. Listed user (Justin Creeden) is not nessesarily author/creator.
# Header End ---


#Define varibles
JFC_directory_output <- "jfc/output/nick/AD_Females_DLPFC/"



# read and rank the KRSA table and use absolute values and descending sorting
read_krsa(krsa_ex, trns = "abs", sort = "desc") -> krsa_table_ranked

# read and rank the UKA table and use absolute values and descending sorting
read_uka(uka_ex, trns = "abs", sort = "desc") -> uka_table_ranked

# combine ranked tables
combine_tools(KRSA_df = krsa_table_ranked, UKA_df = uka_table_ranked) -> combined_df

# save file
write_delim(combined_df, paste0(JFC_directory_output,"ce_combined_ranked_file.txt"), delim = "\t")

# filter out kinases found in quartile 1 or 2 either in KRSA or UKA and use the quartile_figure() for visualization

combined_df %>% filter(Qrt <= 2) %>% pull(Uniprot_Gene) %>% unique() -> sig_kinases

combined_df %>% filter(Uniprot_Gene %in% sig_kinases) %>% quartile_figure()
