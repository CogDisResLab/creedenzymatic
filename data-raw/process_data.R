library(readr)

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

usethis::use_data(uka_db_full,
                  kinome_mp_file_v1,
                  overwrite = T)

