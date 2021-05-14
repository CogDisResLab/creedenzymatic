# Header Start ---------------------------
# User: Justin Creeden / justincreeden@gmail.com
# Date created: 2021_05_13
# Script: import input data
# Purpose: Testing
# Notes (if any): This script is part of test run
#
# Important: This header is automatically generated when user creates new document. Listed user (Justin Creeden) is not nessesarily author/creator.
# Header End ---


#Define variables
JFC_directory_uka_input <- "jfc/input/hinds/2021-03-04-stk-h_f-v-h_c-UKA-Summaryresults 20210304-1330.txt"

JFC_directory_krsa_input <- "jfc/input/hinds/2021-05-14-stk-h_f-v-h_c-acrossChip_KRSA_FullTable_comp1.txt"

#Khaled, I changed your original code from 'median' to 'mean'
JFC_desired_UKA_metric <- "Mean Final Score"
#JFC_desired_UKA_metric <- "Mean Kinase Statistic"

#Import and process
krsa_ex <- read_delim(JFC_directory_krsa_input, delim = "\t")
krsa_ex %>% select(Kinase, Z) %>% rename(Score = Z) -> krsa_ex
uka_ex <- read_delim(JFC_directory_uka_input, delim = "\t")
uka_ex %>% select("Kinase Name", JFC_desired_UKA_metric) %>% rename(Kinase = "Kinase Name", Score = JFC_desired_UKA_metric) -> uka_ex


