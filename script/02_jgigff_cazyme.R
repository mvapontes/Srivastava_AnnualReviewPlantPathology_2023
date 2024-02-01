# Prepare JGI cazyme protein id

# Author: MVAP

# Load all packages ----
library('tidyverse')
library('here')
library('janitor')
library('skimr')
library('dplyr')
library('readxl')

prot_dire <-list.files(path = ".", pattern = "gff", recursive=TRUE, full.names = T)
prot_dire <- grep(prot_dire, pattern='gff', invert=FALSE, value=TRUE)
prot_file <-  str_split(str_split(prot_dire, "/", simplify = T)[,2], "_GeneCatalog", simplify = T)[,1]
fileList <- tibble(original=prot_dire, species=prot_file)

# Load cazyme JGI ----

path <- "input/cazydbJGI.xlsx"
sheets <- excel_sheets(path = path)
df_list <- lapply(excel_sheets(path), function(x)
  read_excel(path, sheet = x)
)

names(df_list) <- sheets

cazyme <- df_list%>% 
  map(~ as_tibble(.x)) %>% 
  map2(names(.), ~ add_column(.x, Name = rep(.y, nrow(.x)))) %>% 
  bind_rows() 

cazyme_df <- cazyme %>% 
  mutate(ProteinId = case_when(is.na(`Protein Id`) ~ `≠≠`, .default = `Protein Id`)) %>% 
  dplyr::select(Name, ProteinId, `CAZy Annotations`) %>% 
  drop_na(ProteinId) %>% 
  mutate(protidNum = str_replace(ProteinId, Name, ""))

# rename Aspfu ----

for (i in 1:nrow(fileList)){
  gtf_file <- tibble::as_tibble(rtracklayer::import("Aspfu1_GeneCatalog_genes_20120808.gff")) %>% 
    drop_na(proteinId) %>% 
    dplyr::select(name, proteinId) %>% 
    unique() %>% 
    mutate(Name = paste(fileList[i,]$species, proteinId, sep = "")) %>% 
    left_join(cazyme_df, by = c("Name" = "ProteinId")) %>% drop_na()

}
