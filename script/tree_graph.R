# Final graph ammey annual review

# Author: MVAP
# Version: 2023-10-22

# Directory tree

directorytree

# Packages ----
# Remember to insert %>% is ctrl + Shift +M
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

library('tidyverse')
library('here')
library('janitor')
library('skimr')
library('dplyr')
library('ggtree')
library('glue')
library('ggtext')
library('readxl')
library('ggtreeExtra')


# Read files ----

# . rooted tree ----
rtree <- read.tree("rotted_tree_newick.txt")

# . species info ----
spname <- read_tsv("speciesList_4.txt", col_names = TRUE) %>% 
  clean_names()

# . proteome data ----

proteome <- read_tsv("proteome_size.txt", col_names = c("filename", "proteome"))

# . secretome data ----

secretome <- read_tsv("secretome_nosignal_nophobius.txt", c("filename", "secretome"))

# . effectorP dsata ----

# dire_effe <-list.files(path="effectorP_seq/", recursive=TRUE, full.names = T)
# dire_effe <- grep(dire_effe, pattern='effectorP.out', invert=FALSE, value=TRUE)
# file_effe <-  str_split(str_split(dire_effe, "/", simplify = T)[,2], "_signalp", simplify = T)[,1]
# fileList_effe<- tibble(original=dire_effe, species=file_effe)

fileList_effe <- readRDS("fileList_effe.rds")

# . cazyme data ----

# . . JGI CAZyme ---- 
path <- "cazydb_3.xlsx"
sheets <- excel_sheets(path = path)
df_list <- lapply(excel_sheets(path), function(x)
  read_excel(path, sheet = x)
)

names(df_list) <- sheets

cazyme <- df_list%>% 
  map(~ as_tibble(.x)) %>% 
  map2(names(.), ~ add_column(.x, Name = rep(.y, nrow(.x)))) %>% 
  bind_rows() 

# . . dbcan CAZyme ----

# prot_dire <-list.files(path="dbcan_run/", recursive=TRUE, full.names = T)
# prot_dire <- grep(prot_dire, pattern='overview.txt', invert=FALSE, value=TRUE)
# prot_dire <- grep(prot_dire, pattern ="old_dbcan", invert = TRUE, value = TRUE)
# prot_file <- str_replace(str_split(str_split(prot_dire, "/", simplify = T)[,2], ".faa|.fasta", simplify = T)[,1], "output_", "")
# fileList <- tibble(original=prot_dire, species=prot_file)

fileList <- readRDS("fileList.rds")
dbcan <- readRDS("dbcan.rds")

# dbcan <- c()
# 
# for (i in 1:nrow(fileList)){
#   lfile <- read_tsv(fileList[i,]$original) %>% clean_names()
#   
#   lcazy <- lfile %>% 
#     mutate(filename = fileList[i,]$species) %>% 
#     filter(number_of_tools >=2) %>% 
#     separate_longer_delim(hmmer,delim = "+") %>% 
#     separate_longer_delim(diamond, delim = "+") %>% 
#    separate( hmmer, into = c("hmmer", NA), sep = "\\(", remove = TRUE, extra = "drop") %>% 
#     dplyr::select(gene_id, hmmer, diamond, filename) %>% 
#     gather(key = "db", value = "cazyme", hmmer, diamond) %>% 
#     dplyr::select(-db) %>% unique() %>% 
#     mutate(across(cazyme, na_if, "-")) %>% drop_na() %>% 
#     mutate(cazyme = str_extract(cazyme, "[aA-zZ]+")) %>% drop_na()
#   
#   dbcan <- rbind(dbcan, lcazy)
# }

# Code ----

# . proteome numbers ----

prot_num <- proteome %>% 
  mutate(filename = str_replace(filename, ".fasta|.faa", "")) %>% 
  pivot_longer(proteome, names_to = "name", values_to = "value")

# . secretome numbers ----

secre_num <- secretome %>% 
  mutate(filename = str_replace(filename, "_signalp.faa_phobius.faa", ""))

# . effector Total numbers ----

# effectorTotal = c()
# 
# for (i in 1:nrow(fileList_effe)){
#   nlines <- countLines(fileList_effe[i,]$original)
#   
#   lfile <- read_tsv(fileList_effe[i,]$original, 
#                     skip = nlines[1]-12, col_names = FALSE) 
#   
#   #lfile$species <- fileList_effe[i,]$species
#   
#   
#   total <- lfile %>%
#     filter(str_detect(X1, "proteins")) %>% 
#     separate(X1, c("TotalProt"), " ", extra = "drop")
#   
#   nc_effector <- lfile %>% 
#     filter(str_detect(X1, "apoplastic effectors:")) %>% 
#     separate(X1, c("algo", "ApoplasticEffectors"), ": ", extra = "drop") %>% 
#     dplyr::select(-algo)
#   
#   effectorTotal <- rbind(effectorTotal,
#                     data.frame(Species = fileList_effe[i,]$species, 
#                                TotalProt = as.numeric(as.character(total$TotalProt)),
#                                ApoplasticEffectors = as.numeric(
#                                  as.character(nc_effector$ApoplasticEffectors)
#                                )
#                     )
#   )
#   
#   #print(effector)
# }

effectorTotal <- readRDS("effectorTotal.rds")

#. Effector Proteins ----

# effectorProteins = c()
# 
# for (i in 1:nrow(fileList_effe)){
#   nlines <- countLines(fileList_effe[i,]$original)
#   
#   lfile <- read_tsv(fileList_effe[i,]$original, skip= 10,
#                     n_max = nlines[1]-12, col_names = TRUE) %>% 
#     clean_names
#   
#   effprot <- lfile %>% drop_na() %>% 
#     filter(apoplastic_effector != "-") %>% 
#     separate(number_identifier, into = c("gene_id",NA), sep = " ") %>% 
#     mutate(filename = fileList_effe[i, ]$species)
#   
#   effectorProteins<- rbind(effectorProteins, effprot)
# 
# }

effectorProteins <- readRDS("effectorProteins.rds")

effectorProteins_jgi <- effectorProteins %>% 
  filter(str_detect(gene_id, "jgi")) %>% 
  separate(gene_id, into = c(NA, "species", "gene", NA), sep = "\\|") %>% 
  mutate(gene_id = paste0(species, gene)) %>% 
  dplyr::select(-species, -gene)

effectorProteins_all <- effectorProteins %>% filter(!str_detect(gene_id, "jgi")) %>% 
  rbind(effectorProteins_jgi)
  


# . cazy numbers ----


# cazyme_df <- cazyme %>% 
#   mutate(ProteinId = case_when(is.na(`Protein Id`) ~ `≠≠`, .default = `Protein Id`)) %>% 
#   dplyr::select(Name, ProteinId, `CAZy Annotations`) %>% 
#   drop_na(ProteinId) %>% 
#   mutate(protidNum = str_replace(ProteinId, Name, ""))


cazyme_df <- cazyme %>% 
  mutate(ProteinId = case_when(is.na(`Protein Id`) ~ `≠≠`, .default = `Protein Id`)) %>% 
  dplyr::select(Name, ProteinId, `CAZy Annotations`) %>% 
  drop_na(ProteinId) %>% 
  mutate_all(funs(str_replace_all(., "GlycosylTransferase Family", "//GT"))) %>% 
  mutate_all(funs(str_replace_all(., "Glycosyltransferase Family", "//GT"))) %>% 
  mutate_all(funs(str_replace_all(., "Carbohydrate-Binding Module Family", "//CBM"))) %>% 
  mutate_all(funs(str_replace_all(., "Glycoside Hydrolase Family", "//GH"))) %>% 
  mutate_all(funs(str_replace_all(., "Class II peroxidase", "//Peroxidase"))) %>% 
  mutate_all(funs(str_replace_all(., "Carbohydrate Esterase Family", "//CE"))) %>% 
  mutate_all(funs(str_replace_all(., "Auxiliary Activity Family", "//AA"))) %>%
  mutate_all(funs(str_replace_all(., "Polysaccharide Lyase Family", "//PL"))) %>%
  mutate_all(funs(str_replace_all(., "Glucooligosaccharide oxidase", "//AA 7"))) %>%
  mutate_all(funs(str_replace_all(., "Multicopper oxidase", "//AA 1"))) %>%
  mutate_all(funs(str_replace_all(., "Peroxidase", "//AA 2"))) %>%
  mutate_all(funs(str_replace_all(., "GMC oxidoreductase", "//AA 3"))) %>%
  mutate_all(funs(str_replace_all(., "PL14_dist", "//PL 14_dist"))) %>%
  mutate_all(funs(str_replace_all(., "PL42_dist", "//PL 42_dist"))) %>%
  mutate_all(funs(str_replace_all(., "AA16_dist", "AA 16"))) %>%
  mutate_all(funs(str_replace_all(., "Distantly related to plant expansins", " //Expansins"))) %>%
  mutate_all(funs(str_replace_all(., "Multicopper oxidase|Copper radical oxidase", "AA 1"))) %>% 
  mutate(CAZy = strsplit(`CAZy Annotations`, "//")) %>% 
  unnest(CAZy) %>% 
  dplyr::select(Name, ProteinId, CAZy) %>% 
  filter(!CAZy == "") %>% 
  filter(!CAZy == " ") %>% 
  separate(CAZy, into = "CAZy", sep = " / ", extra = "drop") %>% 
  separate(CAZy, into = "Familia", sep = " ", extra = "drop", remove = FALSE) %>% 
  group_by(Name) %>% 
  count(Familia) %>% 
  ungroup() %>% 
  pivot_wider( names_from = Familia, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(TotalCazyme = rowSums(.[2:8], na.rm = TRUE)) %>% 
  dplyr::select(-Expansins)

cazyme_enzy <- cazyme_df %>% 
  dplyr::select(-TotalCazyme) %>% 
  pivot_longer(-Name, names_to = "CAZyme", values_to= "value")

cazyme_enzy$CAZyme <- factor(cazyme_enzy$CAZyme, levels = c("GH", "GT", "CE", "PL", "AA", "CBM"))

cazyme_total <- cazyme_df %>% 
  dplyr::select(Name, TotalCazyme)

# . dbcan cazyme ----

dbcan_df <- dbcan %>% 
  dplyr::select(-gene_id) %>% 
  group_by(filename, cazyme) %>% count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = cazyme, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(TotalCazyme = rowSums(.[2:7], na.rm = TRUE)) 


# . cazymeEffectors ----

cazyme_niger <- cazyme %>% filter(str_detect(`Protein Id`, "Aspni_NRRL3_" )) %>% 
  separate(`Protein Id`, into = c("species", "gene") , sep = "_1") %>% 
  mutate(gene = str_pad(gene, 5, "0",side = "left" )) %>% 
  mutate(gene_id = paste0("NRRL3_", gene)) %>% 
  dplyr::select(-species, -gene)


cazymeEffectors <- cazyme %>% 
  rename(gene_id = `Protein Id`) %>% 
  mutate(gene_id = case_when(is.na(gene_id) ~ `≠≠`, .default = gene_id)) %>%
  filter(!str_detect(gene_id, "Aspni_NRRL3_")) %>% 
  rbind(cazyme_niger) %>% 
  drop_na(gene_id) %>% 
  rename(filename = Name) %>% 
  dplyr::select(gene_id, filename) %>%
  rbind(dbcan %>% dplyr::select(gene_id, filename)) %>% 
  left_join(effectorProteins_all, by = c("gene_id" = "gene_id")) %>% 
  drop_na(cytoplasmic_effector)

# . Tree data ----
sptreename <- spname %>% 
  dplyr::select(file_name, species, f_sp, race, var, strain, lifestyle) %>% 

  # mutate(pad = case_when(is.na(f_sp) & is.na(race) & is.na(var) & is.na(strain) ~ strrep(".", 82 - nchar(species)),
  # 
  #                        is.na(f_sp) & is.na(race) & is.na(var) &!is.na(strain) ~ strrep(".", 82 - nchar(paste(species, strain, sep = " "))), 
  #                        .default = ".")) %>% 
# 
#                          !is.na(f_sp) & is.na(race) & is.na(var) & !is.na(strain) ~
#                            glue("italic({backtick(species)})~'f. sp.'~italic({backtick(f_sp)})~{backtick(strain)}"),
# 
#                          !is.na(f_sp) & !is.na(race) & is.na(var) & !is.na(strain) ~
#                            glue("italic({backtick(species)})~'f. sp.'~italic({backtick(f_sp)})~race~{backtick(race)}~{backtick(strain)}"),
# 
#                          is.na(f_sp) & is.na(race) & !is.na(var) & !is.na(strain) ~ glue("italic({backtick(species)})~var.~{backtick(var)}~{backtick(strain)}"),
#                          .default  = file_name))

  mutate(lab = case_when(is.na(f_sp) & is.na(race) & is.na(var) &is.na(strain) ~ glue("italic({backtick(species)})"),
                         
                         is.na(f_sp) & is.na(race) & is.na(var) &!is.na(strain) ~
                           glue("italic({backtick(species)})~{backtick(strain)}"), 
                         
                         !is.na(f_sp) & is.na(race) & is.na(var) & !is.na(strain) ~ 
                           glue("italic({backtick(species)})~'f. sp.'~italic({backtick(f_sp)})~{backtick(strain)}"),
                         
                         !is.na(f_sp) & !is.na(race) & is.na(var) & !is.na(strain) ~ 
                           glue("italic({backtick(species)})~'f. sp.'~italic({backtick(f_sp)})~race~{backtick(race)}~{backtick(strain)}"), 
                         
                          is.na(f_sp) & is.na(race) & !is.na(var) & !is.na(strain) ~ glue("italic({backtick(species)})~var.~{backtick(var)}~{backtick(strain)}"), 
                          .default  = file_name))


  p3 <- ggtree(rtree, aes(color = lifestyle), size =1, show.legend= FALSE) %<+% sptreename  +
  xlim(NA, 6) + 
  geom_tiplab(aes(label=lab), color = "black", 
              align = TRUE, hjust = 1, offset = 1, #aligned to right
              parse = T,  size = 3, show.legend = FALSE) + 
  
  # geom_polygon(aes(fill = lifestyle, x = 0, y = 0)) +
  # scale_fill_discrete(na.translate = F, name = "Lyfe style") +
  # theme_tree(legend.position = c(.1,.88)) +
  
  geom_fruit(data = sptreename %>% dplyr::select(file_name, lifestyle) %>% pivot_longer(-file_name),
             geom = geom_point,
             mapping = aes(y = file_name, x = name , color = value),
             size =3, show.legend = FALSE, offset = 0.225) +
    
  geom_fruit(data = prot_num, geom = geom_text, 
             mapping= (aes(y= filename, x= name, label = value)), 
             size = 3)+

  geom_fruit(data=prot_num, geom=geom_point,
             mapping=aes(y=filename, x=name, size = value), fill = "dodgerblue",
             color = "black", offset = 0.02, shape =21) +
  scale_size(name = "Proteome size")+
  
  geom_fruit(data = secre_num, geom = geom_bar,
             mapping = aes(y = filename, x= secretome),
             pwidth=0.1,
             stat="identity",
             orientation="y", offset = 0.05,
             axis.param=list(axis="x", nbreak=4, vjust=1, text.size =3, text.angle = 45))+

    theme_tree(legend.position = c(.1,.88))

    
# legend lifestyle
  ggplot(sptreename)+
  geom_polygon(aes(fill = lifestyle, x = 0, y = 0)) +
    scale_fill_discrete(na.translate = F, name = "Lyfe style")


p <- p + new_scale_fill() +
  geom_fruit(data=dat2, geom=geom_tile,
             mapping=aes(y=ID, x=Sites, alpha=Abundance, fill=Sites),
             color = "grey50", offset = 0.04,size = 0.02)+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, 
                                            keyheight = 0.3, order=5)) +
  geom_fruit(data=dat3, geom=geom_bar,
             mapping=aes(y=ID, x=HigherAbundance, fill=Sites),
             pwidth=0.38, 
             orientation="y", 
             stat="identity",
  ) +
  scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000",
                             "#800000", "#006400","#800080","#696969"),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4))+
  geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1) +
  theme(legend.position=c(0.93, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  )
p

p + layout_rectangular() + 
  theme(legend.position=c(.05, .7))
