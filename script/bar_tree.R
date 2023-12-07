# Final graph

# Author: MVAP

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
library('R.utils')

# Read files ----

# . rooted tree from OrthoFinder Species tree ----
rtree <- read.tree("input/rooted_tree_newick.txt")

# . species info ----
spname <- read_tsv("input/speciesList.txt", col_names = TRUE) %>% 
  clean_names()

# . proteome data ----

proteome <- read_tsv("input/proteome_size.txt", col_names = c("filename", "proteome"))

# . secretome data ----

secretome <- read_tsv("input/secretome_nosignal_nophobius.txt", c("filename", "secretome"))

# . effectorP data ----

dire_effe <-list.files(path="effectorP_seq/", recursive=TRUE, full.names = T)
dire_effe <- grep(dire_effe, pattern='effectorP.out', invert=FALSE, value=TRUE)
file_effe <-  str_split(str_split(dire_effe, "/", simplify = T)[,2],
                        "_signalp", simplify = T)[,1]
fileList_effe<- tibble(original=dire_effe, species=file_effe)

#fileList_effe <- readRDS("fileList_effe.rds")

# . cazyme data ----

# . . JGI CAZyme ---- 
path <- "input/cazydbJGI.xlsx"
sheets <- excel_sheets(path = path)

df_list <- lapply(excel_sheets(path), function(x)
  read_excel(path, sheet = x)
)

names(df_list) <- sheets

cazyme <- df_list%>% 
  map(~ as_tibble(.x)) %>% 
  map2(names(.), ~ add_column(.x, Name = rep(.y, nrow(.x)))) %>% 
  bind_rows() %>% 
  clean_names %>% 
  filter(!str_detect(name, "Canalb1|Lopnu1|Maggr1|Fusoxcub1"))

# . . dbcan files ----

prot_dire <-list.files(path="dbcan_run/", recursive=TRUE, full.names = T)
prot_dire <- grep(prot_dire, pattern='overview.txt', invert=FALSE, value=TRUE)
prot_dire <- grep(prot_dire, pattern ="old_dbcan", invert = TRUE, value = TRUE)
prot_file <- str_replace(
  str_split(
    str_split(prot_dire, "/", simplify = T)[,2],
    ".faa|.fasta", simplify = T)[,1],
  "output_", "")

# FileList_dbcan
fileList <- tibble(original = prot_dire, species = prot_file)

#fileList <- readRDS("fileList.rds")
#dbcan <- readRDS("dbcan.rds")

# read dbcan output files 
dbcan <- c()

for (i in 1:nrow(fileList)){
  lfile <- read_tsv(fileList[i,]$original) %>% clean_names()

  lcazy <- lfile %>%
    mutate(filename = fileList[i,]$species) %>%
    filter(number_of_tools >=2) %>%
    separate_longer_delim(hmmer,delim = "+") %>%
    separate_longer_delim(diamond, delim = "+") %>%
    separate( hmmer, into = c("hmmer", NA), sep = "\\(", remove = TRUE, extra = "drop") %>%
    dplyr::select(gene_id, hmmer, diamond, filename) %>%
    gather(key = "db", value = "cazyme", hmmer, diamond) %>%
    mutate(across(cazyme, na_if, "-")) %>% drop_na() %>%
    dplyr::select(-db) %>% unique() %>%
    mutate(cazyme = str_extract(cazyme, "[aA-zZ]+")) %>% drop_na()

  dbcan <- rbind(dbcan, lcazy)
}

# Code ----

# . proteome numbers ----

prot_num <- proteome %>% 
  mutate(filename = str_replace(filename, ".fasta|.faa", "")) %>% 
  pivot_longer(proteome, names_to = "name", values_to = "value")

# . secretome numbers ----

secre_num <- secretome %>% 
  mutate(filename = str_replace(filename, "_signalp.faa_phobius.faa", "")) %>% 
  pivot_longer(secretome, names_to = "name", values_to = "value")

prot_scre <- prot_num %>%
  rbind(secre_num)

# . effector Total numbers ----

effectorTotal = c()

for (i in 1:nrow(fileList_effe)){
  nlines <- countLines(fileList_effe[i,]$original)

  # read tsv effector output
  lfile <- read_tsv(fileList_effe[i,]$original,
                    skip = nlines[1]-12, col_names = FALSE) # skip all files except last 12

  total <- lfile %>%
    filter(str_detect(X1, "proteins")) %>%
    separate(X1, c("TotalProt"), " ", extra = "drop")

  nc_effector <- lfile %>%
    mutate(filename = fileList_effe[i, ]$species) %>% # add species name
    filter(str_detect(X1, "Number of ")) %>%
    separate(X1, c("type", "n_effectors"), ": ") %>%
    mutate(type = case_when(type == "Number of predicted cytoplasmic effectors" ~ "Cytoplasmic", 
                            type == "Number of predicted apoplastic effectors" ~ "Apoplastic", 
                            .default = NA)) %>% 
    drop_na()


  effectorTotal <- rbind(effectorTotal, nc_effector)
                    
  #print(effector)
}

effectorTotal <- effectorTotal %>% 
  mutate_at('n_effectors', as.numeric)

# effectorTotal <- readRDS("effectorTotal.rds")

#. Effector Proteins ----

apoplasticProt = c()

for (i in 1:nrow(fileList_effe)){
  nlines <- countLines(fileList_effe[i,]$original)

  # read tsv effector output
  lfile <- read_tsv(fileList_effe[i,]$original, skip= 10, # skip first 10 lines 'blabla'
                    n_max = nlines[1]-12, col_names = TRUE) %>% # skip last 12 lines 
    clean_names() # clean names

  effprot <- lfile %>% drop_na() %>%
    filter(prediction != "Non-effector") %>% # filter out proteins classified as non-effector
    separate(number_identifier, into = c("gene_id",NA), sep = " ") %>%
    mutate(filename = fileList_effe[i, ]$species)

  apoplasticProt<- rbind(apoplasticProt, effprot)

}

#effectorProteins <- readRDS("effectorProteins.rds")

apoplasticProt_jgi <- apoplasticProt %>% 
  filter(str_detect(gene_id, "jgi")) %>% # filter out non jgi identifiers
  separate(gene_id, into = c(NA, "species", "gene", NA), sep = "\\|") %>% # split identifier for cazyme annotation 
  mutate(gene_id = paste0(species, gene)) %>% 
  dplyr::select(-species, -gene)

apoplasticProt_all <- apoplasticProt %>% filter(!str_detect(gene_id, "jgi")) %>% # filter out jgi proteins 
  rbind(apoplasticProt_jgi) # add jgi proteins



# . cazy numbers ----


cazyme_df <- cazyme %>% 
  mutate(protein_id = case_when(is.na(protein_id) ~ x, .default = protein_id)) %>% 
  dplyr::select(name, protein_id, ca_zy_annotations) %>% 
  drop_na(protein_id) %>% 
  # rename families
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
  mutate(CAZy = strsplit(ca_zy_annotations, "//")) %>% 
  unnest(CAZy) %>% 
  dplyr::select(name, protein_id, CAZy) %>% 
  # filter out empty lines and space by cazyme annotation
  filter(!CAZy == "") %>% 
  filter(!CAZy == " ") %>% 
  separate(CAZy, into = "CAZy", sep = " / ", extra = "drop") %>% 
  separate(CAZy, into = "Familia", sep = " ", extra = "drop", remove = FALSE) %>% # separate number
  group_by(name) %>% 
  count(Familia) %>% # count family (type)
  ungroup() %>% 
  pivot_wider( names_from = Familia, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  dplyr::select(-Expansins) %>% # remove expansins, dbcan doesnt identify expansins
  mutate(TotalCazyme = rowSums(.[2:7], na.rm = TRUE))


cazyme_enzy <- cazyme_df %>% 
  dplyr::select(-TotalCazyme) %>% 
  pivot_longer(-name, names_to = "CAZyme", values_to = "value")

# establish order by factor
cazyme_enzy$CAZyme <- factor(cazyme_enzy$CAZyme, levels = c("GH", "GT", "CE", "PL", "AA", "CBM"))

cazyme_total <- cazyme_df %>% 
  dplyr::select(name, TotalCazyme)

# . dbcan cazyme ----

# NCBI species
dbcan_df <- dbcan %>% 
  dplyr::select(-gene_id) %>% 
  group_by(filename) %>% 
  count(cazyme) %>% # count type
  ungroup() %>% 
  pivot_wider(names_from = cazyme, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(TotalCazyme = rowSums(.[2:7], na.rm = TRUE)) 

# JGI species without cazyme analysis
dbcan_jgi <- dbcan %>% 
  filter(str_detect(gene_id, "Canalb1|Fuscu1|Fusla1|Fusma1")) %>% 
  separate(gene_id, into =c(NA, "species", "prot", NA), extra  = "drop") %>% 
  mutate(gene_id = paste0(species, prot)) %>% 
  dplyr::select(gene_id, filename, cazyme)


dbcan_rename <- dbcan %>% 
  filter(!str_detect(gene_id, "Canalb1|Fuscu1|Fusla1|Fusma1")) %>% 
  rbind(dbcan_jgi)

# unify cazyme annotation and add to species info
sptreename_cazymes <- spname %>% 
  left_join(cazyme_df, by = c("jgi_acronims" = "name")) %>% 
  left_join(dbcan_df, by = c("file_name" = "filename")) %>% 
  mutate(AA = case_when(is.na(AA.y)~AA.x, .default = AA.y)) %>% 
  mutate(CE = case_when(is.na(CE.y)~AA.x, .default = CE.y)) %>% 
  mutate(CBM = case_when(is.na(CBM.y)~CBM.x, .default = CBM.y)) %>% 
  mutate(GH = case_when(is.na(GH.y)~GH.x, .default = GH.y)) %>% 
  mutate(GT = case_when(is.na(GT.y)~GT.x , .default = GT.y)) %>% 
  mutate(PL = case_when(is.na(PL.y)~PL.x, .default = PL.y)) %>% 
  mutate(TotalCazyme = case_when(is.na(TotalCazyme.y)~TotalCazyme.x, .default = TotalCazyme.y )) %>% 
  dplyr::select(file_name, GH, CE, GT, CBM, PL, AA, TotalCazyme)
  
prot_scre_cazyme <- sptreename_cazymes %>% 
  dplyr::select(file_name, TotalCazyme) %>% 
  rename(filename = file_name) %>% 
  pivot_longer(-filename, names_to = "name", values_to = "value") %>% 
  rbind(prot_scre)


all_families <- sptreename_cazymes %>% 
  dplyr::select(-TotalCazyme) %>%
  rename(filename = file_name) %>% 
  pivot_longer(-filename, names_to = "name", values_to = "value")

# . cazymeEffectors ----

# fix protein id for A niger
cazyme_niger <- cazyme %>% filter(str_detect(protein_id, "Aspni_NRRL3_" )) %>% 
  separate(protein_id, into = c("species", "gene") , sep = "_1") %>% 
  mutate(gene = str_pad(gene, 5, "0",side = "left" )) %>% 
  mutate(gene_id = paste0("NRRL3_", gene)) %>% 
  dplyr::select(-species, -gene)

# get number cazyme and effectors 
cazymeEffectors <- cazyme %>% 
  mutate(gene_id = case_when(is.na(protein_id) ~ x, .default = protein_id)) %>%
  filter(!str_detect(gene_id, "Aspni_NRRL3_")) %>% dplyr::select(-protein_id) %>% 
  rbind(cazyme_niger) %>% # add to cazyme file
  drop_na(gene_id) %>% 
  rename(filename = name) %>% 
  dplyr::select(gene_id, filename) %>%
  rbind(dbcan_rename %>% dplyr::select(gene_id, filename)) %>% 
  left_join(apoplasticProt_all, by = c("gene_id" = "gene_id")) %>% 
  drop_na(cytoplasmic_effector) %>% 
  separate(apoplastic_effector, into = c("aplo_class", "apo_value"), sep = " ", remove = FALSE) %>% 
  mutate(apo_value= as.numeric(
    str_replace(
      str_replace(apo_value, "\\(", ""),
      "\\)", ""))) %>% 
  separate(cytoplasmic_effector, into = c("cyto_class", "cyto_value"), sep = " ", remove = FALSE) %>% 
  mutate(cyto_value= as.numeric(
    str_replace(
      str_replace(cyto_value, "\\(", ""), 
      "\\)", ""))) %>% 
  mutate(effector = case_when(apo_value > cyto_value ~  "Apoplastic", 
                              cyto_value > apo_value ~ "Cytoplasmic", 
                              .default = "Apoplastic")) %>% 
  dplyr::select(filename.y, effector) %>% 
  rename(filename = filename.y) %>%
  filter(effector == "Apoplastic") %>% 
  group_by(filename) %>% count() %>% 
  rename(n_effectors = n) %>% 
  mutate(type = "CAZy - Apoplastic") %>% 
  rbind(effectorTotal)

# establish order by factor
cazymeEffectors$type<- factor(cazymeEffectors$type, levels = c("Cytoplasmic", "Apoplastic", "CAZy - Apoplastic"))


# . Tree data ----
sptreename <- spname %>%  
  mutate(strain = case_when(is.na(jgi_acronims) ~ paste0(strain, "*"),
                            str_detect(jgi_acronims, "Canalb1") ~ paste0(strain, "*"),
                            .default = strain)) %>% 
  dplyr::select(file_name, species, f_sp, race, var, strain, host2) %>% 
  mutate(lab = case_when(is.na(f_sp) & is.na(race) & is.na(var) &is.na(strain) ~ glue("italic({backtick(species)})"),
                       
                           is.na(f_sp) & is.na(race) & is.na(var) &!is.na(strain) ~ glue("italic({backtick(species)})~{backtick(strain)}"), 
                       !is.na(f_sp) & is.na(race) & is.na(var) & !is.na(strain) ~ 
                         glue("italic({backtick(species)})~'f. sp.'~italic({backtick(f_sp)})~{backtick(strain)}"),
                       
                       !is.na(f_sp) & !is.na(race) & is.na(var) & !is.na(strain) ~ 
                         glue("italic({backtick(species)})~'f. sp.'~italic({backtick(f_sp)})~race~{backtick(race)}~{backtick(strain)}"), 
                       
                       is.na(f_sp) & is.na(race) & !is.na(var) & !is.na(strain) ~ 
                         glue("italic({backtick(species)})~var.~{backtick(var)}~{backtick(strain)}"), 
                       .default  = file_name)
         )

# establish order
sptreename$host2<- factor(sptreename$host2 , level = c("Plant", "Mammal", "Insect", "Saprotroph"))

# Plots ----

# . get tree plot ----
p3 <- ggtree(rtree, 
             aes(color = host2), 
             size = 1, 
             show.legend = FALSE) %<+% sptreename + 
  
  scale_color_manual(values = c("Plant" = "forestgreen", 
                                           "Mammal" = "#FF0000",
                                           "Insect"=  "darkorchid2",
                                           "Saprotroph" = "darkred"), 
                                   na.value = "black", 
                     name = "Host") +
  
  geom_tiplab(aes(label = lab), 
              color = "black", 
              align = TRUE, 
              hjust = 1, 
              offset = 1, #aligned to right
              parse = T, # parse italics
              size = 3, 
              show.legend = FALSE) + 

  geom_fruit(data = sptreename %>% 
               dplyr::select(file_name, host2) %>%
               pivot_longer(-file_name),
             geom = geom_point,
             mapping = aes(y = file_name, 
                           x = name, 
                           color = value),
             size = 3,  
             offset = 0.225) +
  
  scale_fill_manual(values = c("Plant" = "forestgreen", 
                             "Mammal" = "#FF0000",
                             "Insect"=  "darkorchid2",
                             "Saprotroph"="darkred"), 
                    na.translate = F, 
                    name = "Host") + 
  
  theme_tree(legend.position = c(.3,.88), 
             legend.title = element_text(face = "bold"),
             legend.text = element_text(size = 10, face = "bold"))

# . get tree order for geom_col -----

tree_order <- p3b$data[[3]] %>% 
  dplyr::select(y,label) %>% # get tree order
  left_join(sptreename %>% 
              dplyr::select(file_name, lab), 
            by =c("label" = "lab"))

# establish order by tree label order
prot_scre$filename <- factor(prot_scre$filename, level = tree_order$file_name)

# . get proteome/secretome plot ----

ps_ggplot <- ggplot(prot_scre, 
                    aes(x = value, 
                        y = filename, 
                        group = name))+
  
  geom_col(aes(fill = name),
           position = position_stack(reverse = TRUE)) + # reverse order
  
  scale_fill_manual(values = c("proteome" = "darkblue", 
                               "secretome" ="darkorange2"), 
                    name = NULL, 
                    labels = c("Proteome", "Secretome")) +

  scale_x_continuous(position = "top", 
                     name = "Unique proteins", 
                     expand = c(0, 0), # avoid blank space 
                     limits = c(0, NA)) +
  
  theme_classic()+
  
  theme(legend.position = c(0.85, 0.76),
        legend.text = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(face = "bold"), 
        panel.grid.minor.x = element_line())

# establish order by tree label order
all_families$filename <- factor(all_families$filename, level= tree_order$file_name)

# . get cazyme plot ----
cazyme_ggplot <- ggplot(all_families, 
                        aes(x = value, 
                            y = filename, 
                            group = name)) +

  geom_col(aes(fill = name), 
           position = position_stack(reverse = TRUE)) + # reverse order

  scale_fill_manual(values = c("GH" = "darkblue", 
                               "GT" = "darkorange2", 
                               "CE" = "grey50",
                               "PL" = "darkgreen", 
                               "AA" = "darkred", 
                               "CBM" = "yellow3"), 
                    name = NULL) +

  scale_x_continuous(position = "top", 
                     name = "CAZymes",
                     expand = c(0, 0),
                     limits = c(0, NA)) +
  
  theme_classic()+
  
  theme(legend.position = c(0.85, 0.70), 
        axis.title.x = element_text(face = "bold"), 
        panel.grid.major.x = element_line(),
        legend.text = element_text(size = 10, 
                                   face = "bold"))


# establish order by tree label order
cazymeEffectors$filename <- factor(cazymeEffectors$filename, level= tree_order$file_name)

# . get effector plot ----

effector_ggplot <- ggplot(cazymeEffectors, 
                          aes(x= n_effectors, 
                              y = filename, 
                              group = type))+
  geom_col(aes(fill = type), 
           position = position_stack(reverse = TRUE),
           width = 0.9) +
  
  scale_fill_manual(values = c("Cytoplasmic" = "darkblue",
                               "Apoplastic" = "grey50", 
                               "CAZy - Apoplastic" = "forestgreen"), 
                    name = NULL) +

  scale_x_continuous(position = "top", 
                     name = "Effectors", 
                     expand = c(0, 0), 
                     limits = c(0, NA)) +
  
  theme_classic()+
  
  theme(legend.position = c(0.8, 0.74), 
        axis.title.x = element_text(face = "bold"), 
        panel.grid.major.x = element_line(),
        legend.text = element_text(size = 10, face = "bold"))

# . Combined plots ----
g <- grid.arrange(arrangeGrob(NULL, # To combine geom_col with the tree, the tree needs blank space top and bottom to fit plot title and axis
                              p3 + 
                                theme(plot.margin = margin(0,-20,0,0)), 
                              NULL, 
                              heights = c(0.028, 0.9685, 0.0035),
                              nrow = 3), 
             ps_ggplot + theme(axis.text.y = element_blank(), # remove y axis 
                               axis.title.y = element_blank(), # remove axis 
                               axis.ticks.y = element_blank(), axis.line.y = element_line()), 
             cazyme_ggplot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), # remove axis
                                   axis.ticks.y = element_blank(), axis.line.y = element_line()), 
             effector_ggplot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), # remove axis
                                     axis.ticks.y = element_blank(), axis.line.y = element_line()), 
             nrow = 1, ncol = 4, widths = c(40, 12,12,12))



ggsave("phylo_cazyme_effectors.pdf", g, width = 25, height = 15, units = "in", dpi = 300)


# orthogroups -----

ortho <- read_tsv("speciesSelected/OrthoFinder/Results_Oct20/Orthogroups/Orthogroups.tsv")

ortho_prot <- ortho %>% pivot_longer(-Orthogroup, names_to = "species", values_to = "proteinId", values_drop_na = TRUE) %>%
  separate_longer_delim(proteinId, delim = ", ")

cazyme_ortho <- apoplasticProt %>% 
  left_join(ortho_prot, by =c( "gene_id" = "proteinId")) %>%
  drop_na(Orthogroup) %>% 
  group_by(Orthogroup) %>% count()



  # geom_fruit(data = prot_num, geom = geom_text, 
  #            mapping= (aes(y= filename, x= name, label = value)), 
  #            size = 3)+
  # 
  # geom_fruit(data=prot_num, geom=geom_point,
  #            mapping=aes(y=filename, x=name, size = value), fill = "dodgerblue",
  #            color = "black", offset = 0.02, shape =21) +
  # scale_size(name = "Proteome size")+
  # 
  # # geom_fruit(data = secre_num, geom = geom_bar,
  # #            mapping = aes(y = filename, x= secretome),
  # #            pwidth=0.1,
  # #            stat="identity",
  # #            orientation="y", offset = 0.05,
  # #            axis.param=list(axis="x", nbreak=4, vjust=1, text.size =3, text.angle = 45))+
  # 
  # theme_tree(legend.position = c(.1,.88))


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
