# Get EffectorP statistics

# Author: MVAP
# Version: 2023-10-16

# Directory tree

directorytree

# Packages ----
# Remember to insert %>% is ctrl + Shift +M

library('tidyverse')
library('here')
library('janitor')
library('skimr')
library('dplyr')
library('readxl')
lapply(c("Biostrings", "GenomicRanges", "GenomicFeatures", "Rsamtools", "rtracklayer","data.table"), require, character.only = TRUE)
# Read files ----

# . Read effectorP output -----
prot_dire <-list.files(path="effectorP_seq/", recursive=TRUE, full.names = T)
prot_dire <- grep(prot_dire, pattern='effectorP.out', invert=FALSE, value=TRUE)
prot_file <-  str_split(str_split(prot_dire, "/", simplify = T)[,2], "_signalp", simplify = T)[,1]
fileList <- tibble(original=prot_dire, species=prot_file)

# . Read protein files -----

prot_dire <- list.files(path = "newSpecies/", recursive = TRUE, full.names = T)
#prot_dire <- list.files(path = ".")
prot_dire <- grep(prot_dire, pattern = 'protein_longest.faa$', invert = FALSE, value = TRUE)
prot_file <- str_split(prot_dire, "/", simplify = T)[,3]
#prot_file <- paste("../amyannualReview/", prot_dire, sep = "")
protList <- tibble(original=prot_dire, species=prot_file)
#fileList <- tibble(original = prot_file, species = prot_dire)


prot_dire <- list.files(path="../folLPMOs/proteomes_longestprot/",pattern = "faa|fasta", recursive = F, full.names = T)
prot_file <- str_split(prot_dire, "/", simplify = T)[,4]
protList <- protList %>% add_row(original = prot_dire, species = prot_file)

# . Dictionary ----

speciesList <- read_delim("speciesList_2.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

# Load cazyme JGI ----

path <- "cazydb.xlsx"
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
  mutate(TotalCazyme = rowSums(.[2:8], na.rm = TRUE))

# Code ----
# . get unique protome

secretome = c()
for (line in 1:nrow(protList)){
  
 pep_file <- Biostrings::readAAStringSet(protList[line,]$original)
  
 proteins <- data.table(as.data.frame(pep_file), keep.rownames = TRUE) 
 
 secretome <- rbind(secretome , 
                   data.frame(Species = protList[line,]$species, 
                              UniqueProteome = nrow(proteins)))
}

rm(pep_file, proteins)

# . get Effector output
effector = c()

for (line in 1:nrow(fileList)){
  nlines <- countLines(fileList[line,]$original)
  
  lfile <- read_tsv(fileList[line,]$original, 
                    skip = nlines[1]-10, col_names = FALSE) 
                    
  #lfile$species <- fileList[line,]$species
  
  
  total <- lfile %>%
    filter(str_detect(X1, "proteins")) %>% 
    separate(X1, c("TotalProt"), " ", extra = "drop")
  
  nc_effector <- lfile %>% 
    filter(str_detect(X1, "apoplastic effectors:")) %>% 
    separate(X1, c("algo", "ApoplasticEffectors"), ": ", extra = "drop") %>% 
    dplyr::select(-algo)
  
  effector <- rbind(effector,
                    data.frame(Species = fileList[line,]$species, 
                                         TotalProt = as.numeric(as.character(total$TotalProt)),
                                         ApoplasticEffectors = as.numeric(
                                           as.character(nc_effector$ApoplasticEffectors)
                                           )
                               )
                    )
  
  #print(effector)
}

#effector <- read.csv("effector_output.csv")

effector_long <- effector %>% 
  pivot_longer(c(-Species)) %>% 
  #mutate(value = as.numeric(value)) %>% 
  mutate(Species = str_replace(Species, "_protein_longest.faa|_CSFG_models_2014-01-21.faa|_JGI_FilteredModelsv2.0.proteins.fasta|_refseq|_genbank|_GeneCatalog_proteins_20200124.aa.fasta|_JGI_GeneCatalog_proteins_20130412.aa.fasta|_GeneCatalog_proteins_20131310.aa.fasta", ""))%>%
  mutate(Species = str_replace(Species,"_protein_longest.faa", "" )) %>% 
  left_join(speciesList, by = c("Species" = "File")) %>% 
  drop_na(value.y) %>% 
  left_join(cazyme_df, by =c(`JGI Acronims` = "Name")) %>% 
  dplyr::select(Species, name, value.x, lifestyle, GH, GT, CE, PL, AA, CBM, Expansins, TotalCazyme) %>% 
  replace(is.na(.), 0) %>% 
  arrange(lifestyle) %>% 
  mutate(row = row_number())

# ggplot(effector_long, aes(fill=name, x=value.x, y=reorder(Species, row))) + #, label = lifestyle)) + 
#   geom_bar(position="dodge", stat="identity") +
#   scale_x_continuous(name = "Number of unique proteins")+
#   theme(legend.title=element_blank())+
#   theme(legend.position = "top", legend.direction = "horizontal") +
#   facet_wrap(~lifestyle, 
#              scales = "free_y", 
#              switch = "y")


ggplot(effector_long, aes(fill=name, x=value.x, y=reorder(Species, row), group = lifestyle)) + 
  geom_bar(position="dodge", stat="identity") +
  #geom_text(position = position_dodge(width = 1), aes(y=lifestyle, x=0)) +
  scale_x_continuous(name = "Number of unique proteins")+
  theme(legend.title=element_blank())+
  theme(legend.position = "top", legend.direction = "horizontal") 


# heatmap Cazyme

tile_cazyme <- effector_long[c(1,5:11)] %>% 
  pivot_longer(-Species)

tile_cazyme$name <- factor(tile_cazyme$name, levels= c("GH", "GT", "PL", "CE", "AA", "CBM", "Expansins"))
tile_cazyme$Species <- factor(tile_cazyme$Species, levels= unique(tile_cazyme$Species))

hm <- ggplot(data = tile_cazyme, aes(x =factor(name), y = Species, fill = value)) + geom_tile() + 
  scale_fill_distiller(name = "Num. of proteins", palette = "Blues", direction = 1, na.value = "transparent") +
  scale_x_discrete(breaks = unique(tile_cazyme$name),
                   labels = unique(tile_cazyme$name), 
                   name = "CAZy Families") + 
  #theme_gray() +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.title = element_text(size = 10, ), 
        legend.title.align = 0,
        legend.key.size = unit(1,"cm"),
        legend.text = element_text(size = 8),)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

#hm

## Get legend
tmp <- ggplot_gtable(ggplot_build(hm))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]


#After this we will produce new heatmap without legend and axes by applying appropriate theme:
  
# Remove legend from heatmap
hm.clean <- hm +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position="none")
hm.clean

# total heatmap

total <- effector_long[c(1,12)] %>% 
  pivot_longer(-Species)

total$Species <- factor(total$Species, levels= unique(total$Species))

tm <- ggplot(data = total, aes(x = name, y = Species, fill = value)) + geom_tile() + 
  scale_fill_distiller(name = "Total num. \nof CAZymes", palette = "Greens", direction = 1, na.value = "transparent") +
  scale_x_discrete(breaks = unique(tile_cazyme$name),
                   labels = unique(tile_cazyme$name), 
                  name = "Total\nCAZymes") + 
  #theme_gray() +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.title = element_text(size = 10, ), 
        legend.title.align = 0,
        legend.key.size = unit(1,"cm"),
        legend.text = element_text(size = 8),)+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

tm

## Get legend
tmp3 <- ggplot_gtable(ggplot_build(tm))
leg3 <- which(sapply(tmp3$grobs, function(x) x$name) == "guide-box")
legend3 <- tmp3$grobs[[leg3]]


#After this we will produce new heatmap without legend and axes by applying appropriate theme:

# Remove legend from heatmap
tm.clean <- tm +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank(), 
        #axis.ticks.x = element_blank(),
        legend.position="none")
tm.clean

## Get bargraph

effect_prot <- effector_long[1:3] %>% 
  mutate (name = str_replace(name, "ApoplasticEffectors", "Apoplastic effectors")) %>% 
  mutate (name = str_replace(name,"TotalProt", "Secretome"))


effect_prot$Species <- factor(effect_prot$Species, levels= unique(effect_prot$Species))


bp.y <- ggplot(data = effect_prot, aes(fill= name , x = value.x, y = factor(Species))) + 
  geom_bar(stat = "identity", position = "dodge") + 
  #coord_flip() + theme_gray() +
  theme(#axis.title.x = element_blank(), axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(), axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 15, margin = margin(0,10,0,0)),
        legend.position="none") +
  scale_x_continuous(name = "Number of unique proteins")+
  theme(legend.title=element_blank())+
  theme(legend.position = "top", legend.direction = "vertical") 

#bp.y

## Get legend
tmp2 <- ggplot_gtable(ggplot_build(bp.y))
leg2 <- which(sapply(tmp2$grobs, function(x) x$name) == "guide-box")
legend.bp <- tmp2$grobs[[leg2]]


#After this we will produce new heatmap without legend and axes by applying appropriate theme:

# Remove legend from heatmap
bp.clean <- bp.y +
  theme(legend.position="none")
#bp.clean


legends <- plot_grid(legend.bp, legend, legend3, nrow = 3, ncol =1, 
                        align = "v")#, rel_heights = c(3, 1, 1))


grob.title <- textGrob("Main Title", hjust = 0.5, vjust = 0.5, gp = gpar(fontsize = 20))
grid.arrange(bp.clean, hm.clean, tm.clean, legends, nrow = 1, ncol = 4, 
             #widths = c(20, 20,2, 5),
             widths = c(22, 9,2, 5),
             #heights = c(40, 60), 
             top = grob.title)

g <- arrangeGrob(bp.clean, hm.clean, tm.clean, legends, nrow = 1, ncol = 4, 
                  widths = c(15, 9,1.5, 4),
                  #heights = c(40, 60), 
                  top = grob.title)

ggsave("proteome_cazyme.pdf", g, width = 12, height = 16, units = "in", dpi = 300)
