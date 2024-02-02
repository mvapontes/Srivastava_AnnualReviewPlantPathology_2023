#!/usr/bin/env Rscript

#  Script to obtain the longest protein per gene

# Author: MVAP

## To keep only the longest aminoacid sequence per gene, not the longest transcript
## if two proteins have the same length I keep one interchangeably

# Load all packages ----
lapply(c("Biostrings", 
         "GenomicRanges", 
         "GenomicFeatures", 
         "Rsamtools",
         "rtracklayer",
         "tidyr",
         "dplyr", 
         "data.table"), 
       require, 
       character.only = TRUE)

# Get all files ----
prot_dire <- list.files(path = "newSpecies/", 
                        recursive = TRUE)

# . get only gtf and faa ----
prot_dire <- grep(prot_dire, pattern = 'genomic.gtf|protein.faa',
                  invert=FALSE, 
                  value=TRUE)

fileList <- tibble(original = prot_dire, species = prot_dire)

# cool way to split anything after the last "." REMEMBER some of your species contain "." in their name
fileList <- fileList %>% 
  separate(species, into = c("species", "filename"), "/", remove = T) %>%
  mutate(filename = gsub(".gz", "", filename)) %>% 
  separate(filename, into = c(NA, "type"), "\\.(?=[a-z+])", remove = FALSE)

# get gtf and fasta in the same line with the species info
fileList <- fileList %>% 
  dplyr::select(-filename) %>%
  spread(type, original)

for (line in 1:nrow(fileList)){ #for each line 
	
  print("loading gtf file")
	
  # load gtf as a tible through rtracklayer that way is easier to wrangle your data
	gtf_file <- tibble::as_tibble(rtracklayer::import(paste("newSpecies/", fileList[line,]$gtf, sep = ""))) 
	
	print("loading fasta file")
	
	# load protein data as AAStringSet, easy loading easy writing
	pep_file <- Biostrings::readAAStringSet(paste("newSpecies/", fileList[line,]$faa, sep = ""))
	
	# get seq name and length provided from AAStringSet automaticaly. 
	# Split protein name, header contain protein_id species_id with a " " in between
	proteins <- tibble(longname = names(pep_file), shortname = names(pep_file) , width = width(pep_file)) %>% 
	  separate(shortname, into = c("shortname", NA), " ", extra = "drop")
	
	# select only lines with gene_id, transcript_id and protein_id AKA CDS
	# remove duplicated (in case any exon has the same length)
	# get gene_id frequency and add length from proteins file
	gtf_info <- gtf_file %>% 
	  filter(type == "CDS") %>% 
	  dplyr::select(gene_id, transcript_id, protein_id) %>% 
	  unique() %>% 
	  group_by(gene_id) %>% 
	  mutate(freq = n()) %>% 
	  left_join(proteins %>% 
	              select(shortname, width), 
	            by = c("protein_id" = "shortname"))

	# get the genes with the longest proteins
	gtf_max <- gtf_info %>% 
	  group_by(gene_id) %>% 
	  slice(which.max(width))
	
	# get only the proteins that have being descarted
	dup_proteins <- proteins %>% 
	  filter(!shortname %in% gtf_max$protein_id)
	
	print ("load table dataframe pepfile")
	
	# To avoid runing out of memory transform AAStringSet file to data.frame and to data.table
	set1.dt <- data.table(as.data.frame(pep_file), keep.rownames = TRUE)
	
	# get only the longest protein by discarting anything that is not in my duplitated list
	sel_fasta <- set1.dt[!set1.dt$rn %in% dup_proteins$longname,]
	
	# Now back to AAStringSet. 1 join in only line the aa seq
	set12 <- apply(sel_fasta[,!"rn"], 1, paste, collapse = "")
	
	# transformt data.table into AAStringSet
	set12 <- AAStringSet(set12)
	
	# add the names of the proteins
	names(set12) <- sel_fasta$rn
	
	# compose output file
	tmpfile <- paste("newSpecies/", fileList[line,]$species, "/", fileList[line,]$species, "_protein_longest.faa", sep = "")
	
	print (tmpfile)
	
	# write fasta file with the longest proteins per proteome
	writeXStringSet(set12, tmpfile)
	
}
