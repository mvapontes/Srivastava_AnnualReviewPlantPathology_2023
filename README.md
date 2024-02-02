Analysis and supplementary data for Srivastava 

### **1. Download data**

Download proteome and annotation (gtf/gff) for each species in the Supplementary file genomes, also found in the input folder. 

### **2. Select one protein per gene**

Run 01_longest_prot.R 

The output is a fasta file per species including the longest amino acid sequence per gene, rename by gene name. 1gene == 1protein

### **3. Ortholog analysis**

We used **[OrthoFinder](https://github.com/davidemms/OrthoFinder)** to identify ortholog groups. The concatenated multiple sequence aligment (MSA) of genes in orthogroups with only one copy and present in all species was used to infered the tree with **[MAFFT](https://mafft.cbrc.jp/alignment/software/)** v7.453 _[MBE 30:772-780 (2013)](https://doi.org/10.1093/molbev/mst010); [NAR 30:3059-3066 (2002)](https://doi.org/10.1093/nar/gkf436)_

### **4. Identify signal peptides with [SignalP 5.0](https://services.healthtech.dtu.dk/services/SignalP-5.0/)**
_[Nat Biotechnol 37:420-423 (2019)](https://doi.org/10.1038/s41587-019-0036-z)_


Use the following command (bash) to run **SignalP-5.0** locally per genome. 

```bash 
for i in `ls `; do signalp5 -fasta $i -format short > ../signalP_output/${i}_summary.signalp5 ; done
```

### **5. Select proteins with signal peptides identified by SignalP-5.0**

Run 03_filter_signalp.R

The output is a fasta file per species of proteins with signal peptide. 

### **6. Identify transmembrane topology and signal peptides with [Phobius ver 1.01](https://phobius.sbc.su.se/)** 
_[J Mol Bio 338(5):1027-1036 (2004)](https://doi.org/10.1016/j.jmb.2004.03.016)_

Use the following command (bash) to run **phobius-1.01** locally per file of signalp identified proteins. 

```bash 
for i in `ls ../signalP_output/*faa`; do phobius -short ${i} > ${i}_phobius.out ; done
```

### **7. Filter out proteins with signalp and transmembrane domains**

Run 04_filter_phobius.R

The output is a fasta file per species of proteins with signal peptide but without transmembrane domains. 

### **8. Identify effectors in extracellular fungal proteins with [EffectorP v3.0](https://effectorp.csiro.au/)** 
_[Mol Plant Microbe Interact 35:146-56 (2022)](https://doi.org/10.1094/MPMI-08-21-0201-R)_

Use the following command (bash) to run **EffectorP-3.0** locally per file

```python
python3 /mnt/e/software/EffectorP-3.0-main/EffectorP.py -f -i ../phobius/*_signalp.faa_phobius.faa > *.aa_signalp.faa_phobius.faa_fungal_effectorP.out
```

### **9a. Obtain CAZyme protein list from JGI**

Per JGI genomes, obtain the list of CAZymes and save it in an file, 1 file/sheet == 1 species. 

Run 02_jgigff_cazyme.R to make sure all CAZy protein ID are in gff file. 

### **9b. Predict CAZymes with [dbCAN](https://bcb.unl.edu/dbCAN2/)**
_[NAR 51:W115-W121 (2023)](https://doi.org/10.1093/nar/gkad328)_

Predict CAZyme proteins with **dbCAN3** for NCBI proteomes and those species from JGI without a list of CAZymes. 

### **10. Calculate number of proteins per species, 1gene == 1 protein**

```bash 
for i in `ls speciesSelected/*fa*`; do grep -iHEc '>' $i | sed 's/:/\t/g'; done > proteome_size.txt
```

### **11. Calculate number of proteins with signal protein and without transmembrane domain**, aka the secretome

```bash
for i in `ls tmhmm_phobius/*phobius.faa*`; do grep -iHEc '>' $i | sed 's/:/\t/g'; done > secretome_nosignal_nophobius.txt
```

### **12. Make the summary graphic**

Run 05_graph.R

![Figure 2](/plots/phylo_cazyme_effectors.tif)
