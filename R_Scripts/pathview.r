# load necessary packages
library(tidyr)
library(dplyr)
library(pathview)
library(RColorBrewer)
library(knitr)
library(tidyverse)
library(pheatmap)

# read dataset and rename columns
ko <- read.table("data/all.ko.cleaned.txt") %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(ko = V2)

metat_rpkm <-read.csv("data/_MAG_ORFs_RPKM.csv",
                        header = TRUE) %>% 
  dplyr::rename(orf = v1) %>% 
  dplyr::rename(rpkm = v2)

prokka_mag_map <-read.csv("data/Prokka_MAG_map.csv",
                      header = TRUE) %>% 
  dplyr::rename(prokka_id = v1) %>% 
  dplyr::rename(mag = v2)

arc_class <- read.table("data/gtdbtk.ar122.classification_pplacer.tsv", sep="\t")
bac_class <- read.table("data/gtdbtk.bac120.classification_pplacer.tsv", sep="\t")
gtdb_dat <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", 
                               "Order", "Family", "Genus", "Species"))

checkm_dat <- read.table("data/MetaBAT2_SaanichInlet_100m_min1500_checkM_stdout.tsv",
                         header=TRUE,
                         sep="\t",
                         comment.char = '') %>% 
  dplyr::rename(mag = Bin.Id) %>% 
  dplyr::select(mag, Completeness, Contamination)

# rename the bins
metag_rpkm <- read.table("data/SaanichInlet_100m_binned.rpkm.csv",
                         header=T, sep=',')
metag_rpkm <- metag_rpkm %>%
  mutate(Sequence_name = gsub('m_', 'm.', Sequence_name)) %>% 
  mutate(Sequence_name = gsub('Inlet_', 'Inlet.', Sequence_name))%>% 
  separate(col=Sequence_name, into=c("mag", "contig"), 
           sep='_', extra="merge") %>% 
  group_by(Sample_name, mag) %>% 
  summarise(g_rpkm = sum(RPKM)) %>% 
  mutate(mag = gsub('Inlet.', 'Inlet_', mag))

# check the tables
head(ko) %>% kable()
head(metat_rpkm) %>% kable()

# phylum count
gtdb_dat %>% 
  group_by(Phylum) %>% 
  summarise(count = n_distinct(mag)) %>% 
  kable()
gtdb_dat <- dplyr::select(gtdb_dat, mag, Kingdom, Phylum, Class,
                          Order, Family)

# join tables
rpkm_dat <- left_join(ko, metat_rpkm, by = "orf") %>%
  separate(orf, into=c("prokka_id", "orf_id")) %>% # Split the Prokka ORF names into MAG identifier and ORF number for joining
  left_join(prokka_mag_map, by="prokka_id") %>% 
  left_join(gtdb_dat, by="mag") %>% 
  left_join(checkm_dat, by="mag")

# check the merged table
head(rpkm_dat) %>% kable()

# subset by proterobacteria
rpkm_dat_pro <- filter(rpkm_dat, 
                       Phylum == "p__Proteobacteria")
ko_rpkm_pro <- rpkm_dat_pro %>% 
  group_by(mag, ko) %>% 
  summarise(t_rpkm = sum(rpkm)) %>% 
  spread(key = mag, value = t_rpkm)
pv_mat_pro <- dplyr::select(ko_rpkm_pro, -ko)
rownames(pv_mat_pro) <- ko_rpkm_pro$ko

# subset alpha/gamma proteobacteria
rpkm_dat_alpha <- rpkm_dat %>% 
  filter(Class == "c__Alphaproteobacteria")
rpkm_dat_gamma <- rpkm_dat %>% 
  filter(Class == "c__Gammaproteobacteria")

# convert table for alpha/gamma proteobacteria
ko_rpkm_alpha <- rpkm_dat_alpha %>% 
  group_by(mag, ko) %>% 
  summarise(t_rpkm = sum(rpkm)) %>% 
  spread(key = mag, value = t_rpkm)
pv_mat_alpha <- dplyr::select(ko_rpkm_alpha, -ko)
rownames(pv_mat_alpha) <- ko_rpkm_alpha$ko

ko_rpkm_gamma <- rpkm_dat_gamma %>% 
  group_by(mag, ko) %>% 
  summarise(t_rpkm = sum(rpkm)) %>% 
  spread(key = mag, value = t_rpkm)
pv_mat_gamma <- dplyr::select(ko_rpkm_gamma, -ko)
rownames(pv_mat_gamma) <- ko_rpkm_gamma$ko

# subset individual mags, MAG#168 as example
rpkm_dat_168 <- filter(rpkm_dat, 
                        mag == "SaanichInlet_100m.168")
ko_rpkm_168 <- rpkm_dat_168 %>% 
  group_by(mag, ko) %>% 
  summarise(t_rpkm = sum(rpkm)) %>% 
  spread(key = mag, value = t_rpkm)
pv_mat_168 <- dplyr::select(ko_rpkm_168, -ko)
rownames(pv_mat_168) <- ko_rpkm_168$ko

# subset & spread by completeness and contamination
ko_rpkm_905 <- rpkm_dat %>% 
  filter(Completeness >= 90 & Contamination < 5) %>% 
  group_by(mag, ko) %>% 
  summarise(t_rpkm = sum(rpkm)) %>% 
  spread(key = mag, value = t_rpkm)
pv_mat_905 <- dplyr::select(ko_rpkm_905, -ko)
rownames(pv_mat_905) <- ko_rpkm_905$ko

# Nitrogen metabolism of >90% complete, 5% contaminated bins
pv.out_905 <- pathview(gene.data = pv_mat_905,
                       limit = list(gene = c(0,5)),
                   low = list(gene = "#91bfdb"),
                   mid = list(gene = "#ffffbf"),
                   high = list(gene = "#fc8d59"),
                   species = "ko",
                   pathway.id="00910",
                   kegg.dir = "Nmetabolism_955/")

# Nitrogen metabolism of Proteobacteria
pv.out_pro <- pathview(gene.data = pv_mat_pro,
                         limit = list(gene = c(0,5)),
                         low = list(gene = "#91bfdb"),
                         mid = list(gene = "#ffffbf"),
                         high = list(gene = "#fc8d59"),
                         species = "ko",
                         pathway.id="00910",
                         kegg.dir = "Nmetabolism_955/")

# Nitrogen metabolism of Alpha/gammaproteobacteria
pv.out_alpha <- pathview(gene.data = pv_mat_alpha,
                       limit = list(gene = c(0,5)),
                       low = list(gene = "#91bfdb"),
                       mid = list(gene = "#ffffbf"),
                       high = list(gene = "#fc8d59"),
                       species = "ko",
                       pathway.id="00910",
                       kegg.dir = "Nmetabolism_955/")

pv.out_gamma <- pathview(gene.data = pv_mat_gamma,
                            limit = list(gene = c(0,5)),
                            low = list(gene = "#91bfdb"),
                            mid = list(gene = "#ffffbf"),
                            high = list(gene = "#fc8d59"),
                            species = "ko",
                            pathway.id="00910",
                            kegg.dir = "Nmetabolism_955/")

# Nitrogen metabolism of MAG#168
pv.out_168 <- pathview(gene.data = pv_mat_168,
                      limit = list(gene = c(0,5)),
                      low = list(gene = "#91bfdb"),
                      mid = list(gene = "#ffffbf"),
                      high = list(gene = "#fc8d59"),
                      species = "ko",
                      pathway.id="00910",
                      kegg.dir = "Nmetabolism_955/")

# other pathway codes:
# oxidative phosphorylation 00190
# sulfur metabolism 00920
# methane metabolism 00680
# carbon fixation 00720
# photosynthesis 00195