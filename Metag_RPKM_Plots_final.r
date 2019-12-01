#Elizabeth Wong
#Nov 26 2019
#Project 2 Contamination vs Completeness graph

#Set up
library(tidyverse)
library(RColorBrewer)
getwd()

#Load taxonomy datasets into gtdb_dat dataset
arc_class <- read.table("gtdbtk.ar122.classification_pplacer.tsv")
bac_class <- read.table("gtdbtk.bac120.classification_pplacer.tsv")
gtdb_dat <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  dplyr::select(mag, Kingdom, Phylum, Class, Order, Family)

addmargins(table(gtdb_dat$Kingdom))
ggplot(gtdb_dat, aes(x=Phylum)) +
  geom_bar(aes(fill=Kingdom)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle = 0),
        axis.text.x = element_text(angle=45, vjust=0.5, size=10)) +
  ggtitle("gtdb Output Taxonomy")

#Generating checkm data and plotting it
checkm_dat <- read.table("MetaBAT2_SaanichInlet_100m_min1500_checkM_stdout.tsv",
                         header=TRUE,
                         sep="\t",
                         comment.char = '') %>% 
  dplyr::rename(mag = Bin.Id) %>% 
  dplyr::select(mag, Completeness, Contamination)

ggplot(checkm_dat, aes(x=Completeness, y=Contamination)) +
  geom_point() +
  ggtitle("Checkm Output")

#Rename bins to mag
metag_rpkm <- read.table("SaanichInlet_100m_binned.rpkm.csv", header=T, sep=',') %>% 
  mutate(Sequence_name = gsub('m_', 'm.', Sequence_name)) %>% 
  mutate(Sequence_name = gsub('Inlet_', 'Inlet.', Sequence_name)) %>% 
  separate(col=Sequence_name, into=c("mag", "contig"), sep='_', extra="merge") %>% 
  group_by(Sample_name, mag) %>% 
  summarise(g_rpkm = sum(RPKM)) %>% 
  mutate(mag = gsub('Inlet.', 'Inlet_', mag))
metag_rpkm

#Join the metagenome RPKM, checkM and GTDB-tk data frames
rpkm_dat <- left_join(metag_rpkm, checkm_dat, by="mag") %>% 
  left_join(gtdb_dat, by="mag") %>% 
  group_by(mag, Kingdom, Phylum, Class, Order, Family, Completeness, Contamination) %>% 
  summarise(g_rpkm = mean(g_rpkm)) #take the mean of each MAG

#Plot Completeness vs Contaminatino for all MAGS
ggplot(rpkm_dat, aes(x=Completeness, y=Contamination, col=Phylum)) +
  geom_point(aes(size=g_rpkm)) +
  scale_size(range=c(1,10)) +
  xlim(c(50,100)) +
  ylim(c(0,100)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank()) +
  ggtitle("Completeness vs. Contamination Graph")

#Plot to visualize relative abundance
rpkm_dat <- left_join(metag_rpkm, checkm_dat, by="mag") %>% 
  left_join(gtdb_dat, by="mag") %>% 
  filter(Completeness > 50 & Contamination < 10) %>% # good quality MAGs
  mutate(mag = reorder(mag, Order, sort)) # sort by their taxonomic Order so everything shows up together

ggplot(rpkm_dat, aes(x=Sample_name, y=mag, col=Phylum)) +
  geom_point(aes(size=g_rpkm)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))

#Export high quality data as csv file
write.csv(as.data.frame(rpkm_dat), file = "rpkm_dat1.csv")

#Filter rpkm_dat to high quality data
high_quality <- select(rpkm_dat, mag, g_rpkm, Completeness, Contamination, Kingdom, Phylum, Class, Order, Family) %>% 
  filter(Completeness > 90) %>% 
  filter(Contamination < 5)
high_quality

#Exploring high quality data
addmargins(table(high_quality$Phylum))

ggplot(high_quality, aes(x=Phylum)) +
  geom_bar(aes(fill=mag)) +
  scale_colour_gradientn(colours=heat) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle = 0),
        axis.text.x = element_text(angle=45, vjust=0.5, size=10)) +
  ggtitle("High Quality MAGs Taxonomy")

#Plot high quality data Completeness vs Contamination graph
ggplot(high_quality, aes(x=Completeness, y=Contamination, col=Phylum)) +
  geom_point(aes(size=g_rpkm)) +
  scale_size(range=c(1,10)) +
  xlim(c(90,100)) +
  ylim(c(0,5)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank()) +
  ggtitle("High Quality MAGs Completeness vs. Contamination Graph")

#Plot high quality data relative abundance
ggplot(high_quality, aes(x=Sample_name, y=mag, col=Class)) +
  geom_point(aes(size=g_rpkm)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))

ggplot(high_quality, aes(x=mag, y=g_rpkm)) +
  geom_bar(stat = "identity", aes(fill=Phylum)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  ggtitle("High Quality MAGs Abundance")

#Export high quality data as csv file
write.csv(as.data.frame(high_quality), file = "high_quality1.csv")

#Mean RPKM
mean(high_quality$g_rpkm)

#Taxonomy data for high quality
Taxonomy_data <- select(high_quality, mag, g_rpkm, Kingdom, Phylum, Class, Order, Family)
write.csv(as.data.frame(Taxonomy_data), file = "Taxonomy_data1.csv")



