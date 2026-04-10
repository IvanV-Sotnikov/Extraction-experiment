################################################################################
# R script for Sotnikov et al. (2026)
# Large soil samples for large soil organisms: modification of the eDNA
# extraction method for investigating soil invertebrate communities
################################################################################

# Original bioinformatic processing with DADA2 followed the tutorial:
# https://benjjneb.github.io/dada2/tutorial.html

# Clustering was performed with core vsearch functions:
# https://github.com/torognes/vsearch

library("ggplot2")
library("ggpubr")
library("ggrepel")
library("ggvegan")
library("ggVennDiagram")
library("ggh4x")
library("vegan")
library("dplyr") 
library("tidyr") 
library("tibble") 
library("reshape2")
library("scales")
library("car")
library("readxl")
library("readr")
library("stringr")
library("microDecon")
library("multcompView")
library("writexl")

rm(list = ls())
# Add paths to the working directories
path_to_source <- 'Path to source data'
path_to_results <- 'Path ro results'

SampleTable_18S <- read_excel(paste0(path_to_source, 'SampleTable_18S.xlsx'))
identification_18S <- read_excel(paste0(path_to_source, 'Identification_18S_PR2.xlsx'))
vsearch_18S <- read.table(paste0(path_to_source, 'asv_18S_99_to_otu.uc'), header = FALSE, comment.char = '#')
SampleTable_COI <- read_excel(paste0(path_to_source, 'SampleTable_COI.xlsx'))
identification_COI <- read_excel(paste0(path_to_source, 'Identification_COI.xlsx'))
vsearch_COI <- read.table(paste0(path_to_source, 'asv_COI_97_to_otu.uc'), header = FALSE, comment.char = '#')
Morphology <- read_excel(paste0(path_to_source, 'Morphology.xlsx'))
Concentration <- read_excel(path = paste0(path_to_source, 'Concentrations.xlsx'))

#
####                              Morphology data                            ####
Morphology_Genus <- Morphology %>%
  filter(!is.na(Genus))

Morphology_Family <- Morphology %>%
  group_by(Phylum, Class, Order, Family) %>%
  summarise(across(c(2:6), sum)) %>%
  filter(!is.na(Family))

Morphology_Order <- Morphology %>%
  group_by(Phylum, Class, Order) %>%
  summarise(across(c(3:7), sum)) %>%
  filter(!is.na(Order))

####                              Rarefaction of sequences                   ####
# Set seed for reproducibility
set.seed(123)

rarefied_coi_df <- SampleTable_COI %>%
  select(-c(length, NC_F_filt.fastq, PO_3_F_filt.fastq, K_filt_F_filt.fastq, K_0.5_F_filt.fastq,
            K_Max_F_filt.fastq)) %>%
  rename_with(~ gsub("_F_filt.fastq", "", .))

rarefied_18s_df <- SampleTable_18S %>%
  # read.csv2('Results/eDNA/Extract_Test/results from rawdata/18S/SampleTable_18S.csv') %>%
  select(-c(length, NC_F_filt.fastq, PO_3_F_filt.fastq, K_filt_F_filt.fastq, K_0.5_F_filt.fastq,
            K_Max_F_filt.fastq)) %>%
  rename_with(~ gsub("_F_filt.fastq", "", .))

# Perform rarefaction for COI
metadata_COI <- rarefied_coi_df[, c("ASVNumber", "sequence")]
counts_COI <- rarefied_coi_df[, !names(rarefied_coi_df) %in% c("ASVNumber", "sequence")]
counts_t_COI <- as.data.frame(t(counts_COI))
sample_depth_COI <- min(rowSums(counts_t_COI))
counts_rarified_COI <- rrarefy(counts_t_COI, sample = sample_depth_COI)
counts_rarified_t_COI <- as.data.frame(t(counts_rarified_COI))
keep_rows_COI <- rowSums(counts_rarified_t_COI) > 0

# Perform rarefaction for 18S
metadata_18S <- rarefied_18s_df[, c("ASVNumber", "sequence")]
counts_18S <- rarefied_18s_df[, !names(rarefied_18s_df) %in% c("ASVNumber", "sequence")]
counts_t_18S <- as.data.frame(t(counts_18S))
sample_depth_18S <- min(rowSums(counts_t_18S))
counts_rarified_18S <- rrarefy(counts_t_18S, sample = sample_depth_18S)
counts_rarified_t_18S <- as.data.frame(t(counts_rarified_18S))
keep_rows_18S <- rowSums(counts_rarified_t_18S) > 0

# Create a rarefied table for COI
rarefied_coi_ready <- cbind(
  metadata_COI[keep_rows_COI, ],
  counts_rarified_t_COI[keep_rows_COI, ])

# Create a rarefied table for 18S
rarefied_18s_ready <- cbind(
  metadata_18S[keep_rows_18S, ],
  counts_rarified_t_18S[keep_rows_18S, ])

####                              OTU preparation                            ####

# Create OTUs for COI
uc_file_COI <- vsearch_COI

centroid_mapping_COI <- uc_file_COI %>%
  filter(V1 == "S") %>%
  select(centroid_ASV = V9) %>%
  mutate(OTU_name = paste0("OTU_", row_number()))

asv_to_otu_final_COI <- uc_file_COI %>%
  filter(V1 %in% c("S", "H")) %>%
  select(type = V1, ASVNumber = V9, centroid_ASV = V10) %>%
  mutate(centroid_ASV = ifelse(type == "S", ASVNumber, centroid_ASV)) %>%
  left_join(centroid_mapping_COI, by = "centroid_ASV") %>%
  select(OTU_name, ASVNumber, centroid_ASV) %>%
  filter(ASVNumber %in% rarefied_coi_ready$ASVNumber)

# The final 97% OTU table will be prepared after ASV decontamination in the
# Taxonomic table section

# Clustering for 18S
uc_file_18S <- vsearch_18S

centroid_mapping_18S <- uc_file_18S %>%
  filter(V1 == "S") %>%
  select(centroid_ASV = V9) %>%
  mutate(OTU_name = paste0("OTU_", row_number()))

asv_to_otu_final_18S <- uc_file_18S %>%
  filter(V1 %in% c("S", "H")) %>%
  select(type = V1, ASVNumber = V9, centroid_ASV = V10) %>%
  mutate(centroid_ASV = ifelse(type == "S", ASVNumber, centroid_ASV)) %>%
  left_join(centroid_mapping_18S, by = "centroid_ASV") %>%
  select(OTU_name, ASVNumber, centroid_ASV) %>%
  filter(ASVNumber %in% c(rarefied_18s_ready$ASVNumber)) %>%
  rename(OTU_Number = OTU_name) %>%
  merge(identification_18S, by = 'ASVNumber') %>%
  merge(SampleTable_18S, by = 'ASVNumber') %>%
  group_by(OTU_Number, Phylum, Class, Order, Family, Genus, Species) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup() %>%
  distinct(OTU_Number, .keep_all = TRUE)

#
####                              MicroDecon                                 ####
for_Decon_18S <- asv_to_otu_final_18S %>%
  select(-c(length, PO_3_F_filt.fastq)) %>%
  rename_with(~ str_remove(., "_F_filt\\.fastq"), everything()) %>%
  tidyr::unite(lineage, Phylum, Class, Order, Family, Genus, Species, sep = ';') %>%
  relocate(OTU_Number, K_0.5, K_filt, K_Max, NC,
           A1_0.5, A2_0.5, A3_0.5, A4_0.5, B1_0.5, B2_0.5, B3_0.5, B4_0.5, C1_0.5, C2_0.5, C3_0.5, C4_0.5, D1_0.5, D2_0.5, D3_0.5, D4_0.5,
           A1_Max, A2_Max, A3_Max, A4_Max, B1_Max, B2_Max, B3_Max, B4_Max, C1_Max, C2_Max, C3_Max, C4_Max, D1_Max, D2_Max, D3_Max, D4_Max,
           E1_Filt, A2_Filt, A3_Filt, A5_Filt, B1_Filt, B2_Filt, B3_Filt, B4_Filt, C1_Filt, C2_Filt, C3_Filt, C4_Filt, D1_Filt, D2_Filt, D3_Filt, D4_Filt,
           lineage) %>%
  as.data.frame()

MicroDec_18S <- decon(data = for_Decon_18S, numb.blanks = 4, numb.ind = c(16, 16, 16), taxa = TRUE, runs = 2)

Decontaminant_18S <- MicroDec_18S$decon.table %>%
  select(-Mean.blank) %>%
  mutate(across(where(is.numeric), ~ ifelse(. <10, 0, .))) %>%
  filter(rowSums(across(where(is.numeric), ~ .x != 0)) > 0) %>%
  separate(lineage, into = c("Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  filter(Phylum != 'Metazoa')

for_Decon_COI <- merge(SampleTable_COI, identification_COI, by = 'ASVNumber') %>%
  select(-PO_3_F_filt.fastq) %>%
  filter(ASVNumber %in% c(rarefied_coi_ready$ASVNumber)) %>%
  pivot_longer(-c('ASVNumber', 'sequence', "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  separate(name, into = c('Sample', 'Method'), sep = '_') %>%
  tidyr::unite(SAMPLE, Sample, Method) %>%
  group_by(ASVNumber, SAMPLE, Species) %>%
  summarise(Reads = sum(value)) %>%
  select(ASVNumber, SAMPLE, Species, Reads) %>%
  dcast(ASVNumber+Species~SAMPLE, value.var = 'Reads') %>%
  select(-length_NA, -Sim_NA, -Cover_NA, - length_vsearch) %>%
  relocate(ASVNumber, K_0.5, K_filt, K_Max, NC_F,
           A1_0.5, A2_0.5, A3_0.5, A4_0.5, B1_0.5, B2_0.5, B3_0.5, B4_0.5, C1_0.5, C2_0.5, C3_0.5, C4_0.5, D1_0.5, D2_0.5, D3_0.5, D4_0.5,
           A1_Max, A2_Max, A3_Max, A4_Max, B1_Max, B2_Max, B3_Max, B4_Max, C1_Max, C2_Max, C3_Max, C4_Max, D1_Max, D2_Max, D3_Max, D4_Max,
           E1_Filt, A2_Filt, A3_Filt, A5_Filt, B1_Filt, B2_Filt, B3_Filt, B4_Filt, C1_Filt, C2_Filt, C3_Filt, C4_Filt, D1_Filt, D2_Filt, D3_Filt, D4_Filt,
           Species)

MicroDec_COI <- decon(data = for_Decon_COI, numb.blanks = 4, numb.ind = c(16, 16, 16), taxa = TRUE, runs = 2)

Decontaminant_COI <- MicroDec_COI$decon.table %>%
  select(-Mean.blank, -Species) %>%
  mutate(across(where(is.numeric), ~ ifelse(. == 1, 0, .))) %>%
  filter(rowSums(across(where(is.numeric), ~ .x != 0)) > 0)

####                              Taxonomic tables of COI                    ####

# Final COI OTU table
OTU97_COI_df <- merge(asv_to_otu_final_COI, Decontaminant_COI, by = 'ASVNumber') %>%
  group_by(OTU_name) %>%
  summarise(across(where(is.numeric), sum)) %>%
  rename_with(~ str_remove(., "_F_filt\\.fastq"), everything()) %>%
  # Remove rows with only zero values
  filter(rowSums(across(where(is.numeric), ~ .x != 0)) > 0)

identification_COI <- identification_COI %>%
  filter(!Class %in% c('Ostracoda', 'Polychaeta'))

COI_Species <- merge(identification_COI, Decontaminant_COI, by = 'ASVNumber') %>%
  filter(round(Sim) >= 97,
         !Species %in% c('NA NA', 'NA')) %>%
  separate(Species, into = c('gen', 'sp'), sep = ' ') %>%
  filter(sp != 'sp.') %>%
  tidyr::unite(Species, gen, sp, sep = ' ') %>%
  filter(!Order %in% c('Littorinimorpha')) %>%
  filter(!is.na(Species)) %>%
  filter(!str_detect(Species, "sp\\.$")) %>%
  group_by(Phylum, Class, Order, Family, Genus, Species) %>%
  summarise(across(c(5:52), sum), .groups = 'drop') %>%
  filter(rowSums(across(where(is.numeric), ~ .x != 0)) > 0)

COI_Genus <- merge(identification_COI, Decontaminant_COI, by = 'ASVNumber') %>%
  filter(round(Sim) >= 95,
         !Family %in% c('NA NA', 'NA', 'Nemesiidae', 'Eupterotidae', 'Syllidae'),
         !Order %in% c('Littorinimorpha', 'Cumacea', 'Cyclopoida', 'Decapoda', 'Eunicida',
                       'Haplotaxida', 'Octopoda', 'Cumacea', 'Phyllodocida'),
         !Genus %in% c('NA NA', 'NA'),
         !Order %in% c('Littorinimorpha')) %>%
  group_by(Phylum, Class, Order, Family, Genus) %>%
  summarise(across(c(6:53), sum), .groups = 'drop') %>%
  filter(rowSums(across(where(is.numeric), ~ .x != 0)) > 0) %>%
  merge(Morphology_Genus, by = c('Phylum', 'Class', 'Order', 'Family', 'Genus'), all = TRUE) %>%
  mutate(across(everything(), ~ replace_na(., 0)))

COI_Family <- merge(identification_COI, Decontaminant_COI, by = 'ASVNumber') %>%
  filter(Sim >= 90,
         !Family %in% c('NA NA', 'NA', 'Nemesiidae', 'Eupterotidae', 'Syllidae'),
         !Order %in% c('Littorinimorpha', 'Cumacea', 'Cyclopoida', 'Decapoda', 'Eunicida',
                       'Haplotaxida', 'Octopoda', 'Cumacea', 'Phyllodocida')) %>%
  group_by(Phylum, Class, Order, Family) %>%
  summarise(across(c(7:54), sum), .groups = 'drop') %>%
  filter(rowSums(across(where(is.numeric), ~ .x != 0)) > 0) %>%
  merge(Morphology_Family, by = c('Phylum', 'Class', 'Order', 'Family'), all = TRUE) %>%
  mutate(across(everything(), ~ replace_na(., 0)))

COI_Order <- merge(identification_COI, Decontaminant_COI, by = 'ASVNumber') %>%
  filter(round(Sim) >= 90,
         Order != "NA") %>%
  filter(!Class %in% c('Cephalopoda'),
         !Order %in% c('Decapoda', 'Littorinimorpha', 'Calanoida', 'Cumacea', 'Phyllodocida')) %>%
  group_by(Phylum, Class, Order) %>%
  summarise(across(c(8:55), sum), .groups = 'drop') %>%
  filter(rowSums(across(where(is.numeric), ~ .x != 0)) > 0) %>%
  merge(Morphology_Order, by = c('Phylum', 'Class', 'Order'), all = TRUE) %>%
  mutate(across(everything(), ~ replace_na(., 0)))

####                              Rarefaction across samples (Figure 1)      ####

# For COI
abund_table_COI <- COI_Genus %>%
  merge(Morphology, all = TRUE) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
  select(-Phylum, -Class, -Order, -Family, -Genus)

comm_matrix_COI <- as.data.frame(t(abund_table_COI))
colnames(comm_matrix_COI) <- COI_Genus$Genus

samples_COI_05 <- grep("_0.5", rownames(comm_matrix_COI), value = TRUE)
samples_COI_Max <- grep("_Max", rownames(comm_matrix_COI), value = TRUE)
samples_COI_Filt <- grep("_Filt", rownames(comm_matrix_COI), value = TRUE)
samples_COI_Morph <- grep("_Morphology", rownames(comm_matrix_COI), value = TRUE)

comm_COI_05 <- comm_matrix_COI[samples_COI_05, ]
comm_COI_Max <- comm_matrix_COI[samples_COI_Max, ]
comm_COI_Filt <- comm_matrix_COI[samples_COI_Filt, ]
comm_COI_Morph <- comm_matrix_COI[samples_COI_Morph, ]

specaccum_COI_05 <- specaccum(comm_COI_05, method = "random")
specaccum_COI_Max <- specaccum(comm_COI_Max, method = "random")
specaccum_COI_Filt <- specaccum(comm_COI_Filt, method = "random")
specaccum_COI_Morph <- specaccum(comm_COI_Morph, method = "random")

df_05_COI <- tibble(
  Samples = specaccum_COI_05$sites,
  Richness = specaccum_COI_05$richness,
  SD = specaccum_COI_05$sd,
  Method = "0.5"
)

df_Max_COI <- tibble(
  Samples = specaccum_COI_Max$sites,
  Richness = specaccum_COI_Max$richness,
  SD = specaccum_COI_Max$sd,
  Method = "Max"
)

df_Filt_COI <- tibble(
  Samples = specaccum_COI_Filt$sites,
  Richness = specaccum_COI_Filt$richness,
  SD = specaccum_COI_Filt$sd,
  Method = "Filt"
)

df_Morph_COI <- tibble(
  Samples = specaccum_COI_Morph$sites,
  Richness = specaccum_COI_Morph$richness,
  SD = specaccum_COI_Morph$sd,
  Method = "Morphology"
)

df_all_COI <- bind_rows(df_05_COI, df_Max_COI, df_Filt_COI, df_Morph_COI) %>%
  mutate(Method = replace(Method, Method == '0.5', 'Magen'),
         Method = replace(Method, Method == 'Max', 'PowerMax'),
         Method = replace(Method, Method == 'Filt', 'Filtration')

)

df_all_COI$Method <- factor(df_all_COI$Method, levels = c('Magen', 'PowerMax',
                                                  'Filtration', 'Morphology'))

Figure1A <- ggplot(df_all_COI, aes(x = Samples, y = Richness, color = Method, fill = Method)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.2, color = NA) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Number of samples",
    y = "Accumulated genera",
) +
  ggtitle("A") +
  scale_color_manual(values = c("Magen" = "blue", "PowerMax" = "red", "Filtration" = "green", 'Morphology' = 'orange')) +
  scale_fill_manual(values =  c("Magen" = "blue", "PowerMax" = "red", "Filtration" = "green", 'Morphology' = 'orange')) +
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        plot.title = element_text(
          size = 24, hjust = 0, vjust = 0, margin = margin(t = 5, r = 0, b = 7, l = 5)),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.title.x = element_text(size = 22, colour = "black"),
        axis.title.y = element_text(size = 22, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20))
Figure1A

abund_table_18S <- Decontaminant_18S

comm_matrix_18S <- as.data.frame(t(abund_table_18S[2:49]))
colnames(comm_matrix_18S) <- abund_table_18S$Phylum

samples_18S_05 <- grep("_0.5", rownames(comm_matrix_18S), value = TRUE)
samples_18S_Max <- grep("_Max", rownames(comm_matrix_18S), value = TRUE)
samples_18S_Filt <- grep("_Filt", rownames(comm_matrix_18S), value = TRUE)

comm_18S_05 <- comm_matrix_18S[samples_18S_05, ]
comm_18S_Max <- comm_matrix_18S[samples_18S_Max, ]
comm_18S_Filt <- comm_matrix_18S[samples_18S_Filt, ]

specaccum_18S_05 <- specaccum(comm_18S_05, method = "random")
specaccum_18S_Max <- specaccum(comm_18S_Max, method = "random")
specaccum_18S_Filt <- specaccum(comm_18S_Filt, method = "random")

df_18S_05 <- tibble(
  Samples = specaccum_18S_05$sites,
  Richness = specaccum_18S_05$richness,
  SD = specaccum_18S_05$sd,
  Method = "0.5"
)

df_18S_Max <- tibble(
  Samples = specaccum_18S_Max$sites,
  Richness = specaccum_18S_Max$richness,
  SD = specaccum_18S_Max$sd,
  Method = "Max"
)

df_18S_Filt <- tibble(
  Samples = specaccum_18S_Filt$sites,
  Richness = specaccum_18S_Filt$richness,
  SD = specaccum_18S_Filt$sd,
  Method = "Filt"
)

df_all_18S <- bind_rows(df_18S_05, df_18S_Max, df_18S_Filt) %>%
  mutate(Method = replace(Method, Method == '0.5', 'Magen'),
         Method = replace(Method, Method == 'Max', 'PowerMax'),
         Method = replace(Method, Method == 'Filt', 'Filtration')

)

df_all_18S$Method <- factor(df_all_18S$Method, levels = c('Magen', 'PowerMax',
                                                  'Filtration'))

Figure1B <- ggplot(df_all_18S, aes(x = Samples, y = Richness, color = Method, fill = Method)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.2, color = NA) +
  theme_minimal(base_size = 14) +
  ggtitle("B") +
  labs(
    x = "Number of samples",
    y = "Accumulated OTUs",
) +
  scale_color_manual(values = c("Magen" = "blue", "PowerMax" = "red", "Filtration" = "green", 'Morphology' = 'orange')) +
  scale_fill_manual(values =  c("Magen" = "blue", "PowerMax" = "red", "Filtration" = "green", 'Morphology' = 'orange'))+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        plot.title = element_text(
          size = 24, hjust = 0, vjust = 0, margin = margin(t = 5, r = 0, b = 7, l = 5)),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.title.x = element_text(size = 22, colour = "black"),
        axis.title.y = element_text(size = 22, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20))
Figure1B
Fig1 <- ggarrange(Figure1A, Figure1B, common.legend = TRUE, legend = "bottom")
Fig1
ggsave(Fig1, file = paste0(path_to_results, 'Fig 1.png'),
       width = 18, height = 7, dpi = 300)

####                              Richness (Figure 2)                        ####
# Richness across COI taxonomic levels
COI_Richness_ASV <- Decontaminant_COI %>%
  pivot_longer(-ASVNumber) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(value>0) %>%
  group_by(Method, Rep) %>%
  summarise(Number_of_ASVs = n()) %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g\nMagen'),
         Method = replace(Method, Method == 'Filt', '10 g\nFiltration+Magen'),
         Method = replace(Method, Method == 'Max', '10g\nPowerMax'))

COI_Richness_OTU <- OTU97_COI_df %>%
  pivot_longer(-OTU_name) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(value>0) %>%
  group_by(Method, Rep) %>%
  summarise(Number_of_OTUs = n()) %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g\nMagen'),
         Method = replace(Method, Method == 'Filt', '10 g\nFiltration+Magen'),
         Method = replace(Method, Method == 'Max', '10g\nPowerMax'))

COI_Richness_Species <- COI_Species %>%
  select(-Phylum, -Class, -Order, -Family, -Genus) %>%
  pivot_longer(-Species) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(value>0) %>%
  group_by(Method, Rep) %>%
  summarise(Number_of_Species = n()) %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g\nMagen'),
         Method = replace(Method, Method == 'Filt', '10 g\nFiltration+Magen'),
         Method = replace(Method, Method == 'Max', '10g\nPowerMax'))

COI_Richness_Genus <- COI_Genus %>%
  select(-Phylum, -Class, -Order, -Family) %>%
  pivot_longer(-Genus) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  group_by(Genus, Rep, Method) %>%
  summarise(Reads = sum(value)) %>%
  filter(Reads>0) %>%
  group_by(Method, Rep) %>%
  summarise(Number_of_Genus = n()) %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g\nMagen'),
         Method = replace(Method, Method == 'Filt', '10 g\nFiltration+Magen'),
         Method = replace(Method, Method == 'Max', '10g\nPowerMax'))

COI_Richness_Family <- COI_Family %>%
  select(-Phylum, -Class, -Order) %>%
  pivot_longer(-Family) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  group_by(Family, Rep, Method) %>%
  summarise(Reads = sum(value)) %>%
  filter(Reads>0) %>%
  group_by(Method, Rep) %>%
  summarise(Number_of_Family = n()) %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g\nMagen'),
         Method = replace(Method, Method == 'Filt', '10 g\nFiltration+Magen'),
         Method = replace(Method, Method == 'Max', '10g\nPowerMax'))

COI_Richness_Order <- COI_Order %>%
  select(-Phylum, -Class) %>%
  pivot_longer(-Order) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  group_by(Order, Rep, Method) %>%
  summarise(Reads = sum(value)) %>%
  filter(Reads>0) %>%
  group_by(Method, Rep) %>%
  summarise(Number_of_Orders = n()) %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g\nMagen'),
         Method = replace(Method, Method == 'Filt', '10 g\nFiltration+Magen'),
         Method = replace(Method, Method == 'Max', '10g\nPowerMax'))

# Merge alpha diversity across all taxonomic levels
COI_Richness <- COI_Richness_Order %>%
  left_join(COI_Richness_Family, by = c("Method", "Rep")) %>%
  left_join(COI_Richness_Genus, by = c("Method", "Rep")) %>%
  left_join(COI_Richness_Species, by = c("Method", "Rep")) %>%
  left_join(COI_Richness_OTU, by = c("Method", "Rep")) %>%
  left_join(COI_Richness_ASV, by = c("Method", "Rep")) %>%
  pivot_longer(-c(Method, Rep)) %>%
  rename(Variable = name) %>%
  mutate(value = replace(value, is.na(value), 0.00000001)) %>%
  filter(Method != 'Morphology')

COI_Richness$Method <- factor(COI_Richness$Method, levels = c('0.5 g\nMagen', '10g\nPowerMax',
                                                              '10 g\nFiltration+Magen', 'Morphology'))
Richness_statistics <- c()
for (i in (1:length(unique(COI_Richness$Variable)))){
  # i =2
  df <- COI_Richness %>%
  filter(Variable == unique(COI_Richness$Variable)[i])

  SHAPIRO <- shapiro.test(log(df[[4]]))
  LEVENE <- leveneTest(log(df[[4]])~df[[1]])
  ANOVA <- summary(aov(log(df[[4]])~df[[1]]))
  Tukey <- TukeyHSD(aov(log(df[[4]])~df[[1]]))
  # Create compact letter displays
  Diff = Tukey$`df[[1]]`[, 4] < 0.05
  Names = gsub(" ", "", rownames(Tukey$`df[[1]]`))
  names(rownames(Tukey$`df[[1]]`)) = Names
  CLD = multcompLetters(Diff)
  CLD_res <- data.frame(
    Letters = CLD$Letters) %>%
  mutate(Methods = rownames(CLD$LetterMatrix))

  Rich_values_stat <- df %>%
  group_by(Method) %>%
  summarise(
      AVG = mean(value),
      Median = median(value),
      SD = sd(value),
      SE = SD/sqrt(length(value)),
      min = min(value),
      max = max(value)) %>%
  mutate(
      pval = ANOVA[[1]]$`Pr(>F)`[1],
      pval_name='p',
      pval_name = replace(pval_name, pval<0.01, 'p < 0.01'),
      text = paste0('F[list(2, 16)] == ', round(ANOVA[[1]]$`F value`[1], 0), '*", "*',
                    " ~~ ", pval_name
))

  Rich_stat <- data.frame(
    Variable = rep(unique(COI_Richness$Variable)[i], 3),

    Shapiro_W = c(SHAPIRO$statistic, ' ', ' '),
    Shapiro_p = c(SHAPIRO$p.value, ' ', ' '),

    Levene_F = c(LEVENE$`F value`[1], ' ', ' '),
    Leene_p = c(LEVENE$`Pr(>F)`[1], ' ', ' '),

    ANOVA_F = c(round(ANOVA[[1]]$`F value`[1]), ' ', ' '),
    ANOVA_df = c(round(ANOVA[[1]]$Df[1]), ' ', ' '),
    ANOVA_pval = c(round(ANOVA[[1]]$`Pr(>F)`[1], 4), ' ', ' '),

    Tukey_Compraison_Pair = rownames(Tukey$`df[[1]]`),
    Tukey_padj = Tukey$`df[[1]]`[, 4],

    Method = c(CLD_res$Methods),
    Letters  = c(CLD_res$Letters)) %>%
  merge(Rich_values_stat, by = 'Method') %>%
  mutate(letters_y_coord = max+8,
            letters_y_coord = replace(letters_y_coord, Variable == 'Number_of_Orders' & Method == '0.5 g\nMagen', 22),
            letters_y_coord = replace(letters_y_coord, Variable == 'Number_of_Orders' & Method == '10g\nPowerMax', 22),
            letters_y_coord = replace(letters_y_coord, Variable == 'Number_of_Orders' & Method == '10 g\nFiltration+Magen', 28),

            letters_y_coord = replace(letters_y_coord, Variable == 'Number_of_Family' & Method == '10 g\nFiltration+Magen', 59),

            letters_y_coord = replace(letters_y_coord, Variable == 'Number_of_ASVs' & Method == '0.5 g\nMagen', 190),
            letters_y_coord = replace(letters_y_coord, Variable == 'Number_of_ASVs' & Method == '10g\nPowerMax', 160),
            letters_y_coord = replace(letters_y_coord, Variable == 'Number_of_ASVs' & Method == '10 g\nFiltration+Magen', 310),

            letters_y_coord = replace(letters_y_coord, Variable == 'Number_of_OTUs' & Method == '0.5 g\nMagen', 150),
            letters_y_coord = replace(letters_y_coord, Variable == 'Number_of_OTUs' & Method == '10g\nPowerMax', 140),
            letters_y_coord = replace(letters_y_coord, Variable == 'Number_of_OTUs' & Method == '10 g\nFiltration+Magen', 270),
)

  Richness_statistics <- rbind(Richness_statistics, Rich_stat)
}

COI_Richness2 <- COI_Richness %>%
  merge(Richness_statistics, by = c('Variable', 'Method'), all = TRUE) %>%
  filter(Variable != 'Number_of_Species' | Method != 'Morphology',
         Variable != 'Number_of_ASVs' | Method != 'Morphology',
         Variable != 'Number_of_OTUs' | Method != 'Morphology'
) %>%
  ungroup() %>%
  mutate(
    y_label = 0,
    y_label = replace(y_label, Variable == 'Number_of_Orders', 26),
    y_label = replace(y_label, Variable == 'Number_of_Family', 50),
    y_label = replace(y_label, Variable == 'Number_of_Genus', 45),
    y_label = replace(y_label, Variable == 'Number_of_Species', 35),
    y_label = replace(y_label, Variable == 'Number_of_OTUs', 250),
    y_label = replace(y_label, Variable == 'Number_of_ASVs', 270),

    Facet_name = 'facet_name',
    Facet_name = replace(Facet_name, Variable == 'Number_of_Orders', 'Number of orders'),
    Facet_name = replace(Facet_name, Variable == 'Number_of_Family', 'Number of families'),
    Facet_name = replace(Facet_name, Variable == 'Number_of_Genus', 'Number of genera'),
    Facet_name = replace(Facet_name, Variable == 'Number_of_Species', 'Number of species'),
    Facet_name = replace(Facet_name, Variable == 'Number_of_OTUs', 'Number of OTUs'),
    Facet_name = replace(Facet_name, Variable == 'Number_of_ASVs', 'Number of ASVs'),
    Method = as.character(Method),
    Method = replace(Method, Method == '0.5 g\nMagen', 'Magen'),
    Method = replace(Method, Method == '10g\nPowerMax', 'PowerMax'),
    Method = replace(Method, Method == '10 g\nFiltration+Magen', 'Filtration')) %>%
  filter(Method != 'Morphology')

label_data <- COI_Richness2 %>%
  distinct(Facet_name, y_label, text) %>%
  mutate(x = 1.2)

COI_Richness2$Facet_name <- factor(COI_Richness2$Facet_name, levels = c('Number of orders', 'Number of families',
                                                                        'Number of genera', 'Number of species',
                                                                        'Number of OTUs', 'Number of ASVs'))

label_data$Facet_name <- factor(label_data$Facet_name, levels = c('Number of orders', 'Number of families',
                                                                  'Number of genera', 'Number of species',
                                                                  'Number of OTUs', 'Number of ASVs'))

COI_Richness2$Method <- factor(COI_Richness2$Method, levels = c('Magen', 'PowerMax',
                                                                'Filtration'))
Figure2A <- ggplot(COI_Richness2, aes(x = Method, y = value))+
  geom_boxplot(aes(fill = Method), alpha = 0.7)+
  geom_point()+
  facet_wrap(.~Facet_name, ncol = 2,   drop = FALSE, scales = 'free_y')+
  scale_fill_manual(values = c('Magen' = 'blue',
                               'PowerMax' = 'red',
                               'Filtration' = 'green'))+
  ggtitle("A") +
  geom_label(data = label_data,
              aes(x = x, y = y_label, label = text), parse = TRUE, size = 6,
              inherit.aes = FALSE) +
  xlab(' ')+
  ylab(' ')+
  scale_x_discrete(drop = TRUE)+
  geom_text(COI_Richness2, mapping = aes(x = Method, y = letters_y_coord, label = Letters), size = 10)+
  theme_bw()+
  theme(strip.text = element_text(size = 19),
        plot.title = element_text(
          size = 24, hjust = 0, vjust = 0, margin = margin(t = 5, r = 0, b = 7, l = 5)),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.title.x = element_text(size = 22, colour = "black"),
        axis.title.y = element_text(size = 22, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20))
Figure2A
ggsave(Figure2A, file = paste0(path_to_results, 'Fig 2A.png'),
width = 13, heigh = 9, dpi = 300)

# For 18S
S18_alpha_df <- Decontaminant_18S %>%
  select(-Class, -Order, - Family, -Genus, -Species) %>%
  pivot_longer(-c(OTU_Number, Phylum)) %>%
  mutate(Phylum = replace(Phylum, Phylum!= 'Fungi', 'Unicellular eukaryotes')) %>%
  filter(value>0) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  group_by(Method, Rep, Phylum) %>%
  summarise(Number_of_OTUs = n()) %>%
  mutate(Method = replace(Method, Method == '0.5', 'Magen'),
         Method = replace(Method, Method == 'Filt', 'Filtration'),
         Method = replace(Method, Method == 'Max', 'PowerMax')) %>%
  group_by(Phylum) %>%
  mutate(Kruskal_p =  kruskal.test(Method~as.numeric(Number_of_OTUs))$p.value)

shapiro.test(subset(S18_alpha_df$Number_of_OTUs, S18_alpha_df$Phylum == 'Fungi'))
shapiro.test(subset(S18_alpha_df$Number_of_OTUs, S18_alpha_df$Phylum == 'Unicellular eukaryotes'))
hist(subset(S18_alpha_df$Number_of_OTUs, S18_alpha_df$Phylum == 'Unicellular eukaryotes'))
hist(subset(S18_alpha_df$Number_of_OTUs, S18_alpha_df$Phylum == 'Fungi'))

leveneTest(subset(S18_alpha_df$Number_of_OTUs, S18_alpha_df$Phylum == 'Unicellular eukaryotes')
           ~
             subset(S18_alpha_df$Method, S18_alpha_df$Phylum == 'Unicellular eukaryotes'))

leveneTest(subset(S18_alpha_df$Number_of_OTUs, S18_alpha_df$Phylum == 'Fungi')
           ~
             subset(S18_alpha_df$Method, S18_alpha_df$Phylum == 'Fungi'))

anova_18S_uncell <- summary(aov(subset(S18_alpha_df$Number_of_OTUs, S18_alpha_df$Phylum == 'Unicellular eukaryotes')
                                ~
                                  subset(S18_alpha_df$Method, S18_alpha_df$Phylum == 'Unicellular eukaryotes')))
anova_18S_uncell
anova_18S_fungi <- summary(aov(subset(S18_alpha_df$Number_of_OTUs, S18_alpha_df$Phylum == 'Fungi')
                               ~
                                 subset(S18_alpha_df$Method, S18_alpha_df$Phylum == 'Fungi')))
anova_18S_fungi
label_s18 <- S18_alpha_df %>%
  mutate(
    text = 'text',
    text = replace(text, Phylum == 'Fungi',

                   paste0('F[list(2, 16)] == ', round(anova_18S_fungi[[1]]$`F value`[1], 0), '*", "*',
                          " ~~ ", 'p == ',      round(anova_18S_fungi[[1]]$`Pr(>F)`[1], 3))),

    text = replace(text, Phylum == 'Unicellular eukaryotes',

                   paste0('F[list(2, 16)] == ', round(anova_18S_uncell[[1]]$`F value`[1], 0), '*", "*',
                          " ~~ ", 'p == ',      round(anova_18S_uncell[[1]]$`Pr(>F)`[1], 3))),
    x = 1.2,
    y_coord = 0,
    y_coord =ifelse(Phylum == 'Fungi', 380, 700)
) %>%
  distinct(text, x, y_coord)

S18_alpha_df$Method <- factor(S18_alpha_df$Method, levels = c('Magen', 'PowerMax', 'Filtration'))

Figure2B <- ggplot(S18_alpha_df, aes(x = Method, y = Number_of_OTUs))+
  geom_boxplot(aes(fill = Method), alpha= 0.7)+
  geom_point()+
  facet_wrap(.~Phylum, scales = 'free_y')+
  scale_fill_manual(values = c('Magen' = 'blue',
                               'PowerMax' = 'red',
                               'Filtration' = 'green'), guide = FALSE)+
  ggh4x::facetted_pos_scales(
    y = list(
      Phylum == "Unicellular eukaryotes" ~ scale_y_continuous(limits = c(0, 750)),
      Phylum == "Fungi" ~ scale_y_continuous(limits = c(0, 400))
)
)+
  geom_label(data = label_s18,
              aes(x = x, y = y_coord, label = text), parse = TRUE, size = 6,
              inherit.aes = FALSE
) +
  xlab(' ')+
  ylab('Number of OTUs ')+
  ggtitle("B") +
  # geom_text(COI_Richness2, mapping = aes(x = Method, y = letters_y_coord, label = Letters), size = 10)+
  theme_bw()+
  theme(strip.text = element_text(size = 19),
        plot.title = element_text(
          size = 24, hjust = 0, vjust = 0, margin = margin(t = 5, r = 0, b = 7, l = 5)),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.title.x = element_text(size = 22, colour = "black"),
        axis.title.y = element_text(size = 22, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20))
Figure2B
 ggsave(Figure2B, file = paste0(path_to_results, 'Fig 2B.png'),
        width = 13, heigh = 6, dpi = 300)

####                              Diversity (Figure 3)                       ####

#COI diversity

Diversity_COI <- COI_Genus %>%
  select(-Phylum, -Order, -Family) %>%
  pivot_longer(-c(Class, Genus)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  group_by(Class, Genus, Rep, Method) %>%
  summarise(Reads = sum(value)) %>%
  filter(Reads>0) %>%
  group_by(Class, Genus, Method) %>%
  summarise(Reads = sum(Reads)) %>%
  group_by(Class, Method) %>%
  summarise(Number = n()) %>%
  mutate(Method = replace(Method, Method == '0.5', 'Magen'),
         Method = replace(Method, Method == 'Filt', 'Filtration'),
         Method = replace(Method, Method == 'Max', 'PowerMax')) %>%
  filter(Method != 'Morphology')

Diversity_COI$Method <- factor(Diversity_COI$Method, levels = c('Magen', 'PowerMax', 'Filtration', 'Morphology'))
Diversity_COI$Class <- factor(Diversity_COI$Class, levels = c("Insecta", "Arachnida", "Collembola", "Chilopoda",
                                                              "Clitellata", "Polychaeta", "Gastropoda", "Malacostraca",
                                                              "Ostracoda", "Chromadorea", "Eutardigrada", "Bdelloidea"))

Fig3A <- ggplot(Diversity_COI, aes(x = Method, y = Number, fill = Class))+
  geom_bar(stat = 'identity')+
  ggtitle("A") +
  ylab('Number of genera')+
  xlab(' ')+
  theme_bw()+
  theme(strip.text = element_text(size = 19),
        plot.title = element_text(
          size = 24, hjust = 0, vjust = 0, margin = margin(t = 5, r = 0, b = 7, l = 5)),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.title.x = element_text(size = 22, colour = "black"),
        axis.title.y = element_text(size = 22, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20))
Fig3A
Diversity_COI_reads <- COI_Genus %>%
  select(-Phylum, -Order, -Family) %>%
  pivot_longer(-c(Class, Genus)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  group_by(Class, Method) %>%
  summarise(Reads = sum(value)) %>%
  filter(Reads>0) %>%
  mutate(Method = replace(Method, Method == '0.5', 'Magen'),
         Method = replace(Method, Method == 'Filt', 'Filtration'),
         Method = replace(Method, Method == 'Max', 'PowerMax')) %>%
  filter(Method != 'Morphology')

Diversity_COI_reads$Method <- factor(Diversity_COI_reads$Method, levels = c('Magen', 'PowerMax', 'Filtration', 'Morphology'))
Diversity_COI_reads$Class <- factor(Diversity_COI_reads$Class, levels = c("Insecta", "Arachnida", "Collembola", "Chilopoda",
                                                                          "Clitellata", "Polychaeta", "Gastropoda", "Malacostraca",
                                                                          "Ostracoda", "Chromadorea", "Eutardigrada", "Bdelloidea"))

Fig3B <- ggplot(Diversity_COI_reads, aes(x = Method, y = Reads, fill = Class))+
  geom_bar(stat = 'identity')+
  ylab('Number of reads')+
  ggtitle("B") +
  xlab(' ')+
  theme_bw()+
  scale_y_continuous(labels = scales::label_number(
    big.mark = " "))+
  theme(strip.text = element_text(size = 13),
        plot.title = element_text(
          size = 24, hjust = 0, vjust = 0, margin = margin(t = 5, r = 0, b = 7, l = 5)),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.title.x = element_text(size = 22, colour = "black"),
        axis.title.y = element_text(size = 22, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20))
Fig3B
Diversity_S18_df <- Decontaminant_18S %>%
  select(-Class, -Order, - Family, -Genus, -Species) %>%
  pivot_longer(-c(OTU_Number, Phylum)) %>%
  filter(value>0) %>%
  group_by(name, Phylum) %>%
  summarise(N_OTU = n()) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  mutate(
    Phylum = replace(Phylum, !Phylum %in% c('Apicomplexa', 'Cercozoa', 'Conosa', 'Fungi', 'Ciliophora',
                                            'Metamonada', 'Discoba', 'Lobosa',
                                            'Ochrophyta'), 'Other'),
    Method_full = Method,
    Method_full = replace(Method_full, Method_full == '0.5', 'Magen'),
    Method_full = replace(Method_full, Method_full == 'Max', 'PowerMax'),
    Method_full = replace(Method_full, Method_full == 'Filt', 'Filtration')) %>%
  group_by(Method, Phylum, Method_full) %>%
  summarise(n_OTU = sum(N_OTU))

Diversity_S18_df$Method_full <- factor(Diversity_S18_df$Method_full,
                                      levels = c('Magen', 'PowerMax', 'Filtration'))

Fig3C <- ggplot(Diversity_S18_df, aes(x = Method_full, y = n_OTU, fill = Phylum))+
  geom_bar(stat = 'identity')+
  ylab('Number of OTUs')+
  xlab(' ')+
  ggtitle('C')+
  theme_bw()+
  theme(strip.text = element_text(size = 19),
        plot.title = element_text(
          size = 24, hjust = 0, vjust = 0, margin = margin(t = 5, r = 0, b = 7, l = 5)),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.title.x = element_text(size = 22, colour = "black"),
        axis.title.y = element_text(size = 22, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20))
Fig3C
Diversity_S18_reads_df <- Decontaminant_18S %>%
  select(-Class, -Order, - Family, -Genus, -Species) %>%
  pivot_longer(-c(OTU_Number, Phylum)) %>%
  ungroup() %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  group_by(Method, Phylum) %>%
  summarise(N_reads = sum(value)) %>%
  mutate(
    Phylum = replace(Phylum, !Phylum %in% c('Apicomplexa', 'Cercozoa', 'Conosa', 'Fungi', 'Ciliophora',
                                            'Metamonada', 'Discoba', 'Lobosa',
                                            'Ochrophyta'), 'Other'),
    Method_full = Method,
    Method_full = replace(Method_full, Method_full == '0.5', 'Magen'),
    Method_full = replace(Method_full, Method_full == 'Max', 'PowerMax'),
    Method_full = replace(Method_full, Method_full == 'Filt', 'Filtration'))
Diversity_S18_reads_df
Diversity_S18_reads_df$Method_full <- factor(Diversity_S18_reads_df$Method_full,
                                            levels = c('Magen', 'PowerMax', 'Filtration'))

Fig3D <- ggplot(Diversity_S18_reads_df, aes(x = Method_full, y = N_reads, fill = Phylum))+
  geom_bar(stat = 'identity')+
  ylab('Number of reads')+
  xlab(' ')+
  ggtitle('D')+
  scale_y_continuous(
    labels = scales::label_number(
      big.mark = " "))+     # Use spaces as thousands separators
  theme_bw()+
  theme(strip.text = element_text(size = 19),
        plot.title = element_text(
          size = 24, hjust = 0, vjust = 0, margin = margin(t = 5, r = 0, b = 7, l = 5)),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.title.x = element_text(size = 22, colour = "black"),
        axis.title.y = element_text(size = 22, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20))
Fig3D
Fig3AB <- ggarrange(Fig3A, Fig3B, common.legend = TRUE, legend = "right")
Fig3AB
Fig3CD <- ggarrange(Fig3C, Fig3D, common.legend = TRUE, legend = "right")
Fig3CD
Fig3 <- ggarrange(Fig3AB, Fig3CD,  ncol = 1)
Fig3
ggsave(Fig3, file = paste0(path_to_results, 'Fig 3.png'),
       width = 15, height = 10, dpi = 600)

####                              NMDS (Figure 5, S3)                        ####
#For COI
Beta_COI_ASV <- Decontaminant_COI %>%
  column_to_rownames('ASVNumber') %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  separate(Sample, into = c('Sample', 'Method'), sep = '_') %>%
  filter(Method != '3')

set.seed(7)
nmds_COI_asv <- vegan::metaMDS(Beta_COI_ASV[1:(length(colnames(Beta_COI_ASV))-2)],
                           distance = 'bray', k = 3, autotransform = TRUE, na.rm = TRUE)
nmds_COI_asv
nmds_COI_sites_asv <- fortify(nmds_COI_asv) %>%
  filter(score  == 'sites') %>%
  separate(label, into = c('SAMPLE', 'Method'), sep = '_') %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g Magen'),
         Method = replace(Method, Method == 'Filt', '10 g Filtration+Magen'),
         Method = replace(Method, Method == 'Max', '10 g PowerMax'))

nmds_COI_sites_asv$Method <- factor(nmds_COI_sites_asv$Method, levels = c('0.5 g Magen', '10 g PowerMax', '10 g Filtration+Magen'))

Figure5_ASV <- ggplot(nmds_COI_sites_asv, aes(x = NMDS1, y = NMDS2, fill = Method))+
  geom_point(size = 4, shape = 21, col = 'black')+
  stat_ellipse(aes(fill = Method), geom = "polygon", alpha = 0.2, col = 'black')+
  ggtitle('ASV')+
  scale_fill_manual(values = c('0.5 g Magen' = 'blue',
                               '10 g PowerMax' = 'red',
                               '10 g Filtration+Magen' = 'green'))+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Figure5_ASV
anosim_COI_asv <- anosim(Beta_COI_ASV[1:(length(colnames(Beta_COI_ASV))-2)],
                            Beta_COI_ASV$Method, distance = "bray", permutations = 999)
anosim_COI_asv
Beta_COI_OTU <- OTU97_COI_df %>%
  column_to_rownames('OTU_name') %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  separate(Sample, into = c('Sample', 'Method'), sep = '_')

set.seed(7)
nmds_COI_otu <- vegan::metaMDS(Beta_COI_OTU[1:(length(colnames(Beta_COI_OTU))-2)],
                           distance = 'bray', k = 3, autotransform = TRUE, na.rm = TRUE)
nmds_COI_otu
Disp_COI_OTU <- betadisper(vegdist(Beta_COI_OTU[1:(length(colnames(Beta_COI_OTU))-2)]), c(rep('0.5 g Magen', 16),
                                                                                            rep('10 g PowerMax', 16),
                                                                                            rep('Filtration', 16)))
nmds_COI_sites_otu <- fortify(nmds_COI_otu) %>%
  filter(score  == 'sites') %>%
  separate(label, into = c('SAMPLE', 'Method'), sep = '_') %>%
  mutate(Method = replace(Method, Method == '0.5', 'Magen'),
         Method = replace(Method, Method == 'Filt', 'Filtration'),
         Method = replace(Method, Method == 'Max', 'PowerMax'))

nmds_COI_sites_otu$Method <- factor(nmds_COI_sites_otu$Method, levels = c('Magen', 'PowerMax', 'Filtration'))

Figure5_OTU <- ggplot(nmds_COI_sites_otu, aes(x = NMDS1, y = NMDS2, fill = Method))+
  geom_point(size = 4, shape = 21, col = 'black')+
  stat_ellipse(aes(fill = Method), geom = "polygon", alpha = 0.2, col = 'black')+
  ggtitle('OTU')+
  scale_fill_manual(values = c('Magen' = 'blue',
                               'PowerMax' = 'red',
                               'Filtration' = 'green'))+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Figure5_OTU
anosim_COI_otu <- anosim(Beta_COI_OTU[1:(length(colnames(Beta_COI_OTU))-2)],
                            Beta_COI_OTU$Method, distance = "bray", permutations = 999)
anosim_COI_otu
Beta_COI_Species <- COI_Species %>%
  filter(Species != 'NA') %>%
  mutate(Taxa = paste0(Species, ' (', Family, ')')) %>%
  select(-Phylum, - Class, - Order, - Family, -Genus, - Species) %>%
  column_to_rownames('Taxa') %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  separate(Sample, into = c('Sample', 'Method'), sep = '_') %>%
  filter(Method != 'Morphology') %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g Magen'),
         Method = replace(Method, Method == 'Filt', '10 g Filtration+Magen'),
         Method = replace(Method, Method == 'Max', '10 g PowerMax')) %>%
  filter(rowSums(across(where(is.numeric), ~ .x != 0)) > 0)

set.seed(7)
nmds_COI_species <- vegan::metaMDS(Beta_COI_Species[1:(length(colnames(Beta_COI_Species))-2)],
                               distance = 'bray', k = 3, autotransform = TRUE, na.rm = TRUE)
nmds_COI_species
anosim_COI_species <- anosim(Beta_COI_Species[1:(length(colnames(Beta_COI_Species))-2)],
                         Beta_COI_Species$Method, distance = "bray", permutations = 999)
anosim_COI_species
Disp_COI_species <- betadisper(vegdist(Beta_COI_Species[1:(length(colnames(Beta_COI_Species))-2)]), c(rep('0.5 g Magen', 16),
                                                                                                        rep('10 g PowerMax', 16),
                                                                                                    rep('Filtration', 16)))

nmds_COI_sites_species <- fortify(nmds_COI_species) %>%
  filter(score  == 'sites') %>%
  separate(label, into = c('SAMPLE', 'Method'), sep = '_') %>%
  mutate(Method = replace(Method, Method == '0.5', 'Magen'),
         Method = replace(Method, Method == 'Filt', 'Filtration'),
         Method = replace(Method, Method == 'Max', 'PowerMax'))

nmds_COI_sites_species$Method <- factor(nmds_COI_sites_species$Method, levels = c('Magen', 'PowerMax', 'Filtration'))

Figure5_Species <- ggplot(nmds_COI_sites_species, aes(x = NMDS1, y = NMDS2, fill = Method))+
  geom_point(size = 4, shape = 21, col = 'black')+
  stat_ellipse(aes(fill = Method), geom = "polygon", alpha = 0.2, col = 'black')+
  ggtitle('Species')+
  scale_fill_manual(values = c('Magen' = 'blue',
                               'PowerMax' = 'red',
                               'Filtration' = 'green'))+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Figure5_Species
Beta_COI_Genus <- COI_Genus %>%
  filter(Genus != 'NA') %>%
  mutate(Taxa = paste0(Genus, ' (', Family, ')')) %>%
  select(-Phylum, - Class, - Order, - Family, - Genus) %>%
  column_to_rownames('Taxa') %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  separate(Sample, into = c('Sample', 'Method'), sep = '_') %>%
  filter(Method != 'Morphology') %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g Magen'),
         Method = replace(Method, Method == 'Filt', '10 g Filtration+Magen'),
         Method = replace(Method, Method == 'Max', '10 g PowerMax'))

anosim_COI_genus <- anosim(Beta_COI_Genus[1:(length(colnames(Beta_COI_Genus))-2)],
                       Beta_COI_Genus$Method, distance = "bray", permutations = 999)
anosim_COI_genus
set.seed(7)
nmds_COI_genus <- vegan::metaMDS(Beta_COI_Genus[1:(length(colnames(Beta_COI_Genus))-2)],
                             distance = 'bray', k = 3, autotransform = TRUE, na.rm = TRUE)
nmds_COI_genus
Disp_COI_genus <- betadisper(vegdist(Beta_COI_Genus[1:(length(colnames(Beta_COI_Genus))-2)]), c(rep('0.5 g Magen', 16),
                                                                                                      rep('10 g PowerMax', 16),
                                                                                                        rep('Filtration', 16)))
# Test for multivariate homogeneity of dispersion
anova(Disp_COI_genus)
# Pairwise comparisons
TukeyHSD(Disp_COI_genus)

nmds_COI_sites_genus <- fortify(nmds_COI_genus) %>%
  filter(score  == 'sites') %>%
  separate(label, into = c('SAMPLE', 'Method'), sep = '_') %>%
  mutate(Method = replace(Method, Method == '0.5', 'Magen'),
         Method = replace(Method, Method == 'Filt', 'Filtration'),
         Method = replace(Method, Method == 'Max', 'PowerMax'))

nmds_COI_sites_genus$Method <- factor(nmds_COI_sites_genus$Method, levels = c('Magen', 'PowerMax', 'Filtration'))

Figure5_Genus <- ggplot(nmds_COI_sites_genus, aes(x = NMDS1, y = NMDS2, fill = Method))+
  geom_point(size = 4, shape = 21, col = 'black')+
  stat_ellipse(aes(fill = Method), geom = "polygon", col = 'black', alpha = 0.2)+
  scale_fill_manual(values = c('Magen' = 'blue',
                               'PowerMax' = 'red',
                               'Filtration' = 'green'))+
  ggtitle("Genus") +
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Figure5_Genus
Beta_COI_Family <- COI_Family %>%
  filter(Family != 'NA') %>%
  select(-Phylum, - Class, - Order) %>%
  group_by(Family) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames('Family') %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  separate(Sample, into = c('Sample', 'Method'), sep = '_') %>%
  filter(Method != 'Morphology') %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g Magen'),
         Method = replace(Method, Method == 'Filt', '10 g Filtration+Magen'),
         Method = replace(Method, Method == 'Max', '10 g PowerMax'))

anosim_COI_family <- anosim(Beta_COI_Family[1:(length(colnames(Beta_COI_Family))-2)],
                        Beta_COI_Family$Method, distance = "jaccard", permutations = 999)
anosim_COI_family
set.seed(7)
nmds_COI_family <- vegan::metaMDS(Beta_COI_Family[1:(length(colnames(Beta_COI_Family))-2)],
                              distance = 'bray', k = 3, autotransform = TRUE, na.rm = TRUE)
nmds_COI_family
nmds_COI_sites_family <- fortify(nmds_COI_family) %>%
  filter(score  == 'sites') %>%
  separate(label, into = c('SAMPLE', 'Method'), sep = '_') %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g Magen'),
         Method = replace(Method, Method == 'Filt', '10 g Filtration+Magen'),
         Method = replace(Method, Method == 'Max', '10 g PowerMax'))

nmds_COI_sites_family$Method <- factor(nmds_COI_sites_family$Method, levels = c('0.5 g Magen', '10 g PowerMax', '10 g Filtration+Magen'))

Figure5_Family <- ggplot(nmds_COI_sites_family, aes(x = NMDS1, y = NMDS2, fill = Method))+
  geom_point(size = 4, shape = 21, col = 'black')+
  stat_ellipse(aes(fill = Method), geom = "polygon", col = 'black', alpha = 0.2)+
  scale_fill_manual(values = c('0.5 g Magen' = 'blue',
                               '10 g PowerMax' = 'red',
                               '10 g Filtration+Magen' = 'green'))+
  ggtitle("Family") +
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Figure5_Family
Beta_COI_Order <- COI_Order %>%
  # filter(Class %in% c('Clitellata', 'Arachnida', 'Chilopoda', 'Collembola', 'Insecta', 'Malacostraca'),
  #        Order != 'Enchytraeida') %>%
  filter(Order != 'NA') %>%
  select(-Phylum, - Class) %>%
  group_by(Order) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames('Order') %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  separate(Sample, into = c('Sample', 'Method'), sep = '_') %>%
  filter(Method != 'Morphology') %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g Magen'),
         Method = replace(Method, Method == 'Filt', '10 g Filtration+Magen'),
         Method = replace(Method, Method == 'Max', '10 g PowerMax'))

anosim_COI_order <- anosim(Beta_COI_Order[1:(length(colnames(Beta_COI_Order))-2)],
                       Beta_COI_Order$Method, distance = "jaccard", permutations = 999)
anosim_COI_order
set.seed(7)
nmds_COI_order <- vegan::metaMDS(Beta_COI_Order[1:(length(colnames(Beta_COI_Order))-2)],
                             distance = 'bray', k = 3, autotransform = TRUE, na.rm = TRUE)
nmds_COI_order
nmds_COI_sites_order <- fortify(nmds_COI_order) %>%
  filter(score  == 'sites') %>%
  separate(label, into = c('SAMPLE', 'Method'), sep = '_') %>%
  mutate(Method = replace(Method, Method == '0.5', 'Magen'),
         Method = replace(Method, Method == 'Filt', 'Filtration'),
         Method = replace(Method, Method == 'Max', 'PowerMax'))

nmds_COI_sites_order$Method <- factor(nmds_COI_sites_order$Method, levels = c('Magen', 'PowerMax', 'Filtration'))

Figure5_order <- ggplot(nmds_COI_sites_order, aes(x = NMDS1, y = NMDS2, fill = Method))+
  geom_point(size = 4, shape = 21, col = 'black')+
  stat_ellipse(aes(fill = Method), geom = "polygon", col = 'black', alpha = 0.2)+
  scale_fill_manual(values = c('Magen' = 'blue',
                               'PowerMax' = 'red',
                               'Filtration' = 'green'))+
  ggtitle("Order") +
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Figure5_order
Figure5 <- ggarrange(Figure5_order, Figure5_Family, Figure5_Genus,
                            Figure5_Species, Figure5_OTU, Figure5_ASV,
                            common.legend = TRUE, legend = "right")
Figure5
 ggsave(Figure5, file = paste0(path_to_results, '/Fig 5.png'),
 width = 15, height = 10, dpi = 300)

#For18S

Beta_18S_otu <- Decontaminant_18S %>%
  select(-c(Phylum, Class, Order, Family, Genus, Species)) %>%
  pivot_longer(
    cols = ends_with("_0.5") | ends_with("_Max") | ends_with("_Filt"),
    names_to = "Sample",
    values_to = "Abundance") %>%
  pivot_wider(
    names_from = OTU_Number,
    values_from = 'Abundance',
    values_fill = list(Abundance = 0)) %>%
  separate(Sample, into = c('Rep', 'Method'), sep = '_') %>%
  mutate(Sample = paste(Rep, Method, sep = '_')) %>%
  column_to_rownames("Sample")

set.seed(7)
nmds_18S_OTU <- vegan::metaMDS(Beta_18S_otu[3:length(colnames(Beta_18S_otu))],
                               distance = 'bray', k = 3, autotransform = TRUE, na.rm = TRUE)
nmds_18S_OTU
nmds_sites_18S_OTU <- fortify(nmds_18S_OTU) %>%
  filter(score  == 'sites') %>%
  separate(label, into = c('SAMPLE', 'Method'), sep = '_') %>%
  mutate(Method = replace(Method, Method == '0.5', 'Magen'),
         Method = replace(Method, Method == 'Filt', 'Filtration'),
         Method = replace(Method, Method == 'Max', 'PowerMax'))

nmds_sites_18S_OTU$Method <- factor(nmds_sites_18S_OTU$Method, levels = c('Magen', 'PowerMax',
                                                                          'Filtration'))
Figure_S3 <- ggplot(nmds_sites_18S_OTU, aes(x = NMDS1, y = NMDS2, fill = Method))+
  geom_point(size = 4, shape = 21)+
  stat_ellipse(aes(fill = Method), col = 'black', geom = "polygon", alpha = 0.2)+
  scale_fill_manual(values = c('Magen' = 'blue',
                               'PowerMax' = 'red',
                               'Filtration' = 'green'))+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20))
Figure_S3
 ggsave(Figure_S3, file = paste0(path_to_results, 'Fig S3.png'),
 width = 9, height = 7, dpi = 300)

anosim_18S_otu <- anosim(Beta_18S_otu[3:(length(colnames(Beta_18S_otu)))],
                         Beta_18S_otu$Method, distance = "bray", permutations = 999)
anosim_18S_otu
# Total family richness of different orders
#Only orders which were identify both morphology and filtration

Figure6_df <- COI_Family %>%
  select(-Phylum) %>%
  pivot_longer(-c(Class, Order, Family)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(Method %in% c('Morphology', 'Filt')) %>%
  mutate(value = replace(value, value != 0, 1)) %>%
  group_by(Order, Family, Method) %>%
  filter(value>0) %>%
  summarise(Sample_Number = n()) %>%
  group_by(Order, Family, Method) %>%
  summarise(SumRep = sum(Sample_Number)) %>%
  group_by(Order, Method) %>%
  summarise(Fam_number = n()) %>%
  mutate(Method = replace(Method, Method == 'Filt', 'Filtration')
) %>%
  filter(Order %in% c('Sarcoptiformes', 'Mesostigmata', 'Trombidiformes', 'Entomobryomorpha',
                      'Neelipleona', 'Poduromorpha', 'Symphypleona', 'Araneae', 'Coleoptera', 'Diptera',
                      'Hemiptera', 'Crassiclitellata', 'Geophilomorpha', 'Lithobiomorpha'))

Figure6_df$Order <- factor(Figure6_df$Order, levels = c('Lithobiomorpha', 'Geophilomorpha', 'Crassiclitellata',
                                                                              'Hemiptera', 'Diptera', 'Coleoptera', 'Aranea', 'Symphypleona',
                                                                              'Poduromorpha', 'Neelipleona',  'Entomobryomorpha', 'Araneae',
                                                                              'Trombidiformes', 'Mesostigmata', 'Sarcoptiformes'))

Figure6 <- ggplot(Figure6_df, aes(x = Fam_number, y = Order, fill = Method))+
  geom_bar(stat="identity", position="dodge", colour="black",
           width = 0.7)+
  scale_fill_manual(values = c('Morphology' = 'orange',
                               'Filtration' = 'green'))+
  xlab('Number of familes')+
  ylab('Orders')+
  theme_bw()+
  theme(strip.text = element_text(size = 17),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, colour = "black", angle=0), #, hjust=1, vjust=1
        axis.title.x = element_text(size = 24, colour = "black"),
        axis.title.y = element_text(size = 24, colour = "black"),
        legend.title=element_text(face="bold", size=22),
        legend.text=element_text(face="bold", size=20),
        legend.position = "bottom")
Figure6
ggsave(Figure6, file = paste0(path_to_results, 'Fig 6.png'),
       width = 12, height = 10, dpi = 600)

####                              Venn diagram (Figure 4, 7)                 ####
Venn_COI_Species <- COI_Species %>%
  pivot_longer(-c(Phylum, Class, Order, Family, Genus, Species)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(value>0)

Venn_COI_Species_list <- list(
  `0.5` = unique(subset(Venn_COI_Species, Method == '0.5')$Species),
  Filt= unique(subset(Venn_COI_Species, Method == 'Filt')$Species),
  Max = unique(subset(Venn_COI_Species, Method == 'Max')$Species))

Venn_species_graph <- ggVennDiagram(
  Venn_COI_Species_list,
  label_size = 6,
  category.names = c(" ", " ", " "),
  label_alpha = 0.1,
  set_size = 6,
  label = "both",
  category.label.geom = "text",
  set_color = c("blue", "green", "red")) +
  scale_fill_gradient(low = "white", high = "darkolivegreen1") +
  ggtitle('Species')+
  geom_label(
    data = data.frame(
      x = c(-3.2, 7, 2.1),
      y = c(4.7, 4.7, -8.65),
      label = c("Magen", "Filtration", "PowerMax")),
    aes(x = x, y = y, label = label),
    fill = c("blue", "green", "red"), alpha = 0.3, color = "black", size = 8)+
  theme(strip.text = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.title = element_text(hjust = 0.5, vjust = -1, size = 26, face = "bold"))
Venn_species_graph
Venn_COI_data_Genus <- COI_Genus %>%
  pivot_longer(-c(Phylum, Class, Order, Family, Genus)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(value>0)

Venn_Genuses_list <- list(
  `0.5` = unique(subset(Venn_COI_data_Genus, Method == '0.5')$Genus),
  Filtration = unique(subset(Venn_COI_data_Genus, Method == 'Filt')$Genus),
  PowerMax = unique(subset(Venn_COI_data_Genus, Method == 'Max')$Genus))

Venn_genus_graph <- ggVennDiagram(
  Venn_Genuses_list,
  label_size = 6,
  category.names = c(" ", " ", " "),
  label_alpha = 0.1,
  set_size = 6,
  label = "both",
  category.label.geom = "text",
  set_color = c("blue", "green", "red")) +
  scale_fill_gradient(low = "white", high = "darkolivegreen1") +
  ggtitle("Genera") +
  geom_label(
    data = data.frame(
      x = c(-3.2, 7, 2.1),
      y = c(4.7, 4.7, -8.65),
      label = c("Magen", "Filtration", "PowerMax")),
    aes(x = x, y = y, label = label),
    fill = c("blue", "green", "red"), alpha = 0.3, color = "black", size = 8)+
  theme(strip.text = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.title = element_text(hjust = 0.5, vjust = -1, size = 26, face = "bold"))
Venn_genus_graph

Venn_COI_data_Family <- COI_Family %>%
  pivot_longer(-c(Phylum, Class, Order, Family)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(value>0)

Venn_Families_list <- list(
  `0.5` = unique(subset(Venn_COI_data_Family, Method == '0.5')$Family),
  Filtration= unique(subset(Venn_COI_data_Family, Method == 'Filt')$Family),
  PowerMax = unique(subset(Venn_COI_data_Family, Method == 'Max')$Family)

)

Venn_family_graph <- ggVennDiagram(
  Venn_Families_list,
  label_size = 6,
  category.names = c(" ", " ", " "),
  label_alpha = 0.1,
  set_size = 6,
  label = "both",
  category.label.geom = "text",
  set_color = c("blue", "green", "red")) +
  scale_fill_gradient(low = "white", high = "darkolivegreen1") +
  ggtitle("Families") +
  geom_label(
    data = data.frame(
      x = c(-3.2, 7, 2.1),
      y = c(4.7, 4.7, -8.65),
      label = c("Magen", "Filtration", "PowerMax")),
    aes(x = x, y = y, label = label),
    fill = c("blue", "green", "red"), alpha = 0.3, color = "black", size = 8)+
  theme(strip.text = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.title = element_text(hjust = 0.5, vjust = -1, size = 26, face = "bold"))
Venn_family_graph
Venn_COI_data_Orders <- COI_Order %>%
  pivot_longer(-c(Phylum, Class, Order)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(value>0)

Venn_Orders_list <- list(
  `0.5` = unique(subset(Venn_COI_data_Orders, Method == '0.5')$Order),
  Filtration= unique(subset(Venn_COI_data_Orders, Method == 'Filt')$Order),
  PowerMax = unique(subset(Venn_COI_data_Orders, Method == 'Max')$Order))

Venn_orders_graph <- ggVennDiagram(
  Venn_Orders_list,
  label_size = 6,
  category.names = c(" ", " ", " "),
  label_alpha = 0.1,
  set_size = 6,
  label = "both",
  set_color = c("blue", "green", "red")) +
  scale_fill_gradient(low = "white", high = "darkolivegreen1") +
  ggtitle("Orders") +
  geom_label(
    data = data.frame(
      x = c(-3.2, 7, 2.1),
      y = c(4.7, 4.7, -8.65),
      label = c("Magen", "Filtration", "PowerMax")),
    aes(x = x, y = y, label = label),
    fill = c("blue", "green", "red"), alpha = 0.3, color = "black", size = 8)+
  theme(strip.text = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.title = element_text(hjust = 0.5, vjust = -1, size = 26, face = "bold"))
Venn_orders_graph
Figure4 <- ggarrange(Venn_orders_graph, Venn_family_graph, Venn_genus_graph, Venn_species_graph,
                  common.legend = TRUE, legend = "none")
Figure4
 ggsave(Figure4, file = paste0(path_to_results, 'Fig 4.png'),
        width = 13, height = 15, bg = "white", dpi = 300)

### Comparison between filtration and morphology

Venn_data_Genus_Morph <- Venn_COI_data_Genus %>%
  filter(Class %in% c('Clitellata', 'Arachnida', 'Chilopoda', 'Collembola', 'Insecta', 'Malacostraca', 'Gastropoda'),
         Order != 'Enchytraeida')

Venn_Genuses_list2 <- list(
  `Morphology` = unique(subset(Venn_data_Genus_Morph, Method == 'Morphology')$Genus),
  Filtration = unique(subset(Venn_data_Genus_Morph, Method == 'Filt')$Genus))

Venn_genus_graph2 <- ggVennDiagram(
  Venn_Genuses_list2, label_size = 6,
  category.names = c(" ", " "),
  set_size = 6,
  label_alpha = 0.1,
  label_fill = "white",
  label = "both",  # or "percent"
  set_color = c('orange', "green")) +
  scale_fill_gradient(low = "white", high = "darkolivegreen1")+  # ggplot2 style
  ggtitle("Genera") +
  scale_x_continuous(expand = expansion(mult = .3)) +  # Expand horizontal space for labels
  geom_label(aes(x = c(5.8, 6.2), y = c(4, -0.5),  # Label coordinates
                label = c("Filtration", "Morphology")),
            size = 6,
            fill = c("green", 'orange'),
            alpha = 0.3,
            fontface = "bold",
            color = "black")+
  theme(strip.text = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour = "white"),  # White background
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.title = element_text(
          hjust = 0.45,        # Horizontal centering
          vjust = -1,         # Move the title upward (negative values shift up)
          size = 20,
          face = "bold"),
        legend.position = "none")
Venn_genus_graph2
Venn_data_Family_Morph <- Venn_COI_data_Family %>%
  filter(Class %in% c('Clitellata', 'Arachnida', 'Chilopoda', 'Collembola', 'Insecta', 'Malacostraca', 'Gastropoda'),
         Order != 'Enchytraeida')

Venn_Families_list2 <- list(
  `Morphology` = unique(subset(Venn_data_Family_Morph, Method == 'Morphology')$Family),
  Filtration= unique(subset(Venn_data_Family_Morph, Method == 'Filt')$Family))

library(ggVennDiagram)
Venn_family_graph2 <- ggVennDiagram(
  Venn_Families_list2, label_size = 6,
  category.names = c(" ", " "),
  set_size = 6,
  label_alpha = 0.1,
  label_fill = "white",
  label = "both",  # or "percent"
  set_color = c('orange', "green")) +
  scale_fill_gradient(low = "white", high = "darkolivegreen1")+  # ggplot2 style
  ggtitle("Families") +
  scale_x_continuous(expand = expansion(mult = .3)) +  # Expand horizontal space for labels
  geom_label(aes(x = c(5.8, 6.2), y = c(4, -0.5),  # Label coordinates
                 label = c("Filtration", "Morphology")),
             size = 6,
             fill = c("green", 'orange'),
             alpha = 0.3,
             fontface = "bold",
             color = "black")+
  theme(strip.text = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour = "white"),  # White background
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.title = element_text(
          hjust = 0.45,        # Horizontal centering
          vjust = -1,         # Move the title upward (negative values shift up)
          size = 20,
          face = "bold"),
        legend.position = "none")
Venn_family_graph2
Venn_data_Order_Morph <- Venn_COI_data_Orders %>%
  filter(Class %in% c('Clitellata', 'Arachnida', 'Chilopoda', 'Collembola', 'Insecta', 'Malacostraca', 'Gastropoda'),
         Order != 'Enchytraeida')

Venn_Orders_list2 <- list(
  `Morphology` = unique(subset(Venn_data_Order_Morph, Method == 'Morphology')$Order),
  Filtration= unique(subset(Venn_data_Order_Morph, Method == 'Filt')$Order))

Venn_orders_graph2 <- ggVennDiagram(
  Venn_Orders_list2,
  label_size = 6,
  category.names = c(" ", " "),
  set_size = 6,
  label_alpha = 0.1,
  label_fill = "white",
  label = "both",  # or "percent"
  set_color = c('orange', "green")) +
  scale_fill_gradient(low = "white", high = "darkolivegreen1")+  # ggplot2 style
  ggtitle("Orders") +
  scale_x_continuous(expand = expansion(mult = .3)) +  # Expand horizontal space for labels
  geom_label(aes(x = c(5.8, 6.2), y = c(4, -0.5),  # Label coordinates
                 label = c("Filtration", "Morphology")),
             size = 6,
             fill = c("green", 'orange'),
             alpha = 0.3,
             fontface = "bold",
             color = "black")+
  theme(strip.text = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour = "white"),  # White background
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.title = element_text(
          hjust = 0.45,        # Horizontal centering
          vjust = -1,         # Move the title upward (negative values shift up)
          size = 20,
          face = "bold"),
        legend.position = "none")
Venn_orders_graph2
Figure7 <- ggarrange(Venn_orders_graph2, Venn_family_graph2, Venn_genus_graph2, common.legend = TRUE,
                   legend = "none", nrow = 1)
Figure7
ggsave(Figure7, file = paste0(path_to_results, 'Fig 7.png'),
       width = 20, height = 15, bg = 'white', dpi = 300)

####                              Quantitative metrics Family (Figure 8)     ####
Quant_Morphology_Family <- Morphology_Family %>%
  mutate(group = Class,
         group = replace(group, Order %in% c('Mesostigmata'), 'Mesostigmata'),
         group = replace(group, Order %in% c('Sarcoptiformes'), 'Sarcoptiformes'),
         group = replace(group, Order %in% c('Trombidiformes'), 'Trombidiformes'),
         group = replace(group, Order %in% c('Araneae'), 'Araneae')
)

Quant_COI_family <- COI_Family %>%
  mutate(group = Class,
         group = replace(group, Order %in% c('Mesostigmata'), 'Mesostigmata'),
         group = replace(group, Order %in% c('Sarcoptiformes'), 'Sarcoptiformes'),
         group = replace(group, Order %in% c('Trombidiformes'), 'Trombidiformes'),
         group = replace(group, Order %in% c('Araneae'), 'Araneae')) %>%
  filter(group %in% c(unique(Quant_Morphology_Family$group)))

ABAND <- Quant_Morphology_Family %>%
  pivot_longer(-c(Phylum, Class, Order, Family, group)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  group_by(Family, group, Method) %>%
  summarise(Abundance = mean(value)) %>%
  mutate(Variable = 'Aband') %>%
  select(group, Family, Variable, Abundance) %>%
  dcast(group+Family~Variable, value.var = 'Abundance')

FOO <- Quant_COI_family %>%
  pivot_longer(-c(Phylum, Class, Order, Family, group)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(value>0,
         Method != 'Morphology') %>%
  group_by(group, Family, Method) %>%
  summarise(
    Occurence = n(),
    FOO = Occurence*100/16) %>%
  mutate(Variable = paste('FOO', Method, sep = '_')) %>%
  select(group, Family, Variable, FOO) %>%
  dcast(group+Family~Variable, value.var = 'FOO') %>%
  mutate_all(~ifelse(is.na(.), 0, .))

RRA <- Quant_COI_family %>%
  pivot_longer(-c(Phylum, Class, Order, Family, group)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(value>0,
         Method != 'Morphology') %>%
  group_by(Rep, Method) %>%
  mutate(
    total_reads_in_sample = sum(value, na.rm = TRUE),
    prop_in_sample = value / total_reads_in_sample
) %>%
  group_by(Family, group, Method) %>%
  summarise(
    RRA = mean(prop_in_sample, na.rm = TRUE) * 100,
    .groups = 'drop'
) %>%
  mutate(Variable = paste('RRA', Method, sep = '_')) %>%
  select(group, Family, Variable, RRA) %>%
  dcast(group+Family~Variable, value.var = 'RRA') %>%
  mutate_all(~ifelse(is.na(.), 0, .))

raw_metrics_family <- merge(RRA, FOO, by = c('group', 'Family')) %>%
  merge(ABAND, by = c('group', 'Family'), all = TRUE) %>%
  filter(!is.na(Aband)) %>%
  filter(!Family %in% c('Enchytraeidae',
                        'Lumbricidae',
                        'Geophilidae', 'Lithobiidae')) %>%
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  mutate(size_group = 'Macro',
         size_group = replace(size_group, group %in% c('Mesostigmata', 'Sarcoptiformes',
                                                       'Trombidiformes', 'Collembola'), 'Meso')) %>%
  pivot_longer(-c(group, size_group, Family)) %>%
  separate(name, into = c('Methric', 'Method'), sep = '_') %>%
  filter(Methric %in% c('FOO', 'RRA', 'Aband')) %>%
  tidyr::unite(Variable, Methric, Method, na.rm = TRUE) %>%
  dcast(group+size_group+Family~Variable, value.var = 'value') %>%
  pivot_longer(-c(group, size_group, Family, Aband)) %>%
  mutate(facet_name = 'l',
         facet_name = replace(facet_name, name == 'FOO_0.5', 'FOO 0.5 g Magen'),
         facet_name = replace(facet_name, name == 'FOO_Max', 'FOO 10 g PowerMax'),
         facet_name = replace(facet_name, name == 'FOO_Filt', 'FOO 10 g Filtration + Magen'),

         facet_name = replace(facet_name, name == 'RRA_0.5', 'RRA 0.5 g Magen'),
         facet_name = replace(facet_name, name == 'RRA_Max', 'RRA 10 g PowerMax'),
         facet_name = replace(facet_name, name == 'RRA_Filt', 'RRA 10 g Filtration + Magen')) %>%
  group_by(group, facet_name, name) %>%
  filter(group != 'Araneae') %>%
  
  
  mutate(
    Spearman_coef = cor.test(Aband, value, method = 'spearman')$estimate,
    Spearman_pval = cor.test(Aband, value, method = 'spearman')$p.value,
    Spearman_text = paste0('R = ', round(Spearman_coef,2), ', p = ', round(Spearman_pval,2)),
    Spearman_line = ifelse(Spearman_pval < 0.05, 'solid',
                           ifelse(Spearman_pval < 0.1, 'dashed',
                                  NA)),
    show_se = case_when(
      Spearman_pval < 0.05 ~ 0.5, TRUE ~ 0),
    
    x_coord_group = 0,
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'RRA_0.5', 0.015),
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'FOO_0.5', 6),
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'FOO_Max', 15),
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'RRA_Max', 2),
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'RRA_Filt', 1.3),
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'FOO_Filt', 60),

    #
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'FOO_0.5', 30),
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'RRA_0.5', 1.3),
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'FOO_Max', 20),
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'RRA_Max', 1),
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'RRA_Filt', 0.6),
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'FOO_Filt', 60),

    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'FOO_0.5', 75),
    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'RRA_0.5', 0.75),
    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'FOO_Max', 70),
    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'RRA_Max', 0.65),
    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'RRA_Filt', 1),
    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'FOO_Filt', 110),

    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'FOO_0.5', 35),
    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'RRA_0.5', 7.5),
    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'FOO_Max', 100),
    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'RRA_Max', 1.7),
    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'RRA_Filt', 3.5),
    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'FOO_Filt', 120),

    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'FOO_0.5', 17),
    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'RRA_0.5', 0.3),
    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'FOO_Max', 17),
    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'RRA_Max', 6),
    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'FOO_Filt', 45),
    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'RRA_Filt', 1.5)

)

raw_metrics_family$facet_name <- factor(raw_metrics_family$facet_name,
                                         levels = c('FOO 0.5 g Magen', 'FOO 10 g PowerMax', 'FOO 10 g Filtration + Magen',
                                                    'RRA 0.5 g Magen', 'RRA 10 g PowerMax', 'RRA 10 g Filtration + Magen'))

facet_labels <- c(
  "FOO 0.5 g Magen" = "Magen",
  "FOO 10 g PowerMax" = "PowerMax",
  "FOO 10 g Filtration + Magen" = "Filtration",
  "RRA 0.5 g Magen" = "Magen",
  "RRA 10 g PowerMax" = "PowerMax",
  "RRA 10 g Filtration + Magen" = "Filtration"
)

Quant_collembola_family <- ggplot(subset(raw_metrics_family, group == 'Collembola'),
                                  aes(y = subset(raw_metrics_family$value, raw_metrics_family$group == 'Collembola'),
                                      x = subset(raw_metrics_family$Aband, raw_metrics_family$group == 'Collembola'))) +
  geom_point(size = 3) +
  facet_wrap2(~ subset(raw_metrics_family$facet_name, raw_metrics_family$group == 'Collembola'), scales = "free_y", nrow = 2,
              labeller = as_labeller(facet_labels)) +
  coord_cartesian(ylim = c(0, NA))+

  geom_smooth(aes(color = facet_name,
                   linetype = name,
                   alpha = show_se),
              method = "lm", se = TRUE,
              size = 2) +
  scale_alpha_identity() +
  scale_color_manual(values = c("blue", "red", "green",  "blue", "red", "green"), guide = FALSE)+
  #
  scale_linetype_manual(values = setNames(subset(raw_metrics_family$Spearman_line, raw_metrics_family$group == 'Collembola'),
                                          unique(raw_metrics_family$name)),
                        guide = "none")+

  geom_label(
    data = subset(raw_metrics_family, raw_metrics_family$group == 'Collembola'),
    aes(x = 1500, y = subset(raw_metrics_family$x_coord_group, raw_metrics_family$group == 'Collembola'),
        label = subset(raw_metrics_family$Spearman_text, raw_metrics_family$group == 'Collembola')),
    color = "black", size = 5
)+
  ggtitle('Collembola families')+
  labs(tag = "B") +
  ylab('RRA, %               FOO, %')+
  xlab('Abundance ex./m²')+
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        plot.tag = element_text(
          size = 18, face = "bold", hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Quant_collembola_family
Quant_mesostigmata_family <- ggplot(subset(raw_metrics_family, group == 'Mesostigmata'),
                                    aes(y = subset(raw_metrics_family$value, raw_metrics_family$group == 'Mesostigmata'),
                                        x = subset(raw_metrics_family$Aband, raw_metrics_family$group == 'Mesostigmata'))) +
  geom_point(size = 3) +
  facet_wrap(~ subset(raw_metrics_family$facet_name, raw_metrics_family$group == 'Mesostigmata'), scales = "free_y", nrow = 2,
             labeller = as_labeller(facet_labels)) +
  coord_cartesian(ylim = c(0, NA))+

  labs(tag = "C") +
  geom_smooth(aes(color = facet_name,
                   linetype = name,
                   alpha = show_se),
              method = "lm", se = TRUE,
              size = 2) +
  scale_alpha_identity() +
  scale_color_manual(values = c("blue", "red", "green",  "blue", "red", "green"), guide = FALSE)+
  #
  scale_linetype_manual(values = setNames(subset(raw_metrics_family$Spearman_line, raw_metrics_family$group == 'Mesostigmata'),
                                          unique(raw_metrics_family$name)),
                        guide = "none")+

  geom_label(
    data = subset(raw_metrics_family, raw_metrics_family$group == 'Mesostigmata'),
    aes(x = 2000, y = subset(raw_metrics_family$x_coord_group, raw_metrics_family$group == 'Mesostigmata'),
        label = subset(raw_metrics_family$Spearman_text, raw_metrics_family$group == 'Mesostigmata')),
    color = "black", size = 5
)+

  ggtitle('Mesostigmata families')+
  ylab('RRA, %               FOO, %')+
  xlab('Abundance ex./m²')+
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        plot.tag = element_text(
          size = 18, face = "bold", hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Quant_mesostigmata_family
Quant_trombidiformes_family <- ggplot(subset(raw_metrics_family, group == 'Trombidiformes'),
                                      aes(y = subset(raw_metrics_family$value, raw_metrics_family$group == 'Trombidiformes'),
                                          x = subset(raw_metrics_family$Aband, raw_metrics_family$group == 'Trombidiformes'))) +
  geom_point(size = 3) +
  facet_wrap(~ subset(raw_metrics_family$facet_name, raw_metrics_family$group == 'Trombidiformes'), scales = "free_y", nrow = 2,
             labeller = as_labeller(facet_labels)) +
  coord_cartesian(ylim = c(0, NA))+
  labs(tag = "E") +
  geom_smooth(aes(color = facet_name,
                   linetype = name,
                   alpha = show_se),
              method = "lm", se = TRUE,
              size = 2) +
  scale_alpha_identity() +
  scale_color_manual(values = c("blue", "red", "green",  "blue", "red", "green"), guide = FALSE)+
  scale_linetype_manual(values = setNames(subset(raw_metrics_family$Spearman_line,
                                                 raw_metrics_family$group == 'Trombidiformes'),
                                          unique(raw_metrics_family$name)),
                        guide = "none")+

  geom_label(
    data = subset(raw_metrics_family, raw_metrics_family$group == 'Trombidiformes'),
    aes(x = 1000, y = subset(raw_metrics_family$x_coord_group, raw_metrics_family$group == 'Trombidiformes'),
        label = subset(raw_metrics_family$Spearman_text, raw_metrics_family$group == 'Trombidiformes')),
    color = "black", size = 5
)+

  ggtitle('Trombidiformes families')+
  ylab('RRA, %               FOO, %')+
  xlab('Abundance ex./m²')+
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        plot.tag = element_text(
          size = 18, face = "bold", hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Quant_trombidiformes_family
Quant_sarcoptiformes_family <- ggplot(subset(raw_metrics_family, group == 'Sarcoptiformes'),
                                      aes(y = subset(raw_metrics_family$value, raw_metrics_family$group == 'Sarcoptiformes'),
                                          x = subset(raw_metrics_family$Aband, raw_metrics_family$group == 'Sarcoptiformes'))) +
  geom_point(size = 3) +
  facet_wrap(~ subset(raw_metrics_family$facet_name, raw_metrics_family$group == 'Sarcoptiformes'), scales = "free_y", nrow = 2,
             labeller = as_labeller(facet_labels)) +
  coord_cartesian(ylim = c(0, NA))+
  labs(tag = "D") +

  geom_smooth(aes(color = facet_name,
                   linetype = name,
                   alpha = show_se),
              method = "lm", se = TRUE,
              size = 2) +
  scale_alpha_identity() +

  scale_color_manual(values = c("blue", "red", "green",  "blue", "red", "green"), guide = FALSE)+
  #
  scale_linetype_manual(values = setNames(subset(raw_metrics_family$Spearman_line, raw_metrics_family$group == 'Sarcoptiformes'),
                                          unique(raw_metrics_family$name)),
                        guide = "none")+

  geom_label(
    data = subset(raw_metrics_family, raw_metrics_family$group == 'Sarcoptiformes'),
    aes(x = 1000, y = subset(raw_metrics_family$x_coord_group, raw_metrics_family$group == 'Sarcoptiformes'),
        label = subset(raw_metrics_family$Spearman_text, raw_metrics_family$group == 'Sarcoptiformes')),
    color = "black", size = 5
)+

  ggtitle('Sarcoptiformes families')+
  ylab('RRA, %               FOO, %')+
  xlab('Abundance ex./m²')+
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        plot.tag = element_text(
          size = 18, face = "bold", hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold",
        margin = margin(b = 20)))
Quant_sarcoptiformes_family
Quant_Insecta_family <- ggplot(subset(raw_metrics_family, group == 'Insecta'),
                               aes(y = subset(raw_metrics_family$value, raw_metrics_family$group == 'Insecta'),
                                   x = subset(raw_metrics_family$Aband, raw_metrics_family$group == 'Insecta'))) +
  geom_point(size = 3) +
  facet_wrap(~ subset(raw_metrics_family$facet_name, raw_metrics_family$group == 'Insecta'), scales = "free_y", nrow = 2,
             labeller = as_labeller(facet_labels)) +

  geom_smooth(aes(color = facet_name,
                   linetype = name,
                   alpha = show_se),
              method = "lm", se = TRUE,
              size = 2) +
  scale_alpha_identity() +

  scale_color_manual(values = c("blue", "red", "green",  "blue", "red", "green"), guide = FALSE)+
  #
  scale_linetype_manual(values = setNames(subset(raw_metrics_family$Spearman_line, raw_metrics_family$group == 'Insecta'),
                                          unique(raw_metrics_family$name)),
                        guide = "none")+

  labs(tag = "A") +

  geom_label(
    data = subset(raw_metrics_family, raw_metrics_family$group == 'Insecta'),
    aes(x = 150, y = subset(raw_metrics_family$x_coord_group, raw_metrics_family$group == 'Insecta'),
        label = subset(raw_metrics_family$Spearman_text, raw_metrics_family$group == 'Insecta')),
    color = "black", size = 5
)+

  ggtitle('Insecta families')+
  ylab('RRA, %               FOO, %')+
  xlab('Abundance ex./m²')+
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        plot.tag = element_text(
          size = 18, face = "bold", hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Quant_Insecta_family
quant_arrange <- ggarrange(Quant_Insecta_family, Quant_collembola_family, Quant_mesostigmata_family,
                           Quant_sarcoptiformes_family, Quant_trombidiformes_family,
                           common.legend = TRUE, ncol = 1)
quant_arrange
ggsave(quant_arrange, file = paste0(path_to_results, 'Fig 8.png'),
       width = 14, height = 25, dpi = 500)

####                              Quantitative metrics Genera (Figure S4)    ####
ABAND <- Morphology_Genus %>%
  mutate(group = Class,
         group = replace(group, Order %in% c('Mesostigmata'), 'Mesostigmata'),
         group = replace(group, Order %in% c('Sarcoptiformes'), 'Sarcoptiformes'),
         group = replace(group, Order %in% c('Trombidiformes'), 'Trombidiformes'),
         group = replace(group, Order %in% c('Araneae'), 'Araneae')) %>%
  pivot_longer(-c(Phylum, Class, Order, Family, Genus, group)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  mutate(Family = replace(Family, Family == 'Geophylidae', 'Geophilidae')) %>%
  group_by(Genus, group, Method) %>%
  summarise(Abundance = mean(value)) %>%
  mutate(Variable = 'Aband') %>%
  select(group, Genus, Variable, Abundance) %>%
  dcast(group+Genus~Variable, value.var = 'Abundance')

FOO <- COI_Genus %>%
  mutate(group = Class,
         group = replace(group, Order %in% c('Mesostigmata'), 'Mesostigmata'),
         group = replace(group, Order %in% c('Sarcoptiformes'), 'Sarcoptiformes'),
         group = replace(group, Order %in% c('Trombidiformes'), 'Trombidiformes'),
         group = replace(group, Order %in% c('Araneae'), 'Araneae')) %>%
  pivot_longer(-c(Phylum, Class, Order, Family, Genus, group)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(value>0,
         Method != 'Morphology') %>%
  group_by(group, Genus, Method) %>%
  summarise(
    Occurence = n(),
    FOO = Occurence*100/16) %>%
  mutate(Variable = paste('FOO', Method, sep = '_')) %>%
  select(group, Genus, Variable, FOO) %>%
  dcast(group+Genus~Variable, value.var = 'FOO') %>%
  mutate_all(~ifelse(is.na(.), 0, .))

RRA <- COI_Genus %>%
  mutate(group = Class,
         group = replace(group, Order %in% c('Mesostigmata'), 'Mesostigmata'),
         group = replace(group, Order %in% c('Sarcoptiformes'), 'Sarcoptiformes'),
         group = replace(group, Order %in% c('Trombidiformes'), 'Trombidiformes'),
         group = replace(group, Order %in% c('Araneae'), 'Araneae')
) %>%
  pivot_longer(-c(Phylum, Class, Order, Family, Genus, group)) %>%
  separate(name, into = c('Rep', 'Method'), sep = '_') %>%
  filter(value>0,
         Method != 'Morphology') %>%
  group_by(Rep, Method) %>%
  mutate(
    total_reads_in_sample = sum(value, na.rm = TRUE),
    prop_in_sample = value / total_reads_in_sample
) %>%
  group_by(Genus, group, Method) %>%
  summarise(
    RRA = mean(prop_in_sample, na.rm = TRUE) * 100,
    .groups = 'drop'
) %>%
  mutate(Variable = paste('RRA', Method, sep = '_')) %>%
  select(group, Genus, Variable, RRA) %>%
  dcast(group+Genus~Variable, value.var = 'RRA') %>%
  mutate_all(~ifelse(is.na(.), 0, .))

raw_metrics_genus <- merge(RRA, FOO, by = c('group', 'Genus')) %>%
  merge(ABAND, by = c('group', 'Genus'), all = TRUE) %>%
  filter(!is.na(Aband)) %>%
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  mutate(size_group = 'Macro',
         size_group = replace(size_group, group %in% c('Mesostigmata', 'Sarcoptiformes',
                                                       'Trombidiformes', 'Collembola'), 'Meso')) %>%
  pivot_longer(-c(group, size_group, Genus)) %>%
  separate(name, into = c('Methric', 'Method'), sep = '_') %>%
  filter(Methric %in% c('FOO', 'RRA', 'Aband')) %>%
  tidyr::unite(Variable, Methric, Method, na.rm = TRUE) %>%
  dcast(group+size_group+Genus~Variable, value.var = 'value') %>%
  pivot_longer(-c(group, size_group, Genus, Aband)) %>%
  mutate(facet_name = 'l',
         facet_name = replace(facet_name, name == 'FOO_0.5', 'FOO 0.5 g Magen'),
         facet_name = replace(facet_name, name == 'FOO_Max', 'FOO 10 g PowerMax'),
         facet_name = replace(facet_name, name == 'FOO_Filt', 'FOO 10 g Filtration + Magen'),

         facet_name = replace(facet_name, name == 'RRA_0.5', 'RRA 0.5 g Magen'),
         facet_name = replace(facet_name, name == 'RRA_Max', 'RRA 10 g PowerMax'),
         facet_name = replace(facet_name, name == 'RRA_Filt', 'RRA 10 g Filtration + Magen')) %>%
  group_by(group, facet_name, name) %>%
  filter(group != 'Araneae') %>%
  mutate(
    
    # Sperman correlation
    Spearman_coef = cor.test(Aband, value, method = 'spearman')$estimate,
    Spearman_pval = cor.test(Aband, value, method = 'spearman')$p.value,
    Spearman_text = paste0('R = ', round(Spearman_coef,2), ', p = ', round(Spearman_pval,2)),
    Spearman_line = ifelse(Spearman_pval < 0.05, 'solid',
                           ifelse(Spearman_pval < 0.1, 'dashed',
                                  NA)),
    show_se = case_when(
      Spearman_pval < 0.05 ~ 0.5, TRUE ~ 0),

    x_coord_group = 0,
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'RRA_0.5', 0.005),
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'FOO_0.5', 0.005),
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'FOO_Max', 3),
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'RRA_Max', 5),
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'RRA_Filt', 2.5),
    x_coord_group = replace(x_coord_group, group == 'Mesostigmata' & name == 'FOO_Filt', 50),

    #
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'FOO_0.5', 7),
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'RRA_0.5', 50),
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'FOO_Max', 7),
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'RRA_Max', 0.50),
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'RRA_Filt', 1.5),
    x_coord_group = replace(x_coord_group, group == 'Sarcoptiformes' & name == 'FOO_Filt', 60),

    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'FOO_0.5', 7),
    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'RRA_0.5', 0.75),
    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'FOO_Max', 7),
    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'RRA_Max', 0.4),
    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'RRA_Filt', 0.25),
    x_coord_group = replace(x_coord_group, group == 'Trombidiformes' & name == 'FOO_Filt', 30),

    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'FOO_0.5', 8),
    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'RRA_0.5', 15),
    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'FOO_Max', 50),
    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'RRA_Max', 3),
    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'RRA_Filt', 2.5),
    x_coord_group = replace(x_coord_group, group == 'Collembola' & name == 'FOO_Filt', 72),

    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'FOO_0.5', 17),
    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'RRA_0.5', 0.3),
    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'FOO_Max', 17),
    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'RRA_Max', 6),
    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'FOO_Filt', 45),
    x_coord_group = replace(x_coord_group, group == 'Insecta' & name == 'RRA_Filt', 1.5)

)

raw_metrics_genus$facet_name <- factor(raw_metrics_genus$facet_name,
                                         levels = c('FOO 0.5 g Magen', 'FOO 10 g PowerMax', 'FOO 10 g Filtration + Magen',
                                                    'RRA 0.5 g Magen', 'RRA 10 g PowerMax', 'RRA 10 g Filtration + Magen'))

facet_labels <- c(
  "FOO 0.5 g Magen" = "Magen",
  "FOO 10 g PowerMax" = "PowerMax",
  "FOO 10 g Filtration + Magen" = "Filtration",
  "RRA 0.5 g Magen" = "Magen",
  "RRA 10 g PowerMax" = "PowerMax",
  "RRA 10 g Filtration + Magen" = "Filtration"
)

Quant_collembola_genus <- ggplot(subset(raw_metrics_genus, group == 'Collembola'),
                                  aes(y = subset(raw_metrics_genus$value, raw_metrics_genus$group == 'Collembola'),
                                      x = subset(raw_metrics_genus$Aband, raw_metrics_genus$group == 'Collembola'))) +
  geom_point(size = 3) +
  facet_wrap2(~ subset(raw_metrics_genus$facet_name, raw_metrics_genus$group == 'Collembola'), scales = "free_y", nrow = 2,
              labeller = as_labeller(facet_labels)) +
  coord_cartesian(ylim = c(0, NA))+

  geom_smooth(aes(color = facet_name,
                   linetype = name,
                   alpha = show_se),
              method = "lm", se = TRUE,
              size = 2) +
  scale_alpha_identity() +
  scale_color_manual(values = c("blue", "red", "green",  "blue", "red", "green"), guide = FALSE)+
  #
  scale_linetype_manual(values = setNames(subset(raw_metrics_genus$Spearman_line, raw_metrics_genus$group == 'Collembola'),
                                          unique(raw_metrics_genus$name)),
                        guide = "none")+

  geom_label(
    data = subset(raw_metrics_genus, raw_metrics_genus$group == 'Collembola'),
    aes(x = 1500, y = subset(raw_metrics_genus$x_coord_group, raw_metrics_genus$group == 'Collembola'),
        label = subset(raw_metrics_genus$Spearman_text, raw_metrics_genus$group == 'Collembola')),
    color = "black", size = 5
)+
  ggtitle('Collembola genera')+
  labs(tag = "A") +
  ylab('RRA, %               FOO, %')+
  xlab('Abundance ex./m²')+
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        plot.tag = element_text(
          size = 18, face = "bold", hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Quant_collembola_genus
Quant_mesostigmata_genus <- ggplot(subset(raw_metrics_genus, group == 'Mesostigmata'),
                                    aes(y = subset(raw_metrics_genus$value, raw_metrics_genus$group == 'Mesostigmata'),
                                        x = subset(raw_metrics_genus$Aband, raw_metrics_genus$group == 'Mesostigmata'))) +
  geom_point(size = 3) +
  facet_wrap(~ subset(raw_metrics_genus$facet_name, raw_metrics_genus$group == 'Mesostigmata'), scales = "free_y", nrow = 2,
             labeller = as_labeller(facet_labels)) +
  coord_cartesian(ylim = c(0, NA))+

  labs(tag = "B") +

  geom_smooth(aes(color = facet_name,
                   linetype = name,
                   alpha = show_se),
              method = "lm", se = TRUE,
              size = 2) +
  scale_alpha_identity() +
  scale_color_manual(values = c("blue", "red", "green",  "blue", "red", "green"), guide = FALSE)+
  #
  scale_linetype_manual(values = setNames(subset(raw_metrics_genus$Spearman_line, raw_metrics_genus$group == 'Mesostigmata'),
                                          unique(raw_metrics_genus$name)),
                        guide = "none")+

  geom_label(
    data = subset(raw_metrics_genus, group == 'Mesostigmata'),
    aes(x = 2000, y = x_coord_group,
        label = Spearman_text),
     size = 5
)+

  ggtitle('Mesostigmata genera')+
  ylab('RRA, %               FOO, %')+
  xlab('Abundance ex./m²')+
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        plot.tag = element_text(
          size = 18, face = "bold", hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Quant_mesostigmata_genus
Quant_trombidiformes_genus <- ggplot(subset(raw_metrics_genus, group == 'Trombidiformes'),
                                      aes(y = subset(raw_metrics_genus$value, raw_metrics_genus$group == 'Trombidiformes'),
                                          x = subset(raw_metrics_genus$Aband, raw_metrics_genus$group == 'Trombidiformes'))) +
  geom_point(size = 3) +
  facet_wrap(~ subset(raw_metrics_genus$facet_name, raw_metrics_genus$group == 'Trombidiformes'), scales = "free_y", nrow = 2,
             labeller = as_labeller(facet_labels)) +
  coord_cartesian(ylim = c(0, NA))+
  labs(tag = "D") +

  geom_smooth(aes(color = facet_name,
                   linetype = name,
                   alpha = show_se),
              method = "lm", se = TRUE,
              size = 2) +
  scale_alpha_identity() +
  scale_color_manual(values = c("blue", "red", "green",  "blue", "red", "green"), guide = FALSE)+
  #
  scale_linetype_manual(values = setNames(subset(raw_metrics_genus$Spearman_line, raw_metrics_genus$group == 'Trombidiformes'),
                                          unique(raw_metrics_genus$name)),
                        guide = "none")+

  geom_label(
    data = subset(raw_metrics_genus, raw_metrics_genus$group == 'Trombidiformes'),
    aes(x = 1000, y = subset(raw_metrics_genus$x_coord_group, raw_metrics_genus$group == 'Trombidiformes'),
        label = subset(raw_metrics_genus$Spearman_text, raw_metrics_genus$group == 'Trombidiformes')),
    color = "black", size = 5
)+

  ggtitle('Trombidiformes genera')+
  ylab('RRA, %               FOO, %')+
  xlab('Abundance ex./m²')+
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        plot.tag = element_text(
          size = 18, face = "bold", hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Quant_trombidiformes_genus
Quant_sarcoptiformes_genus <- ggplot(subset(raw_metrics_genus, group == 'Sarcoptiformes'),
                                      aes(y = subset(raw_metrics_genus$value, raw_metrics_genus$group == 'Sarcoptiformes'),
                                          x = subset(raw_metrics_genus$Aband, raw_metrics_genus$group == 'Sarcoptiformes'))) +
  geom_point(size = 3) +
  facet_wrap(~ subset(raw_metrics_genus$facet_name, raw_metrics_genus$group == 'Sarcoptiformes'), scales = "free_y", nrow = 2,
             labeller = as_labeller(facet_labels)) +
  coord_cartesian(ylim = c(0, NA))+
  labs(tag = "C") +

  geom_smooth(aes(color = facet_name,
                   linetype = name,
                   alpha = show_se),
              method = "lm", se = TRUE,
              size = 2) +
  scale_alpha_identity() +
  scale_color_manual(values = c("blue", "red", "green",  "blue", "red", "green"), guide = FALSE)+
  scale_linetype_manual(values = setNames(subset(raw_metrics_genus$Spearman_line, raw_metrics_genus$group == 'Sarcoptiformes'),
                                          unique(raw_metrics_genus$name)),
                        guide = "none")+

  geom_label(
    data = subset(raw_metrics_genus, raw_metrics_genus$group == 'Sarcoptiformes'),
    aes(x = 2000, y = subset(raw_metrics_genus$x_coord_group, raw_metrics_genus$group == 'Sarcoptiformes'),
        label = subset(raw_metrics_genus$Spearman_text, raw_metrics_genus$group == 'Sarcoptiformes')),
    color = "black", size = 5
)+

  ggtitle('Sarcoptiformes genera')+
  ylab('RRA, %               FOO, %')+
  xlab('Abundance ex./m²')+
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        plot.tag = element_text(
          size = 18, face = "bold", hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 17, colour = "black"),
        legend.title=element_text(face="bold", size=20),
        legend.text=element_text(face="bold", size=20),
        plot.title = element_text(
          size = 20, hjust = 0.5, face = "bold",
          margin = margin(b = 20)))
Quant_sarcoptiformes_genus
quant_arrange_genus <- ggarrange(Quant_collembola_genus, Quant_mesostigmata_genus,
                           Quant_sarcoptiformes_genus, Quant_trombidiformes_genus, ncol = 1)
quant_arrange_genus
   ggsave(quant_arrange_genus, file = paste0(path_to_results, 'Fig S4.png'),
     width = 12, height = 20, dpi = 500)

####                              Rarefaction curves (Figure S1)             ####

  rarecurve_out_COI <- rarecurve(counts_t_COI,
                            step = 1000,
                            sample = min(rowSums(counts_t_COI)),  tidy = TRUE)

  rarefaction_ggplot_df_COI <- rarecurve_out_COI %>%
  mutate(Method = sub(".*_", "", Site)) %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g Magen'),
           Method = replace(Method, Method == 'Filt', '10 g Filtration'),
           Method = replace(Method, Method == 'Max', '10 g PowerMax'))

  rarefaction_ggplot_df_COI$Method <- factor(rarefaction_ggplot_df_COI$Method, levels = c('0.5 g Magen',
                                                                                '10 g PowerMax',
                                                                                '10 g Filtration'))

  rarefaction_COI_graph <- ggplot(rarefaction_ggplot_df_COI, aes(x = Sample, y = Species, group = Site, color = Method)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = min(rowSums(counts_t_COI)), color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = "A",
      x = "Number of reads",
      y = "ASV richness"
) +

    geom_text(x = 500000, y = 500, label = paste0('Rarefaction depth ', min(rowSums(counts_t_COI))), size = 7, col = 'black')+
    scale_color_manual(
      values = c('0.5 g Magen' = "blue",  '10 g PowerMax' = "red",  '10 g Filtration'= "green"),
      name = "Method"
) +
    scale_x_continuous(
      labels = scales::number)+
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0, size = 25, face = "bold"),
      legend.position = c(0.85, 0.70),
      legend.background = element_rect(fill = "white", color = "grey80"),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 15, colour = "black"),
      axis.text.x = element_text(size = 15, colour = "black"),
      axis.title.x = element_text(size = 17, colour = "black"),
      axis.title.y = element_text(size = 17, colour = "black"),
      legend.title=element_text(face="bold", size=20),
      legend.text=element_text(face="bold", size=20)
)
rarefaction_COI_graph

  rarecurve_out_18S <- rarecurve(counts_t_18S,
                            step = 1000,
                            sample = min(rowSums(counts_t_18S)),  tidy = TRUE)

  # Create a table for ggplot
  rarefaction_ggplot_df_18S <- rarecurve_out_18S %>%
  mutate(Method = sub(".*_", "", Site)) %>%
  mutate(Method = replace(Method, Method == '0.5', '0.5 g Magen'),
           Method = replace(Method, Method == 'Filt', '10 g Filtration'),
           Method = replace(Method, Method == 'Max', '10 g PowerMax'))

  rarefaction_ggplot_df_18S$Method <- factor(rarefaction_ggplot_df_18S$Method, levels = c('0.5 g Magen',
                                                                                '10 g PowerMax',
                                                                                '10 g Filtration'))

  # Plot using ggplot2
  rarefaction_18S_graph <- ggplot(rarefaction_ggplot_df_18S, aes(x = Sample, y = Species, group = Site, color = Method)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = sample_depth_18S, color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = "B",
      x = "Number of reads",
      y = "ASV richness"
) +
    scale_x_continuous(
      labels = scales::number)+
    geom_text(x = 250000, y = 500, label = paste0('Rarefaction depth ', min(rowSums(counts_t_18S))), size = 7, col = 'black')+
    scale_color_manual(
      values = c('0.5 g Magen' = "blue",  '10 g PowerMax' = "red",  '10 g Filtration'= "green"),
      name = "Method"
) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0, size = 25, face = "bold"),
      legend.position = c(0.85, 0.70),
      legend.background = element_rect(fill = "white", color = "grey80"),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 15, colour = "black"),
      axis.text.x = element_text(size = 15, colour = "black"),
      axis.title.x = element_text(size = 17, colour = "black"),
      axis.title.y = element_text(size = 17, colour = "black"),
      legend.title=element_text(face="bold", size=20),
      legend.text=element_text(face="bold", size=20)
)
rarefaction_18S_graph
  FigureS1 <- ggarrange(rarefaction_COI_graph, rarefaction_18S_graph, ncol = 2)
FigureS1
  ggsave(FigureS1, file = paste0(path_to_results, 'Fig S1.png'),
         width = 20, heigh = 7, dpi = 300)

####                              DNA yield (Figure S2)                      ####
  Conc_df <- Concentration %>%
  separate(sample_id, into = c('REP', 'Method'), sep = '_') %>%
  filter(!REP %in% c('K', 'NC', 'PO')) %>%
  mutate(Method_Full =  'Method',
           Method_Full = replace(Method_Full, Method == '0.5', 'Magen'),
           Method_Full = replace(Method_Full, Method == 'Max', 'PowerMax'),
           Method_Full = replace(Method_Full, Method == 'Filt', 'Filtration'),

           text = 'a',
           text = replace(text, Method_Full == 'PowerMax', 'b'),
           text = replace(text, Method_Full == 'Filtration', 'c'),
           y_coord = 127,
           y_coord = replace(y_coord, Method_Full == 'PowerMax', 75),
           y_coord = replace(y_coord, Method_Full == 'Filtration', 120))

  Conc_df$Method_Full <- factor(Conc_df$Method_Full, level = c('Magen', 'PowerMax',
                                                          'Filtration'))

  FigS2 <- ggplot(Conc_df, aes(x = Method_Full, y = Conc, fill = Method_Full))+
    geom_boxplot(alpha = 0.7)+
    geom_point()+
    scale_fill_manual(values = c('Magen' = 'blue',
                                 'PowerMax' = 'red',
                                 'Filtration' = 'green'), guide = FALSE)+
    xlab(' ')+
    ylab('DNA yield µg/ml')+
    # ggtitle()+
    # scale_x_discrete(drop = TRUE)+
    geom_text(Conc_df, mapping = aes(x = Method_Full, y = y_coord, label = text), size = 10)+
    theme_bw()+
    theme(strip.text = element_text(size = 17),
          axis.text.y = element_text(size = 15, colour = "black"),
          axis.text.x = element_text(size = 19, colour = "black"),
          axis.title.x = element_text(size = 19, colour = "black"),
          axis.title.y = element_text(size = 19, colour = "black"),
          legend.title=element_text(face="bold", size=20),
          legend.text=element_text(face="bold", size=20))
FigS2
  ggsave(FigS2, file = paste0(path_to_results, 'Fig S2.png'),
         width = 9, heigh = 7, dpi = 300)

####                              Dispersion analysis                        ####

vegdist_COI_ASV <- vegdist(Beta_COI_ASV[1:(length(colnames(Beta_COI_ASV))-2)])
vegdist_COI_OTU <- vegdist(Beta_COI_ASV[1:(length(colnames(Beta_COI_OTU))-2)])
vegdist_COI_Species <- vegdist(Beta_COI_ASV[1:(length(colnames(Beta_COI_Species))-2)])
vegdist_COI_Genus <- vegdist(Beta_COI_ASV[1:(length(colnames(Beta_COI_Genus))-2)])
vegdist_COI_Family <- vegdist(Beta_COI_ASV[1:(length(colnames(Beta_COI_Family))-2)])
vegdist_COI_Order <- vegdist(Beta_COI_ASV[1:(length(colnames(Beta_COI_Order))-2)])
vegdist_18S_OTU <- vegdist(Beta_18S_otu[3:length(colnames(Beta_18S_otu))])

Disp_COI_ASV <- betadisper(vegdist_COI_ASV, c(rep('0.5 g Magen', 16), rep('10 g PowerMax', 16), rep('Filtration', 16)))
Disp_COI_OTU <- betadisper(vegdist_COI_OTU, c(rep('0.5 g Magen', 16), rep('10 g PowerMax', 16), rep('Filtration', 16)))
disp_COI_Species <- betadisper(vegdist_COI_Species, c(rep('0.5 g Magen', 16), rep('10 g PowerMax', 16), rep('Filtration', 16)))
Disp_COI_Genus <- betadisper(vegdist_COI_Genus, c(rep('0.5 g Magen', 16), rep('10 g PowerMax', 16), rep('Filtration', 16)))
Disp_COI_Family <- betadisper(vegdist_COI_Family, c(rep('0.5 g Magen', 16), rep('10 g PowerMax', 16), rep('Filtration', 16)))
Disp_COI_Order <- betadisper(vegdist_COI_Order, c(rep('0.5 g Magen', 16), rep('10 g PowerMax', 16), rep('Filtration', 16)))
Disp_18S_OTU <- betadisper(vegdist_18S_OTU, c(rep('0.5 g Magen', 16), rep('10 g PowerMax', 16), rep('Filtration', 16)))

aov_COI_ASV_disp <- anova(Disp_COI_ASV)
aov_COI_OTU_disp <- anova(Disp_COI_OTU)
aov_COI_Species_disp <- anova(disp_COI_Species)
aov_COI_Genus_disp <- anova(Disp_COI_Genus)
aov_COI_Family_disp <- anova(Disp_COI_Family)
aov_COI_Order_disp <- anova(Disp_COI_Order)
aov_18S_OTU_disp <- anova(Disp_18S_OTU)

####                              Table S2                                   ####

Table_S2 <- data.frame(
  Rank = c('Order', 'Family', 'Genus', 'Species', '97% OTU', 'ASV', '18S 99% OTU'),
  NMDS_stress = c(round(nmds_COI_order$stress, 2),
                  round(nmds_COI_family$stress, 2),
                  round(nmds_COI_genus$stress, 2),
                  round(nmds_COI_species$stress, 2),
                  round(nmds_COI_otu$stress, 2),
                  round(nmds_COI_asv$stress, 2),
                  round(nmds_18S_OTU$stress, 2)),
  NMDS_Dimensions =c(3, 3, 3, 3, 3, 3, 3),

  ANOSIM_R = c(round(anosim_COI_order$statistic, 2),
               round(anosim_COI_family$statistic, 2),
               round(anosim_COI_genus$statistic, 2),
               round(anosim_COI_species$statistic, 2),
               round(anosim_COI_otu$statistic, 2),
               round(anosim_COI_asv$statistic, 2),
               round(anosim_18S_otu$statistic, 2)),

  ANOSIM_Significance = c(anosim_COI_order$signif,
                          anosim_COI_family$signif,
                          anosim_COI_genus$signif,
                          anosim_COI_species$signif,
                          anosim_COI_otu$signif,
                          anosim_COI_asv$signif,
                          anosim_18S_otu$signif),
  BetaDisperRes = c(
    paste0('F = ', round(aov_COI_Order_disp$`F value`[1]), ', df = ', round(aov_COI_Order_disp$Df[1]), ', p = ', round(aov_COI_Order_disp$`Pr(>F)`[1], 3)),
    paste0('F = ', round(aov_COI_Family_disp$`F value`[1]), ', df = ', round(aov_COI_Family_disp$Df[1]), ', p = ', round(aov_COI_Family_disp$`Pr(>F)`[1], 3)),
    paste0('F = ', round(aov_COI_Genus_disp$`F value`[1]), ', df = ', round(aov_COI_Genus_disp$Df[1]), ', p = ', round(aov_COI_Genus_disp$`Pr(>F)`[1], 3)),
    paste0('F = ', round(aov_COI_Species_disp$`F value`[1]), ', df = ', round(aov_COI_Species_disp$Df[1]), ', p = ', round(aov_COI_Species_disp$`Pr(>F)`[1], 3)),
    paste0('F = ', round(aov_COI_OTU_disp$`F value`[1]), ', df = ', round(aov_COI_OTU_disp$Df[1]), ', p = ', round(aov_COI_OTU_disp$`Pr(>F)`[1], 3)),
    paste0('F = ', round(aov_COI_ASV_disp$`F value`[1]), ', df = ', round(aov_COI_ASV_disp$Df[1]), ', p = ', round(aov_COI_ASV_disp$`Pr(>F)`[1], 3)),
    paste0('F = ', round(aov_18S_OTU_disp$`F value`[1]), ', df = ', round(aov_18S_OTU_disp$Df[1]), ', p = ', round(aov_18S_OTU_disp$`Pr(>F)`[1], 3))

))

write_xlsx(Table_S2, path = 'Results/eDNA/Extract_Test/Worker_graphs/TableS2.xlsx')
