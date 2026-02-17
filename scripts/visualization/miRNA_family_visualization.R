library(ape)
library(ggplot2)
library(dplyr)
library(cowplot)


mir36_pred <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/scripts/results/predictions/filtered_gff/mir36.gff")
wormbasemir35 <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/wormbasemir35family.tsv", col_names = c("seqid","start","end","family","strand")) %>%
  dplyr::mutate(length = end - start)


ggplot() +
  geom_rect(data = mir36_pred, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1, ymax = 2, fill = "MirMachine_predictions")) + 
  geom_rect(data = wormbasemir35, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 2.01, ymax = 3.01, fill = "WormBase_mir35_family")) +
  facet_wrap(~seqid)+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.x = element_line(),
    strip.text = element_text(size = 16, color = "black")) +
  xlab("N2 Genome Coordinates (Mb)")


locus1 <- ggplot() +
  geom_rect(data = mir36_pred, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1, ymax = 2, fill = "MirMachine_predictions")) + 
  geom_rect(data = wormbasemir35, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 2.01, ymax = 3.01, fill = "WormBase_mir35_family")) +
  facet_wrap(~seqid)+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    # axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.x = element_line(),
    strip.text = element_text(size = 16, color = "black")) +
  # xlab("N2 Genome Coordinates (Mb)") +
  coord_cartesian(xlim = c(11.537,11.539))
locus1

locus2 <- ggplot() +
  geom_rect(data = mir36_pred, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1, ymax = 2, fill = "MirMachine_predictions")) + 
  geom_rect(data = wormbasemir35, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 2.01, ymax = 3.01, fill = "WormBase_mir35_family")) +
  facet_wrap(~seqid)+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.x = element_line(),
    strip.text = element_text(size = 16, color = "black")) +
  xlab("N2 Genome Coordinates (Mb)") +
  coord_cartesian(xlim = c(11.8875,11.8925))
locus2



concatplot <- cowplot::plot_grid(
    locus1, locus2,
    rel_heights = c(1, 1),
    ncol = 1,
    align = "v",
    axis = "lr")
concatplot



briggsae <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/briggsae/results/predictions/filtered_gff/QX1410.PRE.gff") %>% 
  dplyr::mutate(length = end-start) %>%
  dplyr::filter(grepl("Mir-36",attributes))

II_left <- ggplot() +
  geom_rect(data = briggsae %>% dplyr::filter(seqid == "II"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1, ymax = 2), fill = "#53886C") + 
  geom_text(data = briggsae %>% dplyr::filter(seqid == "II"), aes(x = ((start + end) / 2) / 1e6, y = 2.05, label = length), size = 8) +
  facet_wrap(~seqid)+
  theme(
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    # axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.x = element_line(),
    strip.text = element_text(size = 16, color = "black")) +
  # xlab("QX1410 Genome Coordinates (Mb)") +
  coord_cartesian(xlim = c(8.166162, 8.167273))
II_left

II_right <- ggplot() +
  geom_rect(data = briggsae %>% dplyr::filter(seqid == "II"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1, ymax = 2), fill = "#53886C") + 
  geom_text(data = briggsae %>% dplyr::filter(seqid == "II"), aes(x = ((start + end) / 2) / 1e6, y = 2.05, label = length), size = 8) +
  facet_wrap(~seqid)+
  theme(
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    # axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.x = element_line(),
    strip.text = element_text(size = 16, color = "black")) +
  # xlab("QX1410 Genome Coordinates (Mb)") +
  coord_cartesian(xlim = c(8.668291, 8.668791))
II_right

IV <- ggplot() +
  geom_rect(data = briggsae %>% dplyr::filter(seqid == "IV"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1, ymax = 2), fill = "#53886C") + 
  geom_text(data = briggsae %>% dplyr::filter(seqid == "IV"), aes(x = ((start + end) / 2) / 1e6, y = 2.05, label = length), size = 8) +
  facet_wrap(~seqid)+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.x = element_line(),
    strip.text = element_text(size = 16, color = "black")) +
  xlab("QX1410 genome coordinates (Mb)") +
  coord_cartesian(xlim = c(3.494880, 3.495320))
IV


X <- ggplot() +
  geom_rect(data = briggsae %>% dplyr::filter(seqid == "X"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1, ymax = 2), fill = "#53886C") + 
  geom_text(data = briggsae %>% dplyr::filter(seqid == "X"), aes(x = ((start + end) / 2) / 1e6, y = 2.02, label = length), size = 5) +
  facet_wrap(~seqid)+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.x = element_line(),
    strip.text = element_text(size = 16, color = "black")) +
  xlab("QX1410 genome coordinates (Mb)") +
  coord_cartesian(xlim = c(15.607495, 15.612189))
X


QX <- cowplot::plot_grid(
  II_left, II_right, IV, X,
  rel_heights = c(1, 1),
  ncol = 2,
  align = "v",
  axis = "lr")
QX


# Now in AF16:
af <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/briggsae/results/predictions/filtered_gff/AF16.PRE.gff") %>% 
  dplyr::mutate(length = end-start) %>%
  dplyr::filter(grepl("Mir-36", attributes))

II_left <- ggplot() +
  geom_rect(data = af %>% dplyr::filter(seqid == "II"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1, ymax = 2), fill = "seagreen3") + 
  geom_text(data = af %>% dplyr::filter(seqid == "II"), aes(x = ((start + end) / 2) / 1e6, y = 2.05, label = length), size = 8) +
  facet_wrap(~seqid)+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.x = element_line(),
    strip.text = element_text(size = 16, color = "black")) + 
  xlab("AF16 Genome Coordinates (Mb)") +
  coord_cartesian(xlim = c(8.377474, 8.378586))
II_left

II_right <- ggplot() +
  geom_rect(data = af %>% dplyr::filter(seqid == "II"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1, ymax = 2), fill = "seagreen3") + 
  geom_text(data = af %>% dplyr::filter(seqid == "II"), aes(x = ((start + end) / 2) / 1e6, y = 2.05, label = length), size = 8) +
  facet_wrap(~seqid)+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.x = element_line(),
    strip.text = element_text(size = 16, color = "black")) +
  xlab("AF16 Genome Coordinates (Mb)") +
  coord_cartesian(xlim = c(8.891273, 8.891535))
II_right


X <- ggplot() +
  geom_rect(data = af %>% dplyr::filter(seqid == "X"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1, ymax = 2), fill = "seagreen3") + 
  geom_text(data = af %>% dplyr::filter(seqid == "X"), aes(x = ((start + end) / 2) / 1e6, y = 2.02, label = length), size = 5) +
  facet_wrap(~seqid)+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = NA),
    axis.line.x = element_line(),
    strip.text = element_text(size = 16, color = "black")) +
  xlab("AF16 genome coordinates (Mb)") +
  coord_cartesian(xlim = c(15.370752, 15.375177))
X


QX <- cowplot::plot_grid(
  II_left, II_right, X,
  rel_heights = c(1),
  ncol = 3,
  align = "v",
  axis = "lr")
QX



# pan-miRNA visualization in elegans
all_predictions_ce <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/elegans/all_strain_predictions/ce_all_142_updated_20260212.tsv", col_names = c("contig","mirmachine","miRNA","start","end","score","strand",'phase',"gene_id","seed","strain")) %>% 
  dplyr::select("contig","start","end","strand","gene_id","seed","strain") #%>%
  dplyr::filter(!seed == "seed=(None)") ############### Removing miRNA family members that do not have seed support

most_mir35 <- all_predictions_ce %>%
  dplyr::filter(gene_id == "Mir-36") %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(count = n())# %>%
  dplyr::ungroup() %>%
  dplyr::filter(count == max(count))
  
ce_strains <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/elegans/all_strain_predictions/ce_142_strains.tsv", col_names = "strain") %>% dplyr::pull()

counts <- all_predictions_ce %>%
  dplyr::group_by(gene_id,strain) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup()

N2_counts <- counts %>% dplyr::filter(strain == "N2") %>% dplyr::select(strain,gene_id,count) %>% dplyr::distinct(gene_id, .keep_all = T) %>% dplyr::mutate(fewest = count, most = count) %>% dplyr::rename(average = count) %>% dplyr::mutate(strain_count = 1)
CGC1_counts <- counts %>% dplyr::filter(strain == 'CGC1') %>% dplyr::select(strain,gene_id,count) %>% dplyr::distinct(gene_id, .keep_all = T) %>% dplyr::mutate(fewest = count, most = count) %>% dplyr::rename(average = count)%>% dplyr::mutate(strain_count = 1)

strain_count <- all_predictions_ce %>%
  dplyr::filter(strain != "N2" & strain != "CGC1") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(n_strain = n_distinct(strain), .groups = "drop")

average_WSs <- counts %>% dplyr::filter(strain != "N2" & strain != "CGC1") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(average = mean(count)) %>%
  dplyr::mutate(fewest = min(count)) %>%
  dplyr::mutate(most = max(count)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(strain.x = "140_WSs") %>%
  dplyr::select(-strain) %>%
  dplyr::rename(strain = strain.x) %>%
  dplyr::select(strain,gene_id,average,fewest,most) %>%
  dplyr::distinct(gene_id, .keep_all = T) %>%
  dplyr::left_join(strain_count, by = "gene_id")

merged <- N2_counts %>% dplyr::bind_rows(CGC1_counts,average_WSs) %>% dplyr::mutate(strain = factor(strain, levels = c("N2","CGC1","140_WSs")))

mir193 <- merged %>% dplyr::filter(gene_id == "Mir-193") %>% dplyr::mutate(strain = "ECA347")

ggplot() + 
  geom_col(data = merged, aes(x = strain, y = average, fill = strain)) +
  geom_errorbar(data = merged %>% dplyr::filter(strain == "140_WSs"), aes(x = strain, ymin = fewest, ymax = most), width = 0.2) +
  geom_text(data = merged %>% dplyr::filter(strain == "140_WSs"), aes(x = strain, y = average + 1.75, label = n_strain), size = 3.5) +
  scale_fill_manual(values = c("N2" = "#DB6333", "CGC1" = "darkorange1", "140_WSs" = "magenta3")) +
  facet_wrap(~gene_id) +
  labs(y = "Count") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    axis.title.x = element_blank())

strain_count_ce <- all_predictions_ce %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(n_strain = n_distinct(strain), .groups = "drop")

ggplot() + 
  geom_col(data = strain_count_ce, aes(x = gene_id, y = n_strain), fill = '#DB6333', alpha = 0.7) +
  geom_hline(yintercept = 142, linetype = 'dashed', color =' black', linewidth = 1.5)+
  labs(y = "Count among all 142 Ce strains") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 12, color = 'black'),
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title.y = element_text(size = 14, color = 'black', face = 'bold')) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0,145))





# pan-miRNA visualization in briggsae
all_predictions <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/briggsae/all_strain_predictions/64_Cb_final_20260212.tsv", col_names = c("contig","mirmachine","miRNA","start","end","score","strand",'phase',"gene_id","seed","strain")) %>% 
  dplyr::select("contig","start","end","strand","gene_id","seed","strain") %>%
  dplyr::filter(!seed == "seed=(None)") ############### Removing miRNA family members that do not have seed support

cb_strains <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/briggsae/all_strain_predictions/cb_64_strains.tsv", col_names = "strain") %>% dplyr::pull()

counts <- all_predictions %>%
  dplyr::group_by(gene_id,strain) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup()

QX_counts <- counts %>% dplyr::filter(strain == "QX1410") %>% dplyr::select(strain,gene_id,count) %>% dplyr::distinct(gene_id, .keep_all = T) %>% dplyr::mutate(fewest = count, most = count) %>% dplyr::rename(average = count) %>% dplyr::mutate(strain_count = 1)
AF_counts <- counts %>% dplyr::filter(strain == 'AF16') %>% dplyr::select(strain,gene_id,count) %>% dplyr::distinct(gene_id, .keep_all = T) %>% dplyr::mutate(fewest = count, most = count) %>% dplyr::rename(average = count)%>% dplyr::mutate(strain_count = 1)
VX_counts <- counts %>% dplyr::filter(strain == 'VX34') %>% dplyr::select(strain,gene_id,count) %>% dplyr::distinct(gene_id, .keep_all = T) %>% dplyr::mutate(fewest = count, most = count) %>% dplyr::rename(average = count)%>% dplyr::mutate(strain_count = 1)
CGC2_counts <- counts %>% dplyr::filter(strain == 'CGC2') %>% dplyr::select(strain,gene_id,count) %>% dplyr::distinct(gene_id, .keep_all = T) %>% dplyr::mutate(fewest = count, most = count) %>% dplyr::rename(average = count)%>% dplyr::mutate(strain_count = 1)


strain_count_cb <- all_predictions %>%
  dplyr::filter(strain != "AF16" & strain != "QX1410" & strain != "VX34" & strain != "CGC2") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(n_strain = n_distinct(strain), .groups = "drop")

average_WSs <- counts %>% dplyr::filter(strain != "QX1410" & strain != "AF16" & strain != "VX34" & strain != "CGC2") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(average = mean(count)) %>%
  dplyr::mutate(fewest = min(count)) %>%
  dplyr::mutate(most = max(count)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(strain.x = "60_WSs") %>%
  dplyr::select(-strain) %>%
  dplyr::rename(strain = strain.x) %>%
  dplyr::select(strain,gene_id,average,fewest,most) %>%
  dplyr::distinct(gene_id, .keep_all = T) %>%
  dplyr::left_join(strain_count_cb, by = "gene_id")

merged <- QX_counts %>% dplyr::bind_rows(AF_counts, VX_counts, CGC2_counts, average_WSs) %>% dplyr::mutate(strain = factor(strain, levels = c("QX1410","AF16","CGC2", "VX34", "60_WSs")))

ggplot() + 
  geom_col(data = merged, aes(x = strain, y = average, fill = strain)) +
  geom_errorbar(data = merged %>% dplyr::filter(strain == "60_WSs"), aes(x = strain, ymin = fewest, ymax = most), width = 0.2) +
  geom_text(data = merged %>% dplyr::filter(strain == "60_WSs"), aes(x = strain, y = average + 8, label = n_strain), size = 3.5) +
  scale_fill_manual(values = c("QX1410" = "#53886C", "AF16" = "forestgreen", "CGC2" = "seagreen3", "VX34" = "purple", "60_WSs" = "violet")) +
  facet_wrap(~gene_id) +
  labs(y = "Count") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    axis.title.x = element_blank())


strain_count_cb_all <- all_predictions %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(n_strain = n_distinct(strain), .groups = "drop")

ggplot() + 
  geom_col(data = strain_count_cb_all, aes(x = gene_id, y = n_strain), fill = '#53886C', alpha = 0.7) +
  geom_hline(yintercept = 64, linetype = 'dashed', color =' black', linewidth = 1.5)+
  labs(y = "Count among all 64 Cb strains") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 12, color = 'black'),
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title.y = element_text(size = 14, color = 'black', face = 'bold')) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0,65))






# Plotting conservation
conserv <- strain_count_ce %>% dplyr::left_join(strain_count_cb_all, by = 'gene_id') %>%
  dplyr::rename(elegans = n_strain.x, briggsae = n_strain.y) %>%
  dplyr::mutate(briggsae = ifelse(is.na(briggsae),0,briggsae)) %>%
  dplyr::mutate(elegans = elegans / 142, briggsae = briggsae / 64) %>%
  dplyr::select(gene_id,elegans,briggsae)

df_long <- conserv %>%
  tidyr::pivot_longer(
    cols = c(elegans, briggsae),
    names_to = "Species",
    values_to = "prop"
  )

ggplot(df_long, aes(x = gene_id, y = prop, color = Species)) +
  geom_point(size = 3.2, position = position_dodge(width = 0.8)) +
  geom_vline(xintercept = seq_along(unique(df_long$gene_id))+0.5, 
             color = "grey80", size = 0.3) +
  scale_color_manual(values = c("elegans" = "#DB6333", 
                               "briggsae" = "#53886C")) +
  labs(x = "Gene ID", y = "Proportion of strains", fill = "Species") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        legend.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 14, color = 'black'),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_text(size = 14, color = 'black', face = 'bold')) 







### Looking to see how seed specificiation classifies miRNA family
seed_grouping_seed <- all_predictions %>%
  dplyr::filter(seed != "seed=(None)") %>%
  dplyr::group_by(seed) %>%
  dplyr::mutate(genes_per_seed = n_distinct(gene_id)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gene_id, seed, genes_per_seed) %>%
  dplyr::filter(genes_per_seed == 2)

seed_grouping_gene_id <- all_predictions %>%
  dplyr::filter(seed != "seed=(None)") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(unique_seeds = n_distinct(seed)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gene_id, seed, unique_seeds) %>% 
  dplyr::filter(unique_seeds > 1)
  
seed_grouping_seed_ce <- all_predictions_ce %>%
  dplyr::filter(seed != "seed=(None)") %>%
  dplyr::group_by(seed) %>%
  dplyr::mutate(genes_per_seed = n_distinct(gene_id)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gene_id, seed, genes_per_seed) %>%
  dplyr::filter(genes_per_seed == 2)

seed_grouping_gene_id_ce <- all_predictions_ce %>%
  dplyr::filter(seed != "seed=(None)") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(unique_seeds = n_distinct(seed)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gene_id, seed, unique_seeds) %>% 
  dplyr::filter(unique_seeds > 1)















