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



briggsae <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/briggsae/results/predictions/filtered_gff/Caendorhabditis_briggsae.PRE.gff") %>% 
  dplyr::mutate(length = end-start)

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









# Looking at alignments of strains that do not have a lys-6 miRNA
############ elegans
nucmer_ce <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-L1,-L2,-LENR,-LENQ) %>% dplyr::filter(strain != "ECA396")

all_predictions <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/elegans/all_strain_predictions/ce_all_142_updated_20260212.tsv", col_names = c("contig","mirmachine","miRNA","start","end","score","strand",'phase',"gene_id","seed","strain")) %>% 
  dplyr::select("contig","start","end","strand","gene_id","seed","strain") 

lsy6_strains <- all_predictions %>% dplyr::filter(gene_id == "Lsy-6") %>% dplyr::distinct(strain) %>% dplyr::pull()
ce_strains <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/elegans/all_strain_predictions/ce_142_strains.tsv", col_names = "strain") %>% dplyr::pull()

no_lsy6_strains <- setdiff(ce_strains,lsy6_strains)

N2_lsy6 <- all_predictions %>% dplyr::rename(strain_ = strain) %>% dplyr::filter(gene_id == "Lsy-6" & strain_ == "N2")

aln_lsy6_ce <- nucmer_ce %>% dplyr::filter(strain %in% no_lsy6_strains & N2_chr == "V")



one <- ggplot(data = aln_lsy6_ce %>% dplyr::filter(strain == "ECA1843")) + 
  geom_rect(data = N2_lsy6, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.5) +
  geom_segment(aes(x = N2S /1e6, xend = N2E /1e6, y = WSS /1e6, yend = WSE /1e6, color = IDY), size = 1.5) + 
  scale_color_gradient(low = "gold", high = "black", limits = c(80, 100), name = "Percent ID") +
  facet_wrap(~strain, nrow = 1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 14, color = 'black'),
    strip.text = element_text(size = 16, color = 'black'),
    legend.position = 'none'
  ) +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)") +
  coord_cartesian(xlim = c(10.647137, 10.647401), ylim = c(11.52475, 11.5251))
one

two <- ggplot(data = aln_lsy6_ce %>% dplyr::filter(strain == "ECA2367")) + 
  geom_rect(data = N2_lsy6, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.5) +
  geom_segment(aes(x = N2S /1e6, xend = N2E /1e6, y = WSS /1e6, yend = WSE /1e6, color = IDY), size = 1.5) + 
  scale_color_gradient(low = "gold", high = "black", limits = c(80, 100), name = "Percent ID") +
  facet_wrap(~strain, nrow = 1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 14, color = 'black'),
    strip.text = element_text(size = 16, color = 'black'),
    legend.position = 'none'
  ) +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)") +
  coord_cartesian(xlim = c(10.647137, 10.647401), ylim = c(11.2386, 11.2389))
two

three <- ggplot(data = aln_lsy6_ce %>% dplyr::filter(strain == "TWN2530")) + 
  geom_rect(data = N2_lsy6, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.5) +
  geom_segment(aes(x = N2S /1e6, xend = N2E /1e6, y = WSS /1e6, yend = WSE /1e6, color = IDY), size = 1.5) + 
  scale_color_gradient(low = "gold", high = "black", limits = c(80, 100), name = "Percent ID") +
  facet_wrap(~strain, nrow = 1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 14, color = 'black'),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 16, color = 'black'),
    legend.text = element_text(size = 12, color = 'black'),
    legend.title = element_text(size = 14, color = 'black')
  ) +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)") +
  coord_cartesian(xlim = c(10.647137, 10.647401), ylim = c(0.81078, 0.8111))
three

cowplot::plot_grid(
  one, two, three,
  nrow = 1
)





############ briggsae
nucmer_cb <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/briggsae/nucmer_aln_WSs/all_61_aln.clean.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","QX_chr","contig","strain")) %>%
  dplyr::select(-L1,-L2,-LENR,-LENQ)


all_predictions_cb <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/briggsae/all_strain_predictions/64_Cb_final_20260212.tsv", col_names = c("contig","mirmachine","miRNA","start","end","score","strand",'phase',"gene_id","seed","strain")) %>% 
  dplyr::select("contig","start","end","strand","gene_id","seed","strain") 

lsy6_strains <- all_predictions_cb %>% dplyr::filter(gene_id == "Lsy-6") %>% dplyr::distinct(strain) %>% dplyr::pull()
cb_strains <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/briggsae/all_strain_predictions/cb_64_strains.tsv", col_names = "strain") %>% dplyr::pull()

no_lsy6_strains <- setdiff(cb_strains,lsy6_strains)

QX_lsy6 <- all_predictions_cb %>% dplyr::rename(strain_ = strain) %>% dplyr::filter(gene_id == "Lsy-6" & strain_ == "QX1410")

aln_lys6_cb <- nucmer_cb %>% dplyr::filter(strain %in% no_lsy6_strains & QX_chr == "V")


cb_ten <- ggplot(data = aln_lys6_cb) + 
  geom_rect(data = QX_lsy6, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = "#53886C", alpha = 0.5) +
  geom_segment(aes(x = N2S /1e6, xend = N2E /1e6, y = WSS /1e6, yend = WSE /1e6, color = IDY), size = 1.5) +
  scale_color_gradient(low = "gold", high = "black", limits = c(80, 100), name = "Percent ID") +
  facet_wrap(~strain, nrow = 2, ncol = 5) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 14, color = 'black'),
    strip.text = element_text(size = 16, color = 'black'),
    legend.text = element_text(size = 12, color = 'black'),
    legend.title = element_text(size = 14, color = 'black')
  ) +
  labs(x = "QX1410 genome position (Mb)", y = "Wild strain contig position (Mb)") +
  coord_cartesian(xlim = c(10.585, 10.585117))
cb_ten









####### Does CGC2 have the same mir35 expansion on chrom V as QX1410??
CGC_QX <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/alignment/CGC2_final.transformed.tsv", col_names = c("QXS","QXE","CGS","CGE","L1","L2","IDY","LENR","LENQ","QX_chrom","CGC2_scaffold")) 

all_predictions <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/briggsae/all_strain_predictions/64_Cb_final_20260212.tsv", col_names = c("contig","mirmachine","miRNA","start","end","score","strand",'phase',"gene_id","seed","strain")) %>% 
  dplyr::select("contig","start","end","strand","gene_id","seed","strain") %>%
  dplyr::filter((strain == "CGC2" | strain == "QX1410") & gene_id == "Mir-36") %>%
  dplyr::filter(contig == "X") %>% 
  dplyr::group_by(strain) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup()

chromV_mir35 <- ggplot(CGC_QX %>% dplyr::filter(QX_chrom == "X" & CGC2_scaffold == "X")) +
  geom_segment(aes(x = QXS / 1e6, xend = QXE / 1e6, y = CGS / 1e6, yend = CGE / 1e6), color = 'black', linewidth = 2) +
  geom_rect(data = all_predictions %>% dplyr::filter(strain == "QX1410"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = "forestgreen", alpha = 0.3) +
  geom_rect(data = all_predictions %>% dplyr::filter(strain == "CGC2"), aes(xmin = -Inf, xmax = Inf, ymin = start / 1e6, ymax = end / 1e6), fill = "magenta", alpha = 0.4) +
  facet_wrap(~QX_chrom, scales = "free") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 20, color = 'black', face = 'bold'),
    legend.text = element_text(size = 12, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    panel.border = element_rect(fill = NA)) +
  labs(y = "CGC2 genome coordinates (Mb)", x = "QX1410 genome coordinates (Mb)") +
  coord_cartesian(xlim = c(15.607, 15.6125), ylim = c(15.789, 15.794))
chromV_mir35


all_predictions <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/briggsae/all_strain_predictions/64_Cb_final_20260212.tsv", col_names = c("contig","mirmachine","miRNA","start","end","score","strand",'phase',"gene_id","seed","strain")) %>% 
  dplyr::select("contig","start","end","strand","gene_id","seed","strain") %>%
  dplyr::filter((strain == "CGC2" | strain == "QX1410") & gene_id == "Mir-36") %>%
  dplyr::filter(contig == "II") %>% 
  dplyr::group_by(strain) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup()

chromII_mir35 <- ggplot(CGC_QX %>% dplyr::filter(QX_chrom == "II" & CGC2_scaffold == "II")) +
  geom_segment(aes(x = QXS / 1e6, xend = QXE / 1e6, y = CGS / 1e6, yend = CGE / 1e6), color = 'black', linewidth = 2) +
  geom_rect(data = all_predictions %>% dplyr::filter(strain == "QX1410"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = "forestgreen", alpha = 0.3) +
  geom_rect(data = all_predictions %>% dplyr::filter(strain == "CGC2"), aes(xmin = -Inf, xmax = Inf, ymin = start / 1e6, ymax = end / 1e6), fill = "magenta", alpha = 0.4) +
  facet_wrap(~QX_chrom, scales = "free") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 20, color = 'black', face = 'bold'),
    legend.text = element_text(size = 12, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    panel.border = element_rect(fill = NA)) +
  labs(y = "CGC2 genome coordinates (Mb)", x = "QX1410 genome coordinates (Mb)") +
  coord_cartesian(xlim = c(8.1660,8.1675), ylim = c(8.203,8.2045))
chromII_mir35






############## briggsae mir-10
read_clustalo_block_pim <- function(fn) {
  lines <- readLines(fn, warn = FALSE)
  
  rows <- list()
  ids  <- character()
  cur_id <- NULL
  
  for (ln in lines) {
    # New row starts with "number:"
    if (grepl("^\\s*\\d+:", ln)) {
      # Remove leading index
      ln_clean <- sub("^\\s*\\d+:\\s*", "", ln)
      
      # Split into tokens
      parts <- strsplit(ln_clean, "\\s+")[[1]]
      
      cur_id <- parts[1]
      ids <- c(ids, cur_id)
      
      nums <- as.numeric(parts[-1])
      rows[[cur_id]] <- nums
    } else {
      # Continuation line
      nums <- as.numeric(strsplit(trimws(ln), "\\s+")[[1]])
      rows[[cur_id]] <- c(rows[[cur_id]], nums)
    }
  }
  
  n <- length(ids)
  mat <- matrix(NA_real_, n, n)
  rownames(mat) <- ids
  colnames(mat) <- ids
  
  for (i in seq_len(n)) {
    mat[i, seq_along(rows[[ids[i]]])] <- rows[[ids[i]]]
  }
  
  mat
}

fn <- "/vast/eande106/projects/Lance/THESIS_WORK/misc/mir35/processed_data/clustalo-I20260212-205150-0917-74873756-p1m.pim"

pid <- read_clustalo_block_pim(fn)

dim(pid)


pid_long <- as.data.frame(pid) |>
  tibble::rownames_to_column("seq1") |>
  tidyr::pivot_longer(-seq1, names_to = "seq2", values_to = "pident")

ggplot(pid_long, aes(seq1, seq2, fill = pident)) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient(low = "gold", high = "darkgreen") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "%ID")

































