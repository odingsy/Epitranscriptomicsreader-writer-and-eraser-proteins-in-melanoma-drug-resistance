# Library, global functions and parameters  ----
library(tidyverse)
library(ggrepel)
library(scales)
base <- '/Users/shiyuanguo/Library/CloudStorage/GoogleDrive-sguo039@ucr.edu/My Drive/PhD_study/IGR_RWE/1.IGR37xp_manuscript'
plot_shape = 20; plot_size = 0.5; plot_alpha = 0.4; plot_pos = '#F8766D'; global_linewidth <- 0.3 # ggploting parameters


# loading SILAC_peptidesRatio between IGR37xp vs IGR37 ----
# obtain the injection info from both replicates
injectionfiles = paste0(file.path(base, 'data', 'injectionLists', '12min_max40txsits'), '_000', 1:5, '.csv')
injectionList <- lapply(injectionfiles, function(n){
  readr::read_csv(n) %>%
    dplyr::mutate(Comment = stringr::str_remove_all(Comment,' \\((light|heavy)\\)') %>% 
                    stringr::str_remove_all('\\[[:graph:][:digit:]{2}\\.[:digit:]{6}\\]')) %>%
    .[['Comment']] %>% 
    unique(.)
})

# skyline document grid > SILAC_peptideRatio > export as 'SILAC_peptideRatio.csv'
ptbl <- read_csv(file.path(base, 'data', 'SILAC_peptideRatio12.csv')) %>% 
  filter(Protein != 'sp|P02769|ALBU_BOVIN-standard') %>%
  ############ filter replicate base on injection number.
  mutate(prmlist = case_when(
    Peptide %in% injectionList[[1]] ~ 1,
    Peptide %in% injectionList[[2]] ~ 2,
    Peptide %in% injectionList[[3]] ~ 3,
    Peptide %in% injectionList[[4]] ~ 4,
    Peptide %in% injectionList[[5]] ~ 5,
    TRUE ~ 0)) %>%
  group_by(prmlist) %>%
  mutate(repli = str_extract(Replicate, '\\d$') %>% as.numeric(),
         fr = str_extract(Replicate, '[:alpha:]')) %>% 
  filter(repli == prmlist) %>% 
  ungroup() %>%
  ############ end 
  mutate(`Total Area` = as.numeric(`Total Area`))%>% 
  dplyr::select(Peptide, `Protein Name`, fr, `Total Area`, `Isotope Label Type`) %>% 
  pivot_wider(names_from = c(fr, `Isotope Label Type`), values_from = `Total Area`)%>%
  mutate(Fwd =   F_heavy / F_light, 
         Rvs =   R_light / R_heavy) %>% 
  mutate(`Protein Name` = case_when(
    str_detect(`Protein Name`, 'ENSP00000409983') ~ paste('IGF2BP3', `Protein Name`, sep = ' '),
    str_detect(`Protein Name`, 'ENSP00000384369') ~ paste('METTL15', `Protein Name`, sep = ' '),
    TRUE ~ `Protein Name`
  ))%>% 
  separate_wider_regex('Protein Name', c(proteinName = ".*?", " |/", ".*")) %>% 
  filter(proteinName != 'METTL9') # METTL9 is not a RWE protein but included in the RWE library because it is in METTL family. 

# 3rd and 4th replicates 
injectionfiles = paste0(file.path(base, 'data', 'injectionLists', 'maxtransition52_9minWindow_3inj'), '_000', 1:3, '.csv')
injectionList <- lapply(injectionfiles, function(n){
  readr::read_csv(n) %>%
    dplyr::mutate(Comment = stringr::str_remove_all(Comment,' \\((light|heavy)\\)') %>% 
                    stringr::str_remove_all('\\[[:graph:][:digit:]{2}\\.[:digit:]{6}\\]')) %>%
    .[['Comment']] %>% 
    unique(.)
})

# load SILAC_peptideRatio34 for 3rd and 4th replicates 
# IV199_15: IGR37xp(L) + IGR37(H) Rvs
# IV199_16: IGR37(L) + IGR37xp(H) Fwd
ptbl34 <- read_csv(file.path(base, 'data', 'SILAC_peptideRatio34.csv')) %>% 
  filter(Protein != 'sp|P02769|ALBU_BOVIN-standard') %>%
  ############ filter replicate base on injection number.
  mutate(prmlist = case_when(
    Peptide %in% injectionList[[1]] ~ 1,
    Peptide %in% injectionList[[2]] ~ 2,
    Peptide %in% injectionList[[3]] ~ 3,
    TRUE ~ 0)) %>%
  group_by(prmlist) %>%
  mutate(repli = str_extract(Replicate, '\\d$') %>% as.numeric(),
         fr = str_extract(Replicate, '1(5|6)'),
         fr = case_when(
           fr == 16 ~ 'F',
           fr == 15 ~ 'R',
         )) %>% 
  filter(repli == prmlist) %>% 
  ungroup() %>%
  ############ end 
  mutate(`Total Area` = as.numeric(`Total Area`))%>% 
  dplyr::select(Peptide, `Protein Name`, fr, `Total Area`, `Isotope Label Type`) %>% 
  pivot_wider(names_from = c(fr, `Isotope Label Type`), values_from = `Total Area`) %>%
  mutate(Fwd =   F_heavy / F_light, 
         Rvs =   R_light / R_heavy) %>% 
  mutate(`Protein Name` = case_when(
    str_detect(`Protein Name`, 'ENSP00000409983') ~ paste('IGF2BP3', `Protein Name`, sep = ' '),
    str_detect(`Protein Name`, 'ENSP00000384369') ~ paste('METTL15', `Protein Name`, sep = ' '),
    str_detect(`Protein Name`, 'ENSP00000358563') ~ paste('DKC1', `Protein Name`, sep = ' '),
    TRUE ~ `Protein Name`
  )) %>% 
  separate_wider_regex('Protein Name', c(proteinName = ".*?", " |/", ".*")) %>% 
  filter(proteinName != 'METTL9') # METTL9 is not a RWE protein but included in the RWE library because it is in METTL family. 


# protein level by averaging all the constituting peptides
prottbl <- bind_rows(
  ptbl %>% rename(F1 = Fwd, R1 = Rvs) %>% 
    select(Peptide, proteinName, F1, R1) %>% 
    pivot_longer(cols = F1:R1, names_to = 'replicates', values_to = 'ratio'),
  ptbl34 %>% rename(F2 = Fwd, R2 = Rvs) %>% 
    select(Peptide, proteinName, F2, R2) %>% 
    pivot_longer(cols = F2:R2, names_to = 'replicates', values_to = 'ratio')) %>% 
  pivot_wider(names_from = 'replicates', values_from = 'ratio') %>% 
  group_by(proteinName) %>% 
  summarise(across(matches('F|R'), ~mean(.x, na.rm = TRUE))) %>% 
  rowwise() %>%
  mutate(log2Fwd = mean(c(F1, F2), na.rm = TRUE) %>% log2(),
         log2Rvs = mean(c(R1, R2), na.rm = TRUE) %>% log2(),
         fc = mean(c(F1, F2, R1, R2), na.rm = TRUE),
         rsd = sd(c(F1, F2, R1, R2), na.rm = TRUE)/fc,
         log2fc = log2(mean(c(F1, F2, R1, R2), na.rm = TRUE))) %>% 
  ungroup() %>% 
  arrange(desc(log2fc))


# output quantifiable protein arranged in desc order protein level fc
pplist <- prottbl %>% 
  dplyr::mutate(rnum = row_number()) %>%
  mutate(shade = case_when(
    rnum %% 2 == 1 ~ 0,
    rnum %% 2 == 0 ~ 1,
  )) %>%
  dplyr::select(-log2Fwd, -log2Rvs)
wb <- openxlsx::createWorkbook();
openxlsx::addWorksheet(wb, 'test');
lsty1 <-  openxlsx::createStyle(bgFill = "#d6d2d2")
openxlsx::conditionalFormatting(wb, 'test', rows = 1:(nrow(pplist)+1), cols = 1:ncol(pplist), rule = paste0("$", LETTERS[which(colnames(pplist) == "shade")], "1==", 1), style = lsty1)
openxlsx::writeData(wb, 'test', pplist);
openxlsx::openXL(wb)

# protein dot plot 
cf <- 1
gt_ls <- prottbl %>% filter(log2fc > cf) %>% na.omit()
lt_ls <- prottbl %>% filter(log2fc < -cf) %>% na.omit()
p <- prottbl %>% 
  ggplot(aes(x = log2Fwd, y = log2Rvs), shape = plot_shape)+
  geom_point(color = 'grey', size = plot_size, alpha = plot_alpha)+
  geom_point(data = bind_rows(gt_ls, lt_ls), color = plot_pos, size = plot_size)+
  ggrepel::geom_text_repel(data = bind_rows(gt_ls, tail(lt_ls,5)), aes(label = proteinName), size= 3, force_pull = 10, force = 1, min.segment.length = 0.1, max.overlaps = Inf)+
  scale_x_continuous(breaks = -3:3)+
  scale_y_continuous(breaks = -3:3)+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(-3, 3), ylim = c(-3, 3), xintercept = 0, yintercept = 0, ratio = 1)+
  labs(x = "Log2(IGR37xp/IGR37), Forward", y = "Log2(IGR37xp/IGR37), Reverse")+
  theme_classic()+ # increase x limit
  theme(legend.position = "none",
        # panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank())
ggsave(filename = file.path(base, 'figures/protein_dotplot.pdf'), device = 'pdf', plot = p, width = 5, height = 5, units = 'in')

# barplot
tb_num <- 10
p <- prottbl %>% 
  na.omit() %>% 
  mutate(proteinName = fct_reorder(proteinName, desc(log2fc))) %>%
  slice(1:tb_num, (n()-tb_num + 1):n()) %>% 
  select(proteinName, F1:R2) %>% 
  pivot_longer(-proteinName, names_to = 'replicates', values_to = 'ratios') %>% 
  mutate(fr = str_extract(replicates, '[:alpha:]'),
         tb = ratios > 1,
         log2ratios = log2(ratios)) %>% 
  ggplot(aes(proteinName, log2ratios))+
  stat_summary(aes(fill = tb), geom = 'bar', fun = 'mean' , na.rm = TRUE, width = .5)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), 
               color = 'black', width = 0.2,linewidth = 0.16)+
  geom_point(aes(color = fr), size = 0.7)+
  ylab('Log2(IGR37xp / IGR37)')+
  ggh4x::coord_axes_inside(labels_inside = TRUE, ylim = c(-3, 3), yintercept = 0)+
  scale_color_manual(values = hcl.colors(hcl.pals('diverging')[2], n = 8)[c(1,8)])+ # dot color
  # scale_fill_manual(values = )+
  guides(fill = 'none', color = guide_legend(override.aes = list(size = 0.5)))+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.2, hjust = 0.2),
        axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.8))
ggsave(filename = file.path(base, 'figures/protein_batplot.pdf'), device = 'pdf', plot = p, width = 6.2, height = 2, units = 'in')


# OXPHOS protein quantification -----
tbl <- readxl::read_excel(file.path(base, 'data', 'ETC_annotation.xlsx'), sheet = 'quantification', range = readxl::cell_rows(1:37)) %>% 
  mutate(LH = ifelse(str_detect(Compound, '\\((K8|R6)\\)$'), 'H', 'L')) 


# shTRMU / shCtrl in IGR37xp cells 
# runs included in the analysis: IV200_6, IV200_9, IV200_13
kd <- tbl %>% 
  select(proteinName, LH, IV200_6, IV200_9, IV200_13) %>%
  pivot_longer(where(is.double), names_to = 'sample', values_to = 'amt') %>% 
  pivot_wider(names_from = 'LH', values_from = 'amt') %>%
  mutate(ratio = case_when(
    str_detect(sample, 'IV200_6') ~ H / L,
    str_detect(sample, 'IV200_(9|13)')  ~ L / H,
    .default = 0.0
  )) %>%
  filter(proteinName != 'ACTB') %>% 
  mutate(proteinName = ifelse(proteinName == 'ACTB_VAP', 'ACTB', proteinName)) %>% # replace back to ACTB to obtain ratio mean later
  filter(!is.na(ratio)) %>% 
  mutate(sample = str_extract(sample, '^([^_-]+_[^_-]+)?')) %>% # from gemini!! 
  group_by(proteinName, sample) %>% 
  summarize(mratio = mean(ratio, na.rm = TRUE),
            rsd = sd(ratio, na.rm = TRUE) / mratio) %>% 
  ungroup() %>% 
  mutate(whereEncode = ifelse(str_detect(proteinName, 'MT-'), 'MT', 'nuc')) %>% 
  dplyr::mutate(proteinName = fct_relevel(proteinName, 
                                          'TRMU', 'MTO1', 'ACTB', 'GAPDH',  #4
                                          'ETFA', 'SDHA', 'MRPL11' ,'TOMM40', #4
                                          'MT-ND6', 'MT-ND5', 'MT-ND4', 'MT-ND3', #4
                                          'MT-CYB', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'MT-CO2',#5
  ))


# MTO1 KO / IGR37xp 
# runs included in the analysis: IV198_8	IV198_12, IV198_7_20250505001525)
mto1ko <- tbl %>% 
  select(proteinName, LH, IV198_8, IV198_12) %>%  
  pivot_longer(where(is.double), names_to = 'sample', values_to = 'amt') %>% 
  pivot_wider(names_from = 'LH', values_from = 'amt') %>%
  mutate(ratio = case_when(
    str_detect(sample, 'IV198_(8|12)')  ~ H / L,
    str_detect(sample, 'IV198_(7|11)') ~ L / H,
    .default = 0.0
  )) %>% 
  # filter(proteinName != 'ACTB') %>% 
  mutate(proteinName = ifelse(proteinName == 'ACTB_VAP', 'ACTB', proteinName)) %>% # replace back to ACTB to obtain ratio mean later
  filter(!is.na(ratio)) %>% 
  mutate(sample = str_extract(sample, '^([^_-]+_[^_-]+)?')) %>% # from gemini!! 
  group_by(proteinName, sample) %>% 
  summarize(mratio = mean(ratio, na.rm = TRUE),
            rsd = sd(ratio, na.rm = TRUE) / mratio) %>% 
  ungroup() %>% 
  mutate(whereEncode = ifelse(str_detect(proteinName, 'MT-'), 'MT', 'nuc')) %>% 
  dplyr::mutate(proteinName = fct_relevel(proteinName, 
                                          'TRMU', 'MTO1', 'ACTB', 'GAPDH',  #4
                                          'ETFA', 'SDHA', 'MRPL11' ,'TOMM40', #4
                                          'MT-ND6', 'MT-ND5', 'MT-ND4', 'MT-ND3', #4
                                          'MT-CYB', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'MT-CO2',#5
  ))



# DONE IGR37xp / IGR37 
# runs included in the analysis: IV198_5_20250502113912, IV198_6_20250502145929, IV198_9
endo <- tbl %>% 
  select(proteinName, LH, IV198_5_20250502113912, IV198_6_20250502145929, IV198_9) %>% 
  pivot_longer(where(is.double), names_to = 'sample', values_to = 'amt') %>% 
  pivot_wider(names_from = 'LH', values_from = 'amt') %>% 
  mutate(ratio = case_when(
    str_detect(sample, 'IV198_(5|9)')  ~ H / L,
    str_detect(sample, 'IV198_(6|10)') ~ L / H, 
    .default = 0.0
  )) %>% 
  filter(proteinName != 'ACTB') %>%
  mutate(proteinName = ifelse(proteinName == 'ACTB_VAP', 'ACTB', proteinName)) %>% # replace back to ACTB to obtain ratio mean later
  filter(!is.na(ratio)) %>% 
  mutate(sample = str_extract(sample, '^([^_-]+_[^_-]+)?')) %>% # regex from gemini!! 
  group_by(proteinName, sample) %>% 
  summarize(mratio = mean(ratio, na.rm = TRUE),
            rsd = sd(ratio, na.rm = TRUE) / mratio) %>% 
  ungroup() %>% 
  mutate(whereEncode = ifelse(str_detect(proteinName, 'MT-'), 'MT', 'nuc'))%>% 
  dplyr::mutate(proteinName = fct_relevel(proteinName, 
                                          'TRMU', 'MTO1', 'ACTB', 'GAPDH',  #4
                                          'ETFA', 'SDHA', 'MRPL11' ,'TOMM40', #4
                                          'MT-ND6', 'MT-ND5', 'MT-ND4', 'MT-ND3', #4
                                          'MT-CYB', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'MT-CO2',#5
  ))


# 1.barplot endogenous
p1 <- endo %>% 
  ggplot(aes(proteinName, mratio))+
  stat_summary(aes(fill = whereEncode), geom = 'bar', fun = 'mean' , na.rm = TRUE, width = .5)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), 
               color = 'black', width = 0.2, linewidth = global_linewidth)+
  geom_hline(yintercept=1, linetype="dashed",color = "red", linewidth = global_linewidth)+
  geom_point(aes(color = sample), size = plot_size)+
  coord_flip()+
  ylab('IGR37xp / IGR37')+
  scale_color_manual(values = hcl.colors(hcl.pals('diverging')[2], n = 8)[c(1,2,7,8)])+ # dot color
  # scale_fill_manual(values = )+
  guides(fill = FALSE, color = FALSE)+
  theme_classic()+
  theme(axis.title.y = element_blank(), 
        legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.4))

# 2.barplot shTRMU
p2 <- kd %>%
  ggplot(aes(proteinName, mratio))+
  stat_summary(aes(fill = whereEncode), geom = 'bar', fun = 'mean' , na.rm = TRUE, width = .5)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), 
               color = 'black', width = 0.2, linewidth = global_linewidth)+
  geom_hline(yintercept=1, linetype="dashed",color = "red", linewidth = .5)+
  geom_point(aes(color = sample), size = plot_size)+
  coord_flip()+
  ylab('shTRMU / shCtrl in IGR37xp')+
  scale_color_manual(values = hcl.colors(hcl.pals('diverging')[2], n = 10)[c(1,3,6,10)])+ # dot color
  # scale_fill_manual(values = c(replication, repair, splicing))+
  guides(fill = FALSE, color = FALSE)+
  theme_classic()+
  theme(axis.title.y = element_blank(), 
        legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.4))

# 3. barplot MTO1ko
p3 <- mto1ko %>% 
  ggplot(aes(proteinName, mratio))+
  stat_summary(aes(fill = whereEncode), geom = 'bar', fun = 'mean' , na.rm = TRUE, width = .5)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), 
               color = 'black', width = 0.2, linewidth = global_linewidth)+
  geom_hline(yintercept=1, linetype="dashed",color = "red", linewidth = .5)+
  geom_point(aes(color = sample), size = plot_size)+
  coord_flip()+
  ylab('MTO1 KO / IGR37xp')+
  scale_color_manual(values = hcl.colors(hcl.pals('diverging')[2], n = 8)[c(1,5,8)])+ # dot color
  # scale_fill_manual(values = c(replication, repair, splicing))+
  guides(fill = FALSE, color = FALSE)+
  theme_classic()+
  theme(axis.title.y = element_blank(), 
        legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.4))



# 4.heatmap
p4 <- bind_rows(kd %>% mutate(trt = 'kd'),
          mto1ko %>% mutate(trt = 'mto1ko'),
          endo %>% mutate(trt = 'endo')) %>% 
  mutate(sample = factor(sample, levels = c('IV198_5', 'IV198_6', 'IV198_9','IV200_6', 'IV200_9', 'IV200_13', 'IV198_8', 'IV198_12'))) %>% 
  ggplot(aes(x = sample, y = proteinName, fill = log2(mratio))) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  scale_x_discrete(position = "top")+
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none") + # Rotate x-axis labels
  labs(title = "",
         x = "",
         y = "",
         fill = "log2Ratios")

p4n <- p4 + theme(legend.position="top")

ggsave(file.path(base, 'figures/1.pdf'), device = 'pdf', plot = p1, width = 3, height = 4, units = 'in')
ggsave(file.path(base, 'figures/2.pdf'), device = 'pdf', plot = p2, width = 3, height = 4, units = 'in')
ggsave(file.path(base, 'figures/3.pdf'), device = 'pdf', plot = p3, width = 3, height = 4, units = 'in')
ggsave(file.path(base, 'figures/4.pdf'), device = 'pdf', plot = p4, width = 4, height = 5, units = 'in')
ggsave(file.path(base, 'figures/5.pdf'), device = 'pdf', plot = p4n, width = 4, height = 4, units = 'in')

# visualizing the coverage of OXPHOS peptide using IRanges -----
library(IRanges)

## IRange Pp2. change xlim for better comparison 
plotRanges_generic <- function(xlim){
  function(x, main=deparse(substitute(x)), col="black", sep=0.5, ...){
    height <- 1
    if (is(xlim, "IntegerRanges") || is(xlim, "GenomicRanges"))
      xlim <- c(min(start(xlim)), max(end(xlim)))
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    plot.new()
    plot.window(xlim, c(0, max(bins)*(height + sep)))
    ybottom <- bins * (sep + height) - height
    rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
    title(main)
    axis(1)
  }
}

plotRanges <- plotRanges_generic(xlim = c(40, 145))

tbl <- readxl::read_excel('/Users/shiyuanguo/Library/CloudStorage/GoogleDrive-sguo039@ucr.edu/My Drive/PhD_study/IGR_RWE/1.IGR37xp_manuscript/data/ETC_annotation.xlsx', sheet = 'quantification', range = readxl::cell_rows(1:37)) %>% 
  mutate(`t start (min)` = ifelse(`t start (min)` < 40, 40, `t start (min)`),
         `t stop (min)` = ifelse(`t stop (min)`> 145, 145, `t stop (min)`)) # only plot from 40 to 145

par(mfrow = c(1,1))
ir <- IRanges(tbl$`t start (min)`, tbl$`t stop (min)`)
plotRanges(ir)

cov <- coverage(ir)
cov <- as.vector(cov)
mat <- cbind(seq_along(cov)-0.5, cov) 
d <- diff(cov) != 0
mat <- rbind(cbind(mat[d,1]+1, mat[d,2]), mat)
mat <- mat[order(mat[,1]), ]

plotRanges_generic(ir)(ir)
lines(mat, col = "red", lwd = 4)
axis(2)

# color picking for spectra
cp <- c(hcl.colors(hcl.pals('diverging')[1], n = 4),
           hcl.colors(hcl.pals('diverging')[2], n = 4),
           hcl.colors(hcl.pals('diverging')[3], n = 4),
           hcl.colors(hcl.pals('diverging')[4], n = 4),
           hcl.colors(hcl.pals('diverging')[5], n = 4),
           hcl.colors(hcl.pals('diverging')[6], n = 4),
           hcl.colors(hcl.pals('diverging')[6], n = 4),
           hcl.colors(hcl.pals('diverging')[7], n = 4),
           hcl.colors(hcl.pals('diverging')[8], n = 4),
           hcl.colors(hcl.pals('diverging')[9], n = 4),
           hcl.colors(hcl.pals('diverging')[10], n = 4),
           hcl.colors(hcl.pals('diverging')[11], n = 4),
           hcl.colors(hcl.pals('diverging')[12], n = 4),
           hcl.colors(hcl.pals('diverging')[13], n = 4),
           hcl.colors(hcl.pals('diverging')[14], n = 4),
           hcl.colors(hcl.pals('diverging')[15], n = 4),
           hcl.colors(hcl.pals('diverging')[16], n = 4)
           )[c(1,2,3,4,8,15,16,21,30,31,35,36,41,44,49,52, 64, 65)]
show_col(cp); cp
