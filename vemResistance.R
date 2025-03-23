# Library, global functions and parameters  ----
library(tidyverse)
library(magrittr)
library(ggrepel)
library(drc)
library(scales)
base <- '/Users/shiyuanguo/Library/CloudStorage/GoogleDrive-sguo039@ucr.edu/My Drive/PhD_study/IGR_RWE/1.IGR37xp_manuscript'
plot_shape = 20; plot_size = 0.5; plot_alpha = 0.4; plot_pos = '#F8766D'# ggploting parameters


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
ptbl <- read_csv(file.path(base, 'data', 'SILAC_peptideRatio.csv')) %>% 
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

# protein level by averaging all the constituting peptides
prottbl <- ptbl %>%
  group_by(proteinName) %>% # using `group_by` followed by `summarise` instead of `nest` to get two list-col for fwd and rvs respectively. 
  summarise(across(matches('Fwd|Rvs'), ~list(.x))) %>% 
  ungroup() %>%
  mutate(across(matches('Fwd|Rvs'), ~map_dbl(., ~mean(.x, na.rm = TRUE)))) %>% 
  rowwise() %>% 
  mutate(log2fc = log2(mean(c(Fwd, Rvs), na.rm = TRUE))) %>% 
  arrange(desc(log2fc))

# shortlist peptides based on mean protein ratios
th <- 2
gt_ls <- prottbl %>% filter(Fwd > th & Rvs > th) %>% pull(proteinName)
gt_tbl <- ptbl %>% 
  group_by(proteinName) %>% 
  filter(proteinName %in% gt_ls) %>% 
  ungroup()
lt_ls <- prottbl %>% filter(Fwd < 1/th & Rvs < 1/th) %>% pull(proteinName)
lt_tbl <- ptbl %>% 
  group_by(proteinName) %>% 
  filter(proteinName %in% lt_ls) %>% 
  ungroup()

# output tables
wb <- openxlsx::createWorkbook(); 
tn <- 'PRM_RWE_fulllist'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, ptbl);
tn <- 'gt'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, gt_tbl); 
tn <- 'lt'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, lt_tbl); 
openxlsx::openXL(wb)

# protein ratio plot 
p <- prottbl %>% 
  ggplot(aes(x = Fwd, y = Rvs), shape = plot_shape)+
  geom_point(color = 'grey', size = plot_size, alpha = plot_alpha)+
  geom_point(data = . %>% filter(proteinName %in% c(gt_ls, lt_ls)), color = plot_pos, size = plot_size)+
  ggrepel::geom_text_repel(data = . %>% filter(proteinName %in% c(gt_ls, tail(lt_ls, 6))), aes(label = proteinName), size= 3, force_pull = 10, force = 1, min.segment.length = 0.1, max.overlaps = Inf)+
  scale_x_continuous(trans = scales::log2_trans(), breaks = c(.2, .5,1 ,2,5))+
  scale_y_continuous(trans = scales::log2_trans(),breaks = c(.2, .5 ,2,5))+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0.2, 5), ylim = c(0.2, 5), xintercept = 0, yintercept = 0, ratio = 1)+
  labs(x = "Log2(IGR37xp/IGR37), Forward", y = "Log2(IGR37xp/IGR37), Reverse")+
  theme_classic()+ # increase x limit
  theme(legend.position = "none",
        # panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank())
ggsave(filename = file.path(base, 'figures/protein_dotplot.svg'), device = 'svg', plot = p, width = 3.5, height = 3.5, units = 'in')



# mitochondrial ETC PRM monitoring in shTRMU -----
# 1: H:shCtrl + L:TRMUsh2 (whole)
# 2: H:TRMUsh2 + L:shCtrl (whole)
# 3: H:shCtrl + L:TRMUsh3 (whole)
# 4: H:TRMUsh3 + L:shCtrl (whole)
# TRMU in whole proteome is not robust, change to mito fractionated TRMU (same cell, 1-4 whole, 5-8 mitofracitonation)
ptbl <- bind_rows(read_csv(file.path(base, 'data', 'SILAC_peptideRatio_241125.csv')) %>% 
                    filter(Peptide == 'TPNPDIVCNK') %>% 
                    mutate(Replicate = case_when(
                      Replicate == 5 ~ 1,
                      Replicate == 6 ~ 2,
                      Replicate == 7 ~ 3,
                      Replicate == 8 ~ 4
                    )),
                  read_csv(file.path(base, 'data', 'SILAC_peptideRatio_whole_241125.csv')) %>% filter(Peptide != 'TPNPDIVCNK')) %>% 
  separate_wider_regex('Protein Name', c(proteinName = ".*?", " |/", ".*")) %>% 
  dplyr::mutate(`Total Area` = as.numeric(`Total Area`)) %>% 
  dplyr::select(Peptide, `proteinName`, `Total Area`, `Isotope Label Type`, Replicate) %>% 
  pivot_wider(names_from = c(`Isotope Label Type`), values_from = `Total Area`) %>%
  dplyr::mutate(`shTRMU/shCtrl` = case_when(
    Replicate == 1 | Replicate == 3 ~ light / heavy,
    Replicate == 2 | Replicate == 4 ~ heavy / light
  )) %>% 
  dplyr::mutate(proteinName = fct_relevel(proteinName, 
                                          'MT-CO2', 'MT-ATP6', 'MT-CYB',
                                          'ETFA', 'ETFB', 'SDHA',  'TOMM40', 'ACTB', 'GAPDH', 'MTO1', 'TRMU'))
# reporting the table
wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', ptbl); openxlsx::openXL(wb)

# dotplot: forward vs reverse
ptbl %>% 
  select(-light, -heavy) %>% 
  mutate(labelType = ifelse(Replicate %in% c(1, 3), 'F', 'R')) %>% 
  group_by(proteinName, labelType) %>% 
  summarise(`k/c` = mean(`shTRMU/shCtrl`)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = c(labelType), values_from = `k/c`)%>%
  ggplot(aes(x = F, y = R))+
  geom_point(size = 0.5,  color = 'red')+
  ggrepel::geom_label_repel(aes(label = proteinName), min.segment.length = 0)+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0,3.5), ylim = c(0,3.5), xintercept = 1, yintercept = 1, ratio = 1)+
  labs(x = "forward = shTRMU(H)/shCtrl(L)", y = "reverse = shTRMU(L)/shCtrl(H)")+
  theme_classic() # increase x limit



# mitochondrial ETC PRM monitoring in resistance vs senistive -----
# F2: H:IGR37xp + L:IGR37 
# R2: L:IGR37xp + H:IGR37
ptbl <- read_csv(file.path(base, 'data', 'SILAC_peptideRatio_250205.csv')) %>% 
  separate_wider_regex('Protein Name', c(proteinName = ".*?", " |/", ".*")) %>% 
  mutate(`Total Area` = as.numeric(`Total Area`))%>% 
  dplyr::select(Peptide, `proteinName`, `Total Area`, `Isotope Label Type`, Replicate) %>% 
  pivot_wider(names_from = c(`Isotope Label Type`), values_from = `Total Area`)%>%
  mutate(`r/s` = case_when(
    Replicate == 'IV167_F2_MTendo' ~  heavy / light,
    Replicate == 'IV167_R2_MTendo' ~  light / heavy
  )) %>% 
  mutate(proteinName = fct_relevel(proteinName, 
                                   'MT-CO2', 'MT-ATP6', 'MT-ATP8', 
                                   'ETFA', 'ETFB', 'SDHA', 'MRPL11','ACTB', 'GAPDH', 'MTO1', 'TRMU')) %>% 
  filter(!(Replicate == 'IV167_R2_MTendo' &  proteinName == 'TRMU'))

# export the table 
wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', ptbl); openxlsx::openXL(wb)

# dotplot: forward vs reverse
ptbl %>% 
  dplyr::select(-light, -heavy) %>% 
  pivot_wider(names_from = c(Replicate), values_from = `r/s`)%>%
  ggplot(aes(x = IV167_F2_MTendo, y = IV167_R2_MTendo))+
  geom_point(size = 0.5,  color = 'red')+
  ggrepel::geom_label_repel(aes(label = proteinName), min.segment.length = 0)+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0,2.5), ylim = c(0,2.5), xintercept = 1, yintercept = 1, ratio = 1)+
  labs(x = "forward = IGR37xp(H)/IGR37(L)", y = "reverse = IGR37xp(L)/IGR37(H)")+
  theme_classic() # increase x limit


# TRMU KD: experiments from following dates 1)247030_shTRMU_IGR37xl 2)240806_shTRMU_IGR37xl 3) 241007_TRMUrescue -------
tbl <- readxl::read_excel(file.path(base, 'supplimentaryTable.xlsx'), sheet = 'S7.survival', range = 'R1C1:R126C9', na = '0') %>%  
  dplyr::filter(timepoint == 2) %>% 
  dplyr::select(-timepoint) %>% 
  pivot_longer(`...6`:`...9`, names_to = 'tr', values_to = 'relative survival') %>%
  group_by(date) %>% 
  mutate(`relative survival` = `relative survival` - mean(`relative survival`[nM == 'BLANK'], na.rm = TRUE),
         `relative survival` = ifelse(`relative survival` < 0, NA, `relative survival`)) %>% # remove incomplete data
  filter(nM != 'BLANK') %>% 
  filter(!is.na(`relative survival`)) %>%
  ungroup() %>% 
  group_by(plate) %>%
  mutate(`relative survival` = `relative survival` / mean(`relative survival`[nM == 1], na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(nM = log10(as.numeric(nM))) %>% 
  mutate(celltype = dplyr::case_when(
    celltype == 'IGR37xp_shCtrl' ~ 'shCtrl',
    celltype == 'IGR37xp_shTRMU' ~ 'shTRMU'))

 
mc <- drm(`relative survival` ~ nM, subset = celltype == 'shCtrl', data = tbl, fct = LN.4()) 
mk <- drm(`relative survival` ~ nM, subset = celltype == 'shTRMU', data = tbl, fct = LN.4()) 
newdata <- expand.grid(nM=c(log10(1000),log10(5000),log10(10000),log10(30000), log10(50000), log10(100000))) # all the confidence interval
cp <- cbind(newdata, predict(mc, newdata = newdata, interval = 'confidence'), celltype = 'shCtrl');names(cp)[2] <- 'relative survival'
kp <- cbind(newdata, predict(mk, newdata = newdata, interval = 'confidence'), celltype = 'shTRMU');names(kp)[2] <- 'relative survival'
ic50_ctrl <- signif(10^(ED(mc, 50, interval = "delta")[1])/1000, digits = 3)
ic50_kd <- signif(10^(ED(mk, 50, interval = "delta")[1])/1000, digits = 3)

global_linewidth <- 0.3
p <- tbl %>%
  ggplot(aes(x = nM, y = `relative survival`, color = celltype, shape = celltype))+
  geom_point(size = 0.1)+
  geom_smooth(method = drm, method.args = list(fct = LN.4()), se = FALSE, linewidth = global_linewidth)+
  geom_errorbar(data = cp, aes(ymin=cp[,'Lower'] , ymax=cp[,'Upper'] ), linewidth = global_linewidth, width = 0.05)+
  geom_errorbar(data = kp, aes(ymin=kp[,'Lower'] , ymax=kp[,'Upper'] ), linewidth = global_linewidth, width = 0.05)+
  scale_x_continuous(breaks = 0:5, labels = c('1 nM', '10 nM', '100 nM', '1 μM', '10 μM', '100 μM'))+
  scale_y_continuous(breaks = seq(0, 1, length.out = 5))+
  theme_classic()+
  coord_cartesian(ylim = c(0, 1))+
  xlab('[vem]')+ ylab('relative survival(%)')+
  # emf cannot take alpha, so i need to change the color to pink and light blue and convert to blue and red in svg. (easier to do for line than points)
  scale_color_manual(values = c("#cdd6ff", "#ffd6cf"), labels = paste0(c('shCtrl', 'shTRMU'),': ', c(ic50_ctrl, ic50_kd), ' μM'))+
  # scale_color_manual(values = c("#0433ff", "#ff2600"), labels = paste0(c('shCtrl', 'shTRMU'),': ', c(ic50_ctrl, ic50_kd), ' μM'))+
  scale_shape_manual(values = c(20, 15), labels = paste0(c('shCtrl', 'shTRMU'),': ', c(ic50_ctrl, ic50_kd), ' μM'), )+
  geom_hline(yintercept=0.5,linetype = "dashed", colour= "grey", linewidth = global_linewidth)+
  guides(color = guide_legend('Treatment (LD50)'), shape = guide_legend('Treatment (LD50)', override.aes = list(size = 2, alpha = 1)))+
  ggbreak::scale_x_break(c(0.2, 2.9))+
  theme(axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        text = element_text(size = 7, family = 'arial'),
        axis.line = element_line(colour = 'black', linewidth = global_linewidth),
        axis.ticks = element_line(colour = "black", linewidth = global_linewidth))
# plotly::ggplotly(p)
noleg <- p + theme(legend.position = "none")

ggsave(filename = file.path(base, 'figures/noleg.svg'), device = 'svg', plot = noleg, width = 2.7, height = 2, units = 'in')
ggsave(filename = file.path(base, 'figures/legend.svg'), device = 'svg', plot = p, width = 5, height = 5, units = 'in')
# subscript, importing legend





