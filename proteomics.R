# Library, global functions and parameters  ----
library(tidyverse)
# library(magrittr)
library(ggrepel)
# library(drc)
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

# # shortlist peptides based on mean protein ratios
# th <- 2
# gt_ls <- prottbl %>% filter(Fwd > th & Rvs > th) %>% pull(proteinName)
# gt_tbl <- ptbl %>% 
#   group_by(proteinName) %>% 
#   filter(proteinName %in% gt_ls) %>% 
#   ungroup()
# lt_ls <- prottbl %>% filter(Fwd < 1/th & Rvs < 1/th) %>% pull(proteinName)
# lt_tbl <- ptbl %>% 
#   group_by(proteinName) %>% 
#   filter(proteinName %in% lt_ls) %>% 
#   ungroup()
# wb <- openxlsx::createWorkbook(); 
# tn <- 'PRM_RWE_fulllist'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, ptbl);
# tn <- 'gt'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, gt_tbl); 
# tn <- 'lt'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, lt_tbl); 
# openxlsx::openXL(wb)

# output quantifiable protein arranged in desc order protein level fc
# pplist <- ptbl %>%
#   rowwise() %>% 
#   mutate(mean_peptidelvl = mean(c(Fwd, Rvs), na.rm = TRUE),
#     log2fc_peptidelvl = log2(mean_peptidelvl)) %>% 
#   ungroup() %>% 
#   group_by(proteinName) %>% # using `group_by` followed by `summarise` instead of `nest` to get two list-col for fwd and rvs respectively. 
#   mutate(across(matches('Fwd|Rvs'), ~list(.x))) %>% 
#   ungroup() %>%
#   rowwise() %>% 
#   mutate(mean_proteinlvl = mean(c(Fwd, Rvs), na.rm = TRUE),
#     log2fc_proteinlvl = log2(mean_proteinlvl)) %>% 
#   dplyr::arrange(desc(log2fc_proteinlvl)) %>% 
#   ungroup() %>%
#   group_by(proteinName) %>%
#   nest() %>%
#   ungroup() %>%
#   mutate(rnum = row_number()) %>%
#   mutate(shade = case_when(
#     rnum %% 2 == 1 ~ 0,
#     rnum %% 2 == 0 ~ 1,
#   )) %>% 
#   unnest(data) %>% 
#   dplyr::select(-Fwd, -Rvs)
# wb <- openxlsx::createWorkbook(); 
# openxlsx::addWorksheet(wb, 'test'); 
# lsty1 <-  openxlsx::createStyle(bgFill = "#d6d2d2")
# openxlsx::conditionalFormatting(wb, 'test', rows = 1:(nrow(pplist)+1), cols = 1:ncol(pplist), rule = paste0("$", LETTERS[which(colnames(pplist) == "shade")], "1==", 1), style = lsty1)
# openxlsx::writeData(wb, 'test', pplist); 
# openxlsx::openXL(wb)

# protein dot plot 
cf <- 1
gt_ls <- prottbl %>% filter(log2fc > cf) %>% na.omit()
lt_ls <- prottbl %>% filter(log2fc < -cf) %>% na.omit()
prottbl %>% 
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
ggsave(filename = file.path(base, 'figures/protein_dotplot.svg'), device = 'svg', plot = p, width = 3.5, height = 3.5, units = 'in')

# TODO split protein barplot into two
tb_num <- 10
prottbl %>% 
  na.omit() %>% 
  mutate(proteinName = fct_reorder(proteinName, log2fc)) %>%
  slice(1:tb_num, (n()-tb_num + 1):n()) %>% 
  select(proteinName, F1:R2) %>% 
  pivot_longer(-proteinName, names_to = 'replicates', values_to = 'ratios') %>% 
  mutate(fr = str_extract(replicates, '[:alpha:]'),
         tb = ratios > 1) %>% 
  ggplot(aes(proteinName, ratios))+
  stat_summary(aes(fill = tb), geom = 'bar', fun = 'mean' , na.rm = TRUE, width = .5)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), color = 'black', width = 0.2 )+
  geom_hline(yintercept=1, linetype="dashed",color = plot_pos, linewidth = .5)+
  geom_point(aes(color = fr), size = 2)+
  coord_flip()+
  ylab('Log2(IGR37xp / IGR37)')+
  scale_y_continuous(trans = scales::log2_trans(), breaks = c(0.25, 0.5, 1, 2, 4, 6))+
  scale_color_manual(values = hcl.colors(hcl.pals('diverging')[2], n = 8)[c(1,8)])+ # dot color
  # scale_fill_manual(values = )+
  guides(fill = FALSE, color = guide_legend(override.aes = list(size = 0.5)))+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.4))


# mitochondrial ETC PRM monitoring in shTRMU -----
# 1: H:shCtrl + L:TRMUsh2 (whole)
# 2: H:TRMUsh2 + L:shCtrl (whole)
# 3: H:shCtrl + L:TRMUsh3 (whole)
# 4: H:TRMUsh3 + L:shCtrl (whole)
# TRMU in whole proteome is not robust, change to mito fractionated TRMU (same cell, 1-4 whole, 5-8 mitofracitonation)
ptblkd <- bind_rows(read_csv(file.path(base, 'data', 'SILAC_peptideRatio_241125.csv')) %>% 
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
wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', ptblkd); openxlsx::openXL(wb)
# cratiotbl <- ptblkd %>% 
#   dplyr::select(Peptide, proteinName, `shTRMU/shCtrl`) %>% 
#   nest(data = -c(Peptide, proteinName)) %>% 
#   mutate(ratiopaste = map_chr(data, ~(paste0(signif(.x$`shTRMU/shCtrl`, digits = 3), collapse = ',') )))
# wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', cratiotbl); openxlsx::openXL(wb)

# dotplot: forward vs reverse
# ptblkd %>% 
#   dplyr::select(-light, -heavy) %>% 
#   mutate(labelType = ifelse(Replicate %in% c(1, 3), 'F', 'R')) %>% 
#   group_by(proteinName, labelType) %>% 
#   summarise(`k/c` = mean(`shTRMU/shCtrl`)) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = c(labelType), values_from = `k/c`)%>%
#   ggplot(aes(x = F, y = R))+
#   geom_point(size = 0.5,  color = 'red')+
#   ggrepel::geom_label_repel(aes(label = proteinName), min.segment.length = 0)+
#   ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0,3.5), ylim = c(0,3.5), xintercept = 1, yintercept = 1, ratio = 1)+
#   labs(x = "forward = shTRMU(H)/shCtrl(L)", y = "reverse = shTRMU(L)/shCtrl(H)")+
#   theme_classic() # increase x limit



# mitochondrial ETC PRM monitoring in resistance vs senistive -----
# F2: H:IGR37xp + L:IGR37 
# R2: L:IGR37xp + H:IGR37
ptblendo <- read_csv(file.path(base, 'data', 'SILAC_peptideRatio_250205.csv')) %>% 
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
wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', ptblendo); openxlsx::openXL(wb)
# cratiotbl <- ptblendo %>% 
#   dplyr::select(Peptide, proteinName, `r/s`) %>% 
#   nest(data = -c(Peptide, proteinName)) %>% 
#   mutate(ratiopaste = map_chr(data, ~(paste0(signif(.x$`r/s`, digits = 3), collapse = ',') )))
# wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', cratiotbl); openxlsx::openXL(wb)

# dotplot: forward vs reverse
# ptblendo %>% 
#   dplyr::select(-light, -heavy) %>% 
#   pivot_wider(names_from = c(Replicate), values_from = `r/s`)%>%
#   ggplot(aes(x = IV167_F2_MTendo, y = IV167_R2_MTendo))+
#   geom_point(size = 0.5,  color = 'red')+
#   ggrepel::geom_label_repel(aes(label = proteinName), min.segment.length = 0)+
#   ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0,2.5), ylim = c(0,2.5), xintercept = 1, yintercept = 1, ratio = 1)+
#   labs(x = "forward = IGR37xp(H)/IGR37(L)", y = "reverse = IGR37xp(L)/IGR37(H)")+
#   theme_classic() # increase x limit


# correlation plot between pblkd and pblendo ----
commprot <- ptblkd$proteinName[ptblkd$proteinName %in% ptblendo$proteinName] |> unique() |> as.character()
p <- bind_rows(
  ptblkd %>% 
  filter(proteinName %in% commprot) %>% 
  mutate(dataname = 'kd') %>% 
  mutate(ratios = `shTRMU/shCtrl`, .keep = 'unused') %>% 
  dplyr::select(proteinName, dataname, ratios), 
  ptblendo %>% 
  filter(proteinName %in% commprot) %>% 
  mutate(dataname = 'endo') %>% 
  mutate(ratios = `r/s`, .keep = 'unused') %>% 
  dplyr::select(proteinName, dataname, ratios)) %>% 
  pivot_wider(names_from = dataname, values_from = ratios, values_fn = ~ mean(.x, na.rm = TRUE)) %>% 
  ggplot(aes(x = kd, y = endo), shape = plot_shape)+
  geom_point(color = plot_pos, size = plot_size)+
  ggrepel::geom_text_repel(aes(label = proteinName), size= 3, force_pull = 10, force = 1, min.segment.length = 0.1, max.overlaps = 5)+
  scale_x_continuous(trans = scales::log2_trans(), breaks = c(.2, .5,1 ,2,5))+
  scale_y_continuous(trans = scales::log2_trans(),breaks = c(.2, .5 ,2,5))+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0.2, 5), ylim = c(0.2, 5), xintercept = 0, yintercept = 0, ratio = 1)+
  labs(x = "Log2(shTRMU / shCtrl)", y = "Log2(IGR37xp/IGR37)")+
  theme_classic()+ # increase x limit
  theme(legend.position = "none",
        # panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank())
ggsave(filename = file.path(base, 'figures/corr_shTRMUvsEndo.svg'), device = 'svg', plot = p, width = 3.5, height = 3.5, units = 'in')

# OXPHOS protein quantification -----
tbl <- readxl::read_excel(file.path(base, 'data', 'ETC_annotation.xlsx'), sheet = 'quantification', range = readxl::cell_rows(1:37)) %>% 
  mutate(LH = ifelse(str_detect(Compound, '\\((K8|R6)\\)$'), 'H', 'L')) 


# shTRMU / shCtrl in IGR37xp cells (runs included in the analysis: )
# before Correction: IV198_13_20250501221836, IV198_3_20250503145202, IV200_2, IV199_5, IV198_14_20250503005930,  IV200_3, IV199_6
# after correction:IV200_8, IV200_6, IV200_16, IV200_12, IV200_9, IV200_13
# after correction corresponding non corrected: IV198_13_20250501221836, IV198_3_20250503145202, IV200_2,IV199_5, IV198_14_20250503005930,IV199_6
ptbl <- tbl %>% 
  # select(proteinName, LH, IV198_13_20250501221836, IV198_3_20250503145202, IV200_2, IV199_5, IV198_14_20250503005930,  IV200_3, IV199_6) %>%  
  select(proteinName, LH, IV200_8, IV200_6, IV200_16, IV200_12, IV200_9, IV200_13) %>%
  pivot_longer(where(is.double), names_to = 'sample', values_to = 'amt') %>% 
  pivot_wider(names_from = 'LH', values_from = 'amt') %>%
  mutate(ratio = case_when(
    str_detect(sample, '(IV198_13|IV198_3|IV200_2|IV199_5|IV200_8|IV200_6|IV200_16|IV200_12)') ~ H / L,
    str_detect(sample, '(IV198_14|IV200_3|IV199_6|IV200_9|IV200_13)')  ~ L / H,
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


# ptbl %>% 
#   filter(sample %in% c('IV198_14', 'IV200_9')) %>% 
#   select(proteinName, sample, mratio) %>% 
#   pivot_wider( names_from = 'sample', values_from = 'mratio') %>% 
# ggplot(aes(x = IV198_14, y = IV200_9))+
# geom_point()+
# ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0, 3), ylim = c(0, 3),  ratio = 1)+
#   theme_classic()



# ratiotable
ratiotab <- ptbl %>%
  select(proteinName, sample, mratio, whereEncode) %>%
  pivot_wider(names_from = 'sample', values_from = 'mratio') %>%
  # select(proteinName, whereEncode, IV198_13, IV198_3, IV200_2,  IV199_5, IV198_14,  IV200_3, IV199_6) %>%
  select(proteinName, whereEncode, IV200_8, IV200_6, IV200_16, IV200_12, IV200_9, IV200_13) %>%
  arrange(proteinName)
# wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', ptbl); openxlsx::openXL(wb)
wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', ratiotab); openxlsx::openXL(wb)

# barplot
ptbl %>%
  filter(!sample %in% c('IV200_12', 'IV200_16', 'IV200_8')) %>% 
  # mutate(peptide_gn = sprintf(paste0('%s', '%',7-nchar(gn)+10, 's'), Sequence, gn) %>% factor(., levels = rev(.))) %>% 
  ggplot(aes(proteinName, mratio))+
  stat_summary(aes(fill = whereEncode), geom = 'bar', fun = 'mean' , na.rm = TRUE, width = .5)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), color = 'black', width = 0.2 )+
  geom_hline(yintercept=1, linetype="dashed",color = "red", linewidth = .5)+
  geom_point(aes(color = sample), size = 2)+
  coord_flip()+
  ylab('shTRMU / shCtrl in IGR37xp')+
  scale_color_manual(values = hcl.colors(hcl.pals('diverging')[2], n = 10)[c(1,3,6,10)])+ # dot color
  # scale_fill_manual(values = c(replication, repair, splicing))+
  guides(fill = 'none', color = guide_legend(override.aes = list(size = 0.5)))+
  theme_classic()+
  theme(axis.title.y = element_blank(), 
        legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.4))

# heatmap
ptbl %>% 
  filter(!sample %in% c('IV200_12', 'IV200_16', 'IV200_8')) %>% 
  ggplot(aes(x = sample, y = proteinName, fill = mratio)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("blue", "white", "red")) + # Choose your color scale
    theme_minimal() + # A clean theme
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate x-axis labels
    labs(title = "Gene Expression Heatmap",
         x = "Sample",
         y = "Gene",
         fill = "Expression Level")


# shTRMU / shCtrl in HEK293T ----
# before Correction: IV198_1, IV198_15_20250430191613, IV199_7, IV198_2_20250503113248, IV198_16_20250501015535, IV198_16, IV199_8
# after correction:IV200_4, IV200_10, IV200_14, IV200_5, IV200_11, IV200_15
ptbl <- tbl %>% 
  # select(proteinName, LH, IV198_1, IV198_15_20250430191613, IV199_7, IV198_2_20250503113248, IV198_16_20250501015535, IV198_16, IV199_8) %>% 
  select(proteinName, LH, IV200_4, IV200_10, IV200_14, IV200_5, IV200_11, IV200_15) %>%  
  pivot_longer(where(is.double), names_to = 'sample', values_to = 'amt') %>% 
  pivot_wider(names_from = 'LH', values_from = 'amt') %>%
  mutate(ratio = case_when(
    str_detect(sample, '(IV198_1|IV198_15|IV199_7|IV200_4|IV200_10|IV200_14)') ~ H / L,
    str_detect(sample, '(IV198_2|IV198_16|IV199_8|IV200_5|IV200_11|IV200_15)')  ~ L / H,
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

# ratiotable
ratiotab <- ptbl %>%
  select(proteinName, sample, mratio, whereEncode) %>%
  pivot_wider(names_from = 'sample', values_from = 'mratio') %>%
  # select(proteinName, whereEncode, IV198_1, IV198_15, IV199_7, IV198_2, IV198_16, IV199_8) %>%
  select(proteinName, whereEncode, IV200_4, IV200_10, IV200_14, IV200_5, IV200_11, IV200_15) %>%
  arrange(proteinName)
# wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', ptbl); openxlsx::openXL(wb)
wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', ratiotab); openxlsx::openXL(wb)

ptbl %>% 
  filter(!sample %in% c('IV200_10', 'IV200_14', 'IV200_5')) %>% 
  # mutate(peptide_gn = sprintf(paste0('%s', '%',7-nchar(gn)+10, 's'), Sequence, gn) %>% factor(., levels = rev(.))) %>% 
  ggplot(aes(proteinName, mratio))+
  stat_summary(aes(fill = whereEncode), geom = 'bar', fun = 'mean' , na.rm = TRUE, width = .5)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), color = 'black', width = 0.2 )+
  geom_hline(yintercept=1, linetype="dashed",color = "red", linewidth = .5)+
  geom_point(aes(color = sample), size = 2)+
  coord_flip()+
  ylab('shTRMU / shCtrl in HEK293T')+
  scale_color_manual(values = hcl.colors(hcl.pals('diverging')[2], n = 10)[c(1,3, 6,10)])+ # dot color
  # scale_fill_manual(values = c(replication, repair, splicing))+
  guides(fill = FALSE, color = guide_legend(override.aes = list(size = 0.5)))+
  theme_classic()+
  theme(axis.title.y = element_blank(), 
        legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.4))

# MTO1 KO / IGR37xp (IV198_8	IV198_12, IV198_7_20250505001525) -----
ptbl <- tbl %>% 
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

# ratiotable 
ratiotab <- ptbl %>% 
  select(proteinName, sample, mratio, whereEncode) %>% 
  pivot_wider(names_from = 'sample', values_from = 'mratio') %>% 
  arrange(proteinName)
wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', ratiotab); openxlsx::openXL(wb)

ptbl %>% 
  # mutate(peptide_gn = sprintf(paste0('%s', '%',7-nchar(gn)+10, 's'), Sequence, gn) %>% factor(., levels = rev(.))) %>% 
  ggplot(aes(proteinName, mratio))+
  stat_summary(aes(fill = whereEncode), geom = 'bar', fun = 'mean' , na.rm = TRUE, width = .5)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), color = 'black', width = 0.2 )+
  geom_hline(yintercept=1, linetype="dashed",color = "red", linewidth = .5)+
  geom_point(aes(color = sample), size = 2)+
  coord_flip()+
  ylab('MTO1 KO / IGR37xp')+
  scale_color_manual(values = hcl.colors(hcl.pals('diverging')[2], n = 8)[c(1,5,8)])+ # dot color
  # scale_fill_manual(values = c(replication, repair, splicing))+
  guides(fill = FALSE, color = guide_legend(override.aes = list(size = 0.5)))+
  theme_classic()+
  theme(axis.title.y = element_blank(), 
        legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.4))

# DONE IGR37xp / IGR37 (runs included in the analysis: IV198_5_20250502113912, IV198_6_20250502145929, IV198_9) ----
ptbl <- tbl %>% 
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
  mutate(whereEncode = ifelse(str_detect(proteinName, 'MT-'), 'MT', 'nuc')) %>% 
  dplyr::mutate(proteinName = fct_relevel(proteinName, 
                                          'TRMU', 'MTO1', 'ACTB', 'GAPDH',  #4
                                          'ETFA', 'SDHA', 'MRPL11' ,'TOMM40', #4
                                          'MT-ND6', 'MT-ND5', 'MT-ND4', 'MT-ND3', #4
                                          'MT-CYB', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'MT-CO2',#5
  ))

# ratiotable 
ratiotab <- ptbl %>% 
  select(proteinName, sample, mratio, whereEncode) %>% 
  pivot_wider(names_from = 'sample', values_from = 'mratio') %>% 
  arrange(proteinName)
wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', ratiotab); openxlsx::openXL(wb)

# barplot
ptbl %>% 
  # filter(!proteinName %in% c('TRMU', 'MTO1', 'ACTB')) %>% 
  # mutate(peptide_gn = sprintf(paste0('%s', '%',7-nchar(gn)+10, 's'), Sequence, gn) %>% factor(., levels = rev(.))) %>% 
  ggplot(aes(proteinName, mratio))+
  stat_summary(aes(fill = whereEncode), geom = 'bar', fun = 'mean' , na.rm = TRUE, width = .5)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), color = 'black', width = 0.2 )+
  geom_hline(yintercept=1, linetype="dashed",color = "red", linewidth = .5)+
  geom_point(aes(color = sample), size = 2)+
  coord_flip()+
  ylab('IGR37xp / IGR37')+
  scale_color_manual(values = hcl.colors(hcl.pals('diverging')[2], n = 8)[c(1,2,7,8)])+ # dot color
  # scale_fill_manual(values = )+
  guides(fill = FALSE, color = guide_legend(override.aes = list(size = 0.5)))+
  theme_classic()+
  theme(axis.title.y = element_blank(), 
        legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.4))
