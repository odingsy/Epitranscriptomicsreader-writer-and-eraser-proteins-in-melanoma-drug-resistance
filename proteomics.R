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
wb <- openxlsx::createWorkbook(); 
tn <- 'PRM_RWE_fulllist'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, ptbl);
tn <- 'gt'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, gt_tbl); 
tn <- 'lt'; openxlsx::addWorksheet(wb, tn); openxlsx::writeData(wb, tn, lt_tbl); 
openxlsx::openXL(wb)

# output quantifiable protein arranged in desc order protein level fc
pplist <- ptbl %>%
  rowwise() %>% 
  mutate(mean_peptidelvl = mean(c(Fwd, Rvs), na.rm = TRUE),
    log2fc_peptidelvl = log2(mean_peptidelvl)) %>% 
  ungroup() %>% 
  group_by(proteinName) %>% # using `group_by` followed by `summarise` instead of `nest` to get two list-col for fwd and rvs respectively. 
  mutate(across(matches('Fwd|Rvs'), ~list(.x))) %>% 
  ungroup() %>%
  rowwise() %>% 
  mutate(mean_proteinlvl = mean(c(Fwd, Rvs), na.rm = TRUE),
    log2fc_proteinlvl = log2(mean_proteinlvl)) %>% 
  dplyr::arrange(desc(log2fc_proteinlvl)) %>% 
  ungroup() %>%
  group_by(proteinName) %>%
  nest() %>%
  ungroup() %>%
  mutate(rnum = row_number()) %>%
  mutate(shade = case_when(
    rnum %% 2 == 1 ~ 0,
    rnum %% 2 == 0 ~ 1,
  )) %>% 
  unnest(data) %>% 
  dplyr::select(-Fwd, -Rvs)

wb <- openxlsx::createWorkbook(); 
openxlsx::addWorksheet(wb, 'test'); 
lsty1 <-  openxlsx::createStyle(bgFill = "#d6d2d2")
openxlsx::conditionalFormatting(wb, 'test', rows = 1:(nrow(pplist)+1), cols = 1:ncol(pplist), rule = paste0("$", LETTERS[which(colnames(pplist) == "shade")], "1==", 1), style = lsty1)
openxlsx::writeData(wb, 'test', pplist); 
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


  