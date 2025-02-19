# The code processes the outputs from skyline 

library(tidyverse)
library(magrittr)
library(ggrepel)


# global functions and parameters  ----
base <- '/Users/shiyuanguo/Library/CloudStorage/GoogleDrive-sguo039@ucr.edu/My Drive/PhD_study/IGR_RWE/1.IGR37xp_manuscript'

plot_shape = 20; plot_size = 0.5; plot_alpha = 0.4; plot_pos = '#F8766D'# ggploting parameters

mm <- list(  
  Mean = ~.x/mean(cur_data()[[cur_column()]], na.rm = TRUE), 
  Median = ~.x/median(cur_data()[[cur_column()]], na.rm = TRUE)
) 
rsd = ~sd(c(.x, .y), na.rm = TRUE)/mean(c(.x, .y),na.rm = TRUE) # calculating the RSD
env <- environment() # assign format to obtain the intermediate dataset. 
lm_equation <- function(y, x, data, ...){
  # annotating the equation on ggplot. 
  formula <- as.formula(paste(y,paste(x, collapse = ' + '),sep = ' ~ '))
  m <- lm(formula, data, ...) 
  substituteList <- list(a = format(abs(unname(coef(m)[1])), digits = 3),
                         b = format(unname(coef(m)[2]), digits = 3),
                         r2 = format(summary(m)$r.squared, digits = 3))
  if (coef(m)[1] > 0){
    eq <- substitute(italic(y) == b %.% italic(x)+a*";"~~italic(r)^2~"="~r2, substituteList)
  } else {
    eq <- substitute(italic(y) == b %.% italic(x)-a*";"~~italic(r)^2~"="~r2, substituteList)
  }
  as.character(as.expression(eq))
}



# obtain the injection info from both replicates ----
injectionfiles = paste0(file.path(base, 'data', 'injectionLists', '12min_max40txsits'), '_000', 1:5, '.csv')
injectionList <- lapply(injectionfiles, function(n){
  readr::read_csv(n) %>%
    dplyr::mutate(Comment = stringr::str_remove_all(Comment,' \\((light|heavy)\\)') %>% 
                    stringr::str_remove_all('\\[[:graph:][:digit:]{2}\\.[:digit:]{6}\\]')) %>%
    .[['Comment']] %>% 
    unique(.)
})


# loading SILAC_peptidesRatio between vemurafenib resistance vs sensitive cell ----
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
  ############
  mutate(`Total Area` = as.numeric(`Total Area`))%>% 
  select(Peptide, `Protein Name`, fr, `Total Area`, `Isotope Label Type`) %>% 
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

th <- 2
gt_tbl <- ptbl %>% 
  group_by(proteinName) %>% 
  filter(proteinName %in% proteinName[Fwd > th & Rvs > th])
lt_tbl <- ptbl %>% 
  group_by(proteinName) %>% 
  filter(proteinName %in% proteinName[Fwd < 1/th & Rvs < 1/th]) 
wb <- openxlsx::createWorkbook(); 
openxlsx::addWorksheet(wb, 'gt'); openxlsx::writeData(wb, 'gt', gt_tbl); 
openxlsx::addWorksheet(wb, 'lt'); openxlsx::writeData(wb, 'lt', lt_tbl); 
openxlsx::openXL(wb)
# wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, 'test'); openxlsx::writeData(wb, 'test', ptbl); openxlsx::openXL(wb)

# peptide level plot ----
# TODO change to log scale
cutoff <- 2
lab_cutoff <- 2
ptbl %>% 
  ggplot(aes(x = Fwd, y = Rvs), shape = plot_shape)+
  geom_point(color = 'grey', size = plot_size, alpha = plot_alpha)+
  geom_point(data = . %>% filter((Fwd > cutoff & Rvs > cutoff) | (Fwd < 1/cutoff & Rvs < 1/cutoff)), color = plot_pos, size = plot_size)+
  ggrepel::geom_text_repel(data = . %>% filter(Fwd > lab_cutoff & Rvs > lab_cutoff), aes(label = proteinName), size= 3, force_pull = 10, force = 1, min.segment.length = 0.1, max.overlaps = Inf)+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0, 4), ylim = c(0, 4), xintercept = 1, yintercept = 1, ratio = 1)+
  labs(x = "Ratio(Resistant/Sensitive), Forward", y = "Ratio(Resistance/Sensitive), Reverse")+
  theme_classic()+ # increase x limit
  theme(legend.position = "none")


# protein level (w\o normalization) ----
prottbl <- ptbl %>%
  group_by(proteinName) %>%
  summarise(across(matches('Fwd|Rvs'), ~list(.x))) %>%
  ungroup() %>%
  mutate(rsd = map2_dbl(Fwd, Rvs, rsd)) %>% 
  mutate(across(matches('Fwd|Rvs'), ~map_dbl(., ~mean(.x, na.rm = TRUE)))) %>% 
  # filter(lengths(Fwd_normByMedian) > 1) %>%
  # mutate(rsd_scaled = ifelse(rsd > 0.5, 1, (rsd-min(rsd, na.rm = TRUE)) / (0.5 - min(rsd, na.rm = TRUE)))) %>%
  rowwise() %>% 
  mutate(log2fc = log2(mean(Fwd, Rvs, na.rm = TRUE))) %>% 
  arrange(desc(log2fc))

# protein plot 
lab_pcut <- pcut <- 2
lab_ncut <- ncut <- 2.5
gt_protein <- gt_tbl %>% 
  summarise(ratioByPeptide = base::mean(c(Fwd, Rvs)))
lt_protein <- lt_tbl %>% 
  summarise(ratioByPeptide = base::mean(c(Fwd, Rvs)))
prottbl %>% 
  ggplot(aes(x = Fwd, y = Rvs), shape = plot_shape)+
  geom_point(color = 'grey', size = plot_size, alpha = plot_alpha)+
  geom_point(data = . %>% filter(proteinName %in% c(gt_protein$proteinName, lt_protein$proteinName)), color = plot_pos, size = plot_size)+
  ggrepel::geom_text_repel(data = . %>% filter(proteinName %in% c(gt_protein$proteinName, lt_protein$proteinName[rank(lt_protein$ratioByPeptide) < 6])), aes(label = proteinName), size= 3, force_pull = 10, force = 1, min.segment.length = 0.1, max.overlaps = Inf)+
  scale_x_continuous(trans = scales::log2_trans(), breaks = c(.1, .2, .5,1 ,2,5,10))+
  scale_y_continuous(trans = scales::log2_trans(),breaks = c(.1, .2, .5 ,2,5,10))+
  ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0.1, 12), ylim = c(0.1, 12), xintercept = 0, yintercept = 0, ratio = 1)+
  labs(x = "Ratio(Resistant/Sensitive), Forward", y = "Ratio(Resistance/Sensitive), Reverse")+
  theme_classic()+ # increase x limit
  theme(legend.position = "none",
        # panel.border = element_rect(colour = "black", fill=NA),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))


  
# GSEA
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

prottbl[prottbl$proteinName == 'HuR', 'proteinName'] <- 'ELAVL1'

# toEndb <- select(org.Hs.eg.db, keys = prottbl$proteinName, keytype = 'SYMBOL', columns = 'ENSEMBL')
ctbl <- clusterProfiler::bitr(geneID = prottbl$proteinName, fromType = 'ALIAS', toType = 'ENSEMBL', OrgDb = org.Hs.eg.db, drop = FALSE)
# s <- clusterProfiler::bitr(geneID = prottbl$proteinName, fromType = 'SYMBOL', toType = 'ENSEMBL', OrgDb = org.Hs.eg.db, drop = FALSE)
# wb <- openxlsx::createWorkbook();
# openxlsx::addWorksheet(wb, 'a'); openxlsx::writeData(wb, 'a', a);
# openxlsx::addWorksheet(wb, 's'); openxlsx::writeData(wb, 's', s);
# openxlsx::openXL(wb)
gtbl <- left_join(prottbl, ctbl, by = c('proteinName' = 'ALIAS'))

gl <- gtbl$log2fc
names(gl) <- gtbl$ENSEMBL


gse <- clusterProfiler::gseGO(geneList=gl,
                              ont ="BP",
                              keyType = "ENSEMBL",
                              minGSSize = 3,
                              maxGSSize = 43,
                              pvalueCutoff = 0.05,
                              verbose = TRUE,
                              OrgDb = org.Hs.eg.db,
                              pAdjustMethod = "none")
enrichplot::dotplot(gse, showCategory=10, split=".sign", label_format = 50) + facet_grid(.~.sign)


#
library(msigdbr)
library(enrichplot)
msigdbr::msigdbr_species()
msigdbr::msigdbr_collections()
cgp <- msigdbr::msigdbr(species = 'human', category = "C2")
msigdbr_t2g = cgp %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
prottbl_filtered <- prottbl %>% filter(Fwd > pcut & Rvs > pcut | Fwd < 1/ncut & Rvs < 1/ncut)

et <- clusterProfiler::enricher(gene = prottbl_filtered$proteinName, TERM2GENE = msigdbr_t2g)
enrichplot::dotplot(et, showCategory=10) #+ facet_grid(.~.sign)
cnetplot(et)
ridgeplot(et)

mettl9 <- cgp %>% filter(gene_symbol == 'METTL9')
trmu <- cgp %>% filter(gene_symbol == 'TRMU')
pus7 <- cgp %>% filter(gene_symbol == 'PUS7')
mto1 <- cgp %>% filter(gene_symbol == 'MTO1')
wb <- openxlsx::createWorkbook();
openxlsx::addWorksheet(wb, 'mettl9'); openxlsx::writeData(wb, 'mettl9', mettl9);
openxlsx::addWorksheet(wb, 'trmu'); openxlsx::writeData(wb, 'trmu', trmu);
openxlsx::addWorksheet(wb, 'pus7'); openxlsx::writeData(wb, 'pus7', pus7);
openxlsx::addWorksheet(wb, 'mto1'); openxlsx::writeData(wb, 'mto1', mto1);
openxlsx::openXL(wb)
  

  
  
  





# protein Level using normalized ratio
# mean_Fwd mean_Rvs
# 0.843     1.18
# ptbl %>% # summarise(mean_Fwd = mean(Fwd, na.rm = TRUE), mean_Rvs = mean(Rvs, na.rm = TRUE)) 
#   # filter((Fwd < 2 & Rvs < 2)|(Fwd > 0.5 & Rvs > 0.5) ) %>% # summarise(mean_Fwd = mean(Fwd, na.rm = TRUE), mean_Rvs = mean(Rvs, na.rm = TRUE)) 
#   mutate(across(c(Fwd,Rvs), mm, .names = '{.col}_normBy{.fn}')) %>% 
#   # mutate(rsd = map2_dbl(Fwd_normByMean, Rvs_normByMean, rsd)) %T>% assign('t1',., envir = env) %>% 
#   # filter(rsd < 0.2) %T>% assign('t2',., envir = env) %>% 
#   group_by(`Protein Gene`) %>% 
#   summarise(across(contains('normByMean'), ~list(.x))) %>% 
#   ungroup() %>% 
#   filter(!is.na(`Protein Gene`)) %>% View() # 13 peptides does not belong to a protein
#   mutate(rsd = map2_dbl(Fwd_normByMean, Rvs_normByMean, rsd)) %T>% assign('t3',., envir = env) %>% 
#   #filter(lengths(Fwd_normByMedian) > 1) %>% 
#   mutate(across(contains('normByMean'), ~map_dbl(., ~mean(.x, na.rm = TRUE)))) %>% 
#   mutate(rsd_scaled = ifelse(rsd > 0.5, 1, (rsd-min(rsd, na.rm = TRUE)) / (0.5 - min(rsd, na.rm = TRUE)))) %>% 
#   separate(col = 'Protein Gene', into = c('proteinName', 'id'), sep = " ") %>% 
#   separate(col = 'proteinName', into = c('proteinName', 'xx'), sep = '/' ) %>% 
#   dplyr::select(-c(id, xx)) %T>% assign('t4',., envir = env)
# 
# t4 %>% 
#   ggplot(aes(x = Fwd_normByMean, y = Rvs_normByMean))+
#   geom_point(data = . %>% filter((Fwd_normByMean > 1/cutoff & Fwd_normByMean < cutoff) | (Rvs_normByMean > 1/cutoff & Rvs_normByMean < cutoff)), color = 'grey', size = 1)+
#   geom_point(data = . %>% filter((Fwd_normByMean > cutoff & Rvs_normByMean > cutoff) | (Fwd_normByMean < 1/cutoff & Rvs_normByMean < 1/cutoff)), aes(alpha = 1-rsd_scaled), color = 'red', size = 1)+
#   ggrepel::geom_label_repel(data = . %>% filter(Fwd_normByMean < 3 & Rvs_normByMean < 3 & Fwd_normByMean > cutoff & Rvs_normByMean > cutoff), aes(label = proteinName), box.padding = 0.5, max.overlaps = Inf)+
#   # geom_point(data = . %>% filter(nchar(zincFinger) !=0 & `num<1.5` <= 2 & mean > 1.5) %>% filter(gn %in% c('POLD1')), aes(color = `Arsenite replaceable zinc-finger domain`), size = 1)+
#   ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0, 3), ylim = c(0, 3), xintercept = 1, yintercept = 1, ratio = 1)+
#   labs(x = "Fwd = Resistant(H) / WT(L)", y = "1/Rvs = Resistant(L) / WT(H)")+
#   theme_classic()+ # increase x limit
#   theme(legend.position = "none")
# 
# t4 %>% 
#   ggplot(aes(x = Fwd_normByMean, y = Rvs_normByMean))+
#   geom_point(data = . %>% filter((Fwd_normByMean > 1/cutoff & Fwd_normByMean < cutoff) | (Rvs_normByMean > 1/cutoff & Rvs_normByMean < cutoff)), color = 'grey', size = 1)+
#   geom_point(data = . %>% filter((Fwd_normByMean > cutoff & Rvs_normByMean > cutoff) | (Fwd_normByMean < 1/cutoff & Rvs_normByMean < 1/cutoff)), aes(alpha = 1-rsd_scaled), color = 'red', size = 1)+
#   ggrepel::geom_label_repel(data = . %>% filter(Fwd_normByMean < 1/2.5 & Rvs_normByMean < 1/2.5), aes(label = proteinName), box.padding = 0.5, max.overlaps = Inf)+
#   # geom_point(data = . %>% filter(nchar(zincFinger) !=0 & `num<1.5` <= 2 & mean > 1.5) %>% filter(gn %in% c('POLD1')), aes(color = `Arsenite replaceable zinc-finger domain`), size = 1)+
#   ggh4x::coord_axes_inside(labels_inside = TRUE, xlim = c(0, 0.7), ylim = c(0, 0.7), xintercept = 1, yintercept = 1, ratio = 1)+
#   labs(x = "Fwd = Resistant(H) / WT(L)", y = "1/Rvs = Resistant(L) / WT(H)")+
#   theme_classic()+ # increase x limit
#   theme(legend.position = "none")



