# Library, global functions and parameters  ----
library(tidyverse)
library(drc)
library(drda)
library(scales)
base <- '/Users/shiyuanguo/Library/CloudStorage/GoogleDrive-sguo039@ucr.edu/My Drive/PhD_study/IGR_RWE/1.IGR37xp_manuscript'
global_linewidth <- 0.3 # ggploting parameters


# TRMU KD: experiments from following dates 1)247030_shTRMU_IGR37xl 2)240806_shTRMU_IGR37xl 3) 241007_TRMUrescue -------
tbl1 <- readxl::read_excel(file.path(base, 'data', 'phenotypicAssays.xlsx'), sheet = 'shTRMU survival', range = 'R1C1:R126C9', na = '0') %>%  
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
  mutate(nM = as.numeric(nM)) %>% 
  mutate(celltype = dplyr::case_when(
    celltype == 'IGR37xp_shCtrl' ~ 'shCtrl',
    celltype == 'IGR37xp_shTRMU' ~ 'shTRMU1'))

tbl23 <- readxl::read_excel(file.path(base, 'data', 'phenotypicAssays.xlsx'), sheet = 'shTRMU otherSeq', range = 'R1C1:R75C7', na = '0') %>%  
  dplyr::filter(timepoint == 2) %>% 
  dplyr::select(-timepoint) %>% 
  pivot_longer(`...4`:`...7`, names_to = 'tr', values_to = 'relative survival') %>%
  filter(!is.na(`relative survival`)) %>%
  mutate(`relative survival` = `relative survival` - mean(`relative survival`[nM == 'BLANK'], na.rm = TRUE),
         `relative survival` = ifelse(`relative survival` < 0, NA, `relative survival`)) %>% # remove incomplete data
  filter(nM != 'BLANK') %>% 
  filter(!is.na(`relative survival`)) %>%
  group_by(celltype) %>%
  mutate(`relative survival` = `relative survival` / mean(`relative survival`[nM == 1], na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(nM = as.numeric(nM)) %>% 
  mutate(celltype = dplyr::case_when(
    celltype == 'IGR37xl_shctrl' ~ 'shCtrl',
    celltype == 'IGR37xl_shTRMUsh2' ~ 'shTRMU2',
    celltype == 'IGR37xl_shTRMUsh3' ~ 'shTRMU3'))


# LC50 of each shTRMU1 vs LC50 of each shCtrl 
lc50_cal<- function(shctrl, shtrt){
  mc <- drda::drda(`relative survival` ~ nM,  data = shctrl)
  mk <-drda::drda(`relative survival` ~ nM,  data = shtrt)
  fittedctrl <- shctrl$nM[as.numeric(names(mc$fitted.values))] 
  fittedtrt <- shtrt$nM[as.numeric(names(mk$fitted.values))] 
  lc50_c <- approx(x = mc$fitted.values, y = fittedctrl, xout = 0.5)$y
  lc50_t <- approx(x = mk$fitted.values, y = fittedtrt, xout = 0.5)$y
  return(c(lc50_c, lc50_t))
}

shctrl <- tbl1 %>%  
  filter(date == '240730') %>% 
  filter(celltype== 'shCtrl')
shtrt <- tbl1 %>% 
  filter(date == '240730') %>% 
  filter(celltype == 'shTRMU1')
lc50_cal(shctrl, shtrt)

shctrl <- tbl1 %>%  
  filter(date == '240806') %>% 
  filter(celltype== 'shCtrl')
shtrt <- tbl1 %>% 
  filter(date == '240806') %>% 
  filter(celltype == 'shTRMU1')
lc50_cal(shctrl, shtrt)

shctrl <- tbl1 %>%  
  filter(date == '241007') %>% 
  filter(celltype== 'shCtrl')
shtrt <- tbl1 %>% 
  filter(date == '241007') %>% 
  filter(celltype == 'shTRMU1')
lc50_cal(shctrl, shtrt)

shctrl <- tbl1 %>%  
  filter(celltype== 'shCtrl')
shtrt <- tbl23 %>% 
  filter(celltype == 'shTRMU2')
lc50_cal(shctrl, shtrt)


# LC50 shTRMU1 vs LC50 shCtrl (need log10 transformed)
tbl <- bind_rows(tbl1 %>% dplyr::select(celltype, nM, tr, `relative survival`),tbl23 %>% filter(celltype == 'shTRMU2')) %>% 
  mutate(nM = log10(nM)) 
mc <- drm(`relative survival` ~ nM, subset = celltype == 'shCtrl', data = tbl, fct = LN.4())
mk1 <- drm(`relative survival` ~ nM, subset = celltype == 'shTRMU1', data = tbl, fct = LN.4())
mk2 <- drm(`relative survival` ~ nM, subset = celltype == 'shTRMU2', data = tbl, fct = LN.4())
newdata <- expand.grid(nM=c(log10(1000),log10(5000),log10(10000),log10(30000), log10(50000), log10(100000))) # all the confidence interval
cp <- cbind(newdata, predict(mc, newdata = newdata, interval = 'confidence'), celltype = 'shCtrl');names(cp)[2] <- 'relative survival'
kp1 <- cbind(newdata, predict(mk1, newdata = newdata, interval = 'confidence'), celltype = 'shTRMU1');names(kp1)[2] <- 'relative survival'
kp2 <- cbind(newdata, predict(mk2, newdata = newdata, interval = 'confidence'), celltype = 'shTRMU2');names(kp2)[2] <- 'relative survival'
ic50_ctrl <- signif(10^(ED(mc, 50, interval = "delta")[1])/1000, digits = 3)
ic50_kd1 <- signif(10^(ED(mk1, 50, interval = "delta")[1])/1000, digits = 3)
ic50_kd2 <- signif(10^(ED(mk2, 50, interval = "delta")[1])/1000, digits = 3)

# survival plot 
p <- tbl %>% 
  ggplot(aes(x = nM, y = `relative survival`, color = celltype, shape = celltype))+
  geom_point(size = 0.1)+
  geom_smooth(method = drm, method.args = list(fct = LN.4()), se = FALSE, linewidth = global_linewidth)+
  geom_errorbar(data = cp, aes(ymin=cp[,'Lower'] , ymax=cp[,'Upper'] ), linewidth = global_linewidth, width = 0.05)+
  geom_errorbar(data = kp1, aes(ymin=kp1[,'Lower'] , ymax=kp1[,'Upper'] ), linewidth = global_linewidth, width = 0.05)+
  geom_errorbar(data = kp2, aes(ymin=kp2[,'Lower'] , ymax=kp2[,'Upper'] ), linewidth = global_linewidth, width = 0.05)+
  scale_x_continuous(breaks = 0:5, labels = c('1 nM', '10 nM', '100 nM', '1 μM', '10 μM', '100 μM'))+
  scale_y_continuous(breaks = seq(0, 1, length.out = 5))+
  theme_classic()+
  coord_cartesian(ylim = c(0, 1))+
  xlab('[vem]')+ ylab('relative survival(%)')+
  # emf cannot take alpha, so i need to change the color to pink and light blue and convert to blue and red in svg. (easier to do for line than points)
  scale_color_manual(values = c("#cdd6ff", "#f1bfd6", "#f6cda9"))+
  # scale_color_manual(values = c("#0433ff",  "#ff0074", "#ff7800"))+
  scale_shape_manual(values = c(20, 15, 17))+
  geom_hline(yintercept=0.5,linetype = "dashed", colour= "grey", linewidth = global_linewidth)+
  guides(color = guide_legend(title="Treatment"), shape = guide_legend('Treatment', override.aes = list(size = 2, alpha = 1)))+
  # ggbreak::scale_x_break(c(0.2, 2.9))+
  theme(axis.text.x.top = element_blank(),
        legend.position = c(0.2, 0.2),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        # text = element_text(size = 7, family = 'arial'),
        axis.line = element_line(colour = 'black', linewidth = global_linewidth),
        axis.ticks = element_line(colour = "black", linewidth = global_linewidth))
# plotly::ggplotly(p)
# noleg <- p + theme(legend.position = "none")

# ggsave(filename = file.path(base, 'figures/noleg.svg'), device = 'svg', plot = noleg, width = 2.7, height = 2, units = 'in')
ggsave(filename = file.path(base, 'figures/legend.pdf'), device = 'pdf', plot = p, width = 3, height = 2, units = 'in')
# subscript, importing legend


# NAC survival ----
tbl <- readxl::read_excel(file.path(base,'data', 'phenotypicAssays.xlsx'), sheet = 'shTRMU NAC survival', range = 'R1C1:R506C9', na = '0') %>%  
  dplyr::filter(timepoint == 3) %>% 
  dplyr::select(-timepoint) %>% 
  pivot_longer(`...4`:`...9`, names_to = 'tr', values_to = 'relative survival') %>%
  mutate(`relative survival` = `relative survival` - mean(`relative survival`[nM == 'BLANK'], na.rm = TRUE),
         `relative survival` = ifelse(`relative survival` < 0, NA, `relative survival`)) %>% # remove incomplete data
  filter(nM != 'BLANK') %>% 
  filter(!is.na(`relative survival`)) %>%
  group_by(celltype) %>%
  mutate(`relative survival` = `relative survival` / mean(`relative survival`[nM == 1], na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(nM = log10(as.numeric(nM))) %>% 
  mutate(celltype = dplyr::case_when(
    celltype == 'IGR37xl_EV' ~ 'shCtrl',
    celltype == 'IGR37xl_EV_NAC' ~ 'shCtrl_NAC',
    celltype == 'IGR37xl_shTRMU' ~ 'shTRMU',
    celltype == 'IGR37xl_shTRMU_NAC' ~ 'shTRMU_NAC'))


mc <- drm(`relative survival` ~ nM, subset = celltype == 'shCtrl', data = tbl, fct = LN.4())
mcn <- drm(`relative survival` ~ nM, subset = celltype == 'shCtrl_NAC', data = tbl, fct = LN.4()) 
mk <- drm(`relative survival` ~ nM, subset = celltype == 'shTRMU', data = tbl, fct = LN.4())
mkn <- drm(`relative survival` ~ nM, subset = celltype == 'shTRMU_NAC', data = tbl, fct = LN.4()) 
newdata <- expand.grid(nM=c(log10(1000),log10(5000),log10(10000),log10(30000), log10(50000))) # all the confidence interval
cp <- cbind(newdata, predict(mc, newdata = newdata, interval = 'confidence'), celltype = 'shCtrl');names(cp)[2] <- 'relative survival'
cpn <- cbind(newdata, predict(mcn, newdata = newdata, interval = 'confidence'), celltype = 'shCtrl_NAC');names(cpn)[2] <- 'relative survival'
kp <- cbind(newdata, predict(mk, newdata = newdata, interval = 'confidence'), celltype = 'shTRMU');names(kp)[2] <- 'relative survival'
kpn <- cbind(newdata, predict(mkn, newdata = newdata, interval = 'confidence'), celltype = 'shTRMU_NAC');names(kpn)[2] <- 'relative survival'
ic50_ctrl <- signif(10^(ED(mc, 50, interval = "delta")[1])/1000, digits = 3)
ic50_ctrl_NAC <- signif(10^(ED(mcn, 50, interval = "delta")[1])/1000, digits = 3)
ic50_kd <- signif(10^(ED(mk, 50, interval = "delta")[1])/1000, digits = 3)
ic50_kd_NAC <- signif(10^(ED(mkn, 50, interval = "delta")[1])/1000, digits = 3)


p <- tbl %>%
  ggplot(aes(x = nM, y = `relative survival`, color = celltype, shape = celltype))+
  geom_point(size = 0.1)+
  geom_smooth(method = drm, method.args = list(fct = LN.4()), se = FALSE, linewidth = global_linewidth)+
  geom_errorbar(data = cp, aes(ymin=cp[,'Lower'] , ymax=cp[,'Upper'] ), linewidth = global_linewidth, width = 0.05)+
  geom_errorbar(data = cpn, aes(ymin=cpn[,'Lower'] , ymax=cpn[,'Upper'] ), linewidth = global_linewidth, width = 0.05)+
  geom_errorbar(data = kp, aes(ymin=kp[,'Lower'] , ymax=kp[,'Upper'] ), linewidth = global_linewidth, width = 0.05)+
  geom_errorbar(data = kpn, aes(ymin=kpn[,'Lower'] , ymax=kpn[,'Upper'] ), linewidth = global_linewidth, width = 0.05)+
  scale_x_continuous(breaks = 0:5, labels = c('1 nM', '10 nM', '100 nM', '1 μM', '10 μM', '100 μM'))+
  scale_y_continuous(breaks = seq(0, 1, length.out = 5))+
  theme_classic()+
  coord_cartesian(ylim = c(0, 1))+
  xlab('[vem]')+ ylab('relative survival(%)')+
  # emf cannot take alpha, so i need to change the color to pink and light blue and convert to blue and red in svg. (easier to do for line than points)
  scale_color_manual(values = c("#b9b6f2", "#aedef2", "#f1bfd6", "#f6cda9"))+
  # scale_color_manual(values = c('#0b00ff', "#04b3ff",  "#ff0074", "#ff7800"))+
  scale_shape_manual(values = c(16, 20, 15, 17))+
  geom_hline(yintercept=0.5,linetype = "dashed", colour= "grey", linewidth = global_linewidth)+
  guides(color = guide_legend('Treatment'), shape = guide_legend('Treatment', override.aes = list(size = 2, alpha = 1)))+
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


# ATP Production Rate of shTRMU -----
tbl <- readxl::read_excel(file.path(base, 'data', 'phenotypicAssays.xlsx'), sheet = 'seahorse', range = 'R9C1:R21C4', na = '0') %>% 
  dplyr::select(-Well) %>% 
  dplyr::rename_with(~ str_replace(.x, pattern = ' .*$', replacement = ''), ends_with('Rate')) 

pval <- tbl %>% 
  pivot_longer(ends_with('ATP'), names_to = 'Source', values_to = 'Rate') %>% 
  group_by(Source) %>%
  rstatix::t_test(Rate ~ Group) %>%
  # rstatix::anova_test(Rate ~ Group * Source) %>% 
  rstatix::add_significance("p")

p <- tbl %>% 
  pivot_longer(ends_with('ATP'), names_to = 'Source', values_to = 'Rate') %>% 
  ggpubr::ggbarplot(.,x = 'Group', y = 'Rate', fill = "Source", add = "mean_sd", position = position_stack(), width = .3, size = 0.1)+
  labs(x = '', y = 'ATP Production Rate (pmol/min)')+
  scale_fill_manual(name = '', values = c('mitoATP' = '#C3DBFD', 'glycoATP' = '#fffaca'))+
  theme(legend.position = c(0.8, 0.7),
        text = element_text(size = 5.3, family = 'arial'),
        legend.text = element_text(size = 5.3, family = 'arial'),
        legend.key.width = unit(0.18, "cm"),
        legend.key.height = unit(0.1, "cm"),
        axis.line = element_line(colour = 'black', linewidth = global_linewidth),
        axis.ticks = element_line(colour = "black", linewidth = global_linewidth))+
  # guides(fill = guide_legend(override.aes = list(size = 0.5)))+
  ggpubr::stat_pvalue_manual(pval[1,], y.position = 1200, step.increase = 0.1, label = "p.signif") # pval cannot be from anova test?? `Error in asserttat_group_columns_exists(data) : data should contain group1 and group2 columns`


ggsave(filename = file.path(base, 'figures/ATPprodRate.svg'), device = 'svg', plot = p, width = 2, height = 1.5, units = 'in')

# percentage of ATP production by mitochondria
pval <- tbl %>% 
  mutate(perc = mitoATP / (glycoATP + mitoATP)) %>% 
  rstatix::t_test(perc ~ Group) %>% 
  rstatix::add_significance("p")

p <- tbl %>% 
  mutate(perc = mitoATP / (glycoATP + mitoATP)) %>% 
  ggplot(aes(x = Group, y = perc)) +
  stat_summary(geom = 'bar', fun = 'mean' , color = 'black', fill = '#C3DBFD', na.rm = TRUE, width = .3, linewidth = 0.1)+
  stat_summary(geom = 'errorbar', fun.min = function(y) mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE), fun.max = function(y) mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE), color = 'black', width = 0.2 )+
  theme_classic()+
  # ggpubr::ggbarplot(., x = 'Group', y = 'perc', add = "mean_se", position = position_identity(), width = .3, size = 0.1)+
  labs(x = '', y = '% of mitoATP') +
  theme(legend.position = c(0.8, 0.7),
        text = element_text(size = 5.3, family = 'arial'),
        legend.text = element_text(size = 5.3, family = 'arial'),
        legend.key.width = unit(0.18, "cm"),
        legend.key.height = unit(0.1, "cm"),
        axis.line = element_line(colour = 'black', linewidth = global_linewidth),
        axis.ticks = element_line(colour = "black", linewidth = global_linewidth))+
  # guides(fill = guide_legend(override.aes = list(size = 0.5)))+
  ggpubr::stat_pvalue_manual(pval, y.position = 0.8, step.increase = 0.1, label = "p.signif")

ggsave(filename = file.path(base, 'figures/mitoPerc.svg'), device = 'svg', plot = p, width = 2, height = 1.5, units = 'in')


  