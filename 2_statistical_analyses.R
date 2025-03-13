

#**********************************************************
#* Code supporting "Lost in space: When scale terms blur 
#* #actual study size in plant community ecology"
#* Authors: (Hung, Perez-Navarro, Brian)
#* Code: Statistical analyses, ms figures and tables
#* ********************************************************

library(dplyr)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(effectsize)
library(effects)
library(lmerTest)
library(RColorBrewer)
library(MuMIn)
library(multcomp)
library(forcats)
library(grid)
library(gridExtra)
library(stringr)


# load and prepare data ------------------------------------------------

 analyses_df <- read.csv("./data/Data_S1b.csv", sep=",", header=T)
# analyses_df <- read.csv("./data/Data_S1.csv", sep=",", header=T)
# analyses_df <- read.csv("./data/Data_S1_div_habitat.csv", sep=",", header=T)

analyses_df <- analyses_df%>%
  filter(n_rows>=7)

# create separated dataframes for methodological and ecological scales 

analyses_df_meth <- analyses_df%>%
  filter(scale_type=="method_scale")

analyses_df_eco <- analyses_df%>%
  filter(scale_type=="ecological_scale")



# statistical analyses ----------------------------------------------------

# 1. models for scale size m2 vs scale name for eco dataframe -----------

# 1.0 general model

mod_lm_eco <- lm(log(scale_size_m2)~scale_name, data=analyses_df_eco)
sum_eco <- summary(mod_lm_eco)

em_cat_eco <- emmeans(mod_lm_eco,  "scale_name")
pairs(em_cat_eco)

letter_clas_eco <- multcomp::cld(em_cat_eco, alpha = 0.05,
                                 Letters = LETTERS)#consider to update letters to 0.05 before submitting
letter_clas_eco


mod_lm_eco_all <- lm(log(scale_size_m2)~scale_name + habitat
                      + region + type_study, data=analyses_df_eco)
summary(mod_lm_eco_all)
car::Anova(mod_lm_eco_all)

# 1.1 scale size m2 vs scale name * habitat

mod_lm_eco_hab <- lm(log(scale_size_m2)~habitat*scale_name,
                     data=analyses_df_eco)
summary(mod_lm_eco_hab)
car::Anova(mod_lm_eco_hab)

# emmeans

em_eco_hab <- emmeans(mod_lm_eco_hab, c("scale_name", "habitat"),
                      by="habitat",contr = "pairwise")

#it can be also written as
# em_eco_hab <- emmeans(mod_lm_eco_hab, ~scale_name|habitat,
#                       contr = "pairwise")

letter_clas_habi<- multcomp::cld(em_eco_hab, alpha = 0.05,
                                 Letters = LETTERS)

letter_clas_habi <-data.frame(letter_clas_habi)

#table S1
write.table(letter_clas_habi, "./output/scale_name_x_habitat_m2_eco.csv", 
            sep=",", dec=".", col.names = T,
            row.names=F)


# 1.2 scale size m2 vs scale name * region

mod_lm_eco_reg <- lm(log(scale_size_m2)~region*scale_name,
                     data=analyses_df_eco)
summary(mod_lm_eco_reg)
car::Anova(mod_lm_eco_reg)

# emmeans

em_eco_reg <- emmeans(mod_lm_eco_reg, c("scale_name", "region"), 
                       by="scale_name",contr = "pairwise")

letter_clas_reg<- multcomp::cld(em_eco_reg, alpha = 0.05,
                                 Letters = LETTERS)

letter_clas_reg <-data.frame(letter_clas_reg)

# table S2
write.table(letter_clas_reg, "./output/scale_name_x_region_m2_eco.csv", 
            sep=",", dec=".", col.names = T,
            row.names=F)

# 1.3 scale size m2 vs scale name * type_study

mod_lm_eco_stu <- lm(log(scale_size_m2)~type_study*scale_name,
                     data=analyses_df_eco)
summary(mod_lm_eco_stu)
car::Anova(mod_lm_eco_stu)

# emmeans

em_eco_stu <- emmeans(mod_lm_eco_stu, c("scale_name", "type_study"), 
                      by="scale_name",contr = "pairwise")

letter_clas_stu<- multcomp::cld(em_eco_stu, alpha = 0.05,
                                Letters = LETTERS)

letter_clas_stu <-data.frame(letter_clas_stu)

# table S3
write.table(letter_clas_stu, "./output/scale_name_x_type_study_m2_eco.csv", 
            sep=",", dec=".", col.names = T,
            row.names=F)


# 2. models for scale name vs scale size m2 for meth dataframe ------

# 2.0 general model
mod_lm_meth <- lm(log(scale_size_m2)~scale_name, data=analyses_df_meth)
summary(mod_lm_meth)

em_cat_meth <- emmeans(mod_lm_meth,  "scale_name")
pairs(em_cat_meth)

letter_clas_meth <- multcomp::cld(em_cat_meth, alpha = 0.05,
                                  Letters = LETTERS)
letter_clas_meth

mod_lm_meth_all <- lm(log(scale_size_m2)~scale_name + habitat
                      + region + type_study, data=analyses_df_meth)
summary(mod_lm_meth_all)
car::Anova(mod_lm_meth_all)


# 2.1 scale size m2 vs scale name * habitat

mod_lm_meth_hab <- lm(log(scale_size_m2)~habitat*scale_name,
                     data=analyses_df_meth)
summary(mod_lm_meth_hab)
car::Anova(mod_lm_meth_hab)

# emmeans

em_meth_hab <- emmeans(mod_lm_meth_hab, c("scale_name", "habitat"),
                      by="habitat",contr = "pairwise")

letter_clas_habi<- multcomp::cld(em_meth_hab, alpha = 0.05,
                                 Letters = LETTERS)

letter_clas_habi <-data.frame(letter_clas_habi)

# table S4
write.table(letter_clas_habi, "./output/scale_name_x_habitat_m2_meth.csv", 
            sep=",", dec=".", col.names = T,
            row.names=F)


# 2.2 scale size m2 vs scale name * region

mod_lm_meth_reg <- lm(log(scale_size_m2)~region*scale_name,
                     data=analyses_df_meth)
summary(mod_lm_meth_reg)
car::Anova(mod_lm_meth_reg)

# emmeans

em_meth_reg <- emmeans(mod_lm_meth_reg, c("scale_name", "region"), 
                      by="scale_name",contr = "pairwise")

letter_clas_reg<- multcomp::cld(em_meth_reg, alpha = 0.05,
                                Letters = LETTERS)

letter_clas_reg <-data.frame(letter_clas_reg)

# table S5
write.table(letter_clas_reg, "./output/scale_name_x_region_m2_meth.csv", 
            sep=",", dec=".", col.names = T,
            row.names=F)

# 2.3 scale size m2 vs scale name * type_study

mod_lm_meth_stu <- lm(log(scale_size_m2)~type_study*scale_name,
                     data=analyses_df_meth)
summary(mod_lm_meth_stu)
car::Anova(mod_lm_meth_stu)

# emmeans

em_meth_stu <- emmeans(mod_lm_meth_stu, c("scale_name", "type_study"), 
                      by="scale_name",contr = "pairwise")

letter_clas_stu<- multcomp::cld(em_meth_stu, alpha = 0.05,
                                Letters = LETTERS)

letter_clas_stu <-data.frame(letter_clas_stu)

# table S6
write.table(letter_clas_stu, "./output/scale_name_x_type_study_m2_meth.csv", 
            sep=",", dec=".", col.names = T,
            row.names=F)

# 3. models for scale name vs scale size km for eco dataframe ----

mod_lm_km_eco <- lm(log(scale_size_km)~scale_name, 
                    data=analyses_df_eco%>%
                      filter(scale_size_km>0))
summary(mod_lm_km_eco)

em_cat_km_eco <- emmeans(mod_lm_km_eco,  c("scale_name"))
pairs(em_cat_km_eco)
letter_clas_km_eco <- multcomp::cld(em_cat_km_eco, alpha = 0.05,
                                    Letters = LETTERS)
letter_clas_km_eco

analyses_df_eco%>%
  filter(!is.na(scale_size_km))%>%
  nrow()# only 33 data

# 4. models for scale name vs scale size km for meth dataframe -----

outlier_rows <- c(204,605)#these are in fact the only two rows of scale=site, rest are transect
mod_lm_km_meth <- lm(log(scale_size_km)~scale_name, 
                     data=analyses_df_meth
                     # data=analyses_df_meth%>%
                     #   slice(-outlier_rows)
                     )#remove outliers

summary(mod_lm_km_meth)

em_cat_km_meth <- emmeans(mod_lm_km_meth,  c("scale_name"))
pairs(em_cat_km_meth)
letter_clas_km_meth <- multcomp::cld(em_cat_km_meth, alpha = 0.05,
                                     Letters = LETTERS)
letter_clas_km_meth

analyses_df_meth%>%
  filter(!is.na(scale_size_km))%>%
  nrow()# only 22 data

# 5. shannon diversity models for eco dataframe --------------------

mod_lm_eco_shannon <- lm(shannon_diver_habitat~scale_name, data=analyses_df_eco)
summary(mod_lm_eco_shannon)

em_cat_eco_shannon <- emmeans(mod_lm_eco_shannon,  c("scale_name"))
pairs(em_cat_eco_shannon)
letter_clas_eco_shannon <- multcomp::cld(em_cat_eco_shannon, alpha = 0.05,
                                         Letters = LETTERS)
letter_clas_eco_shannon

# 6. shannon diversity models for meth dataframe --------------------

mod_lm_meth_shannon <- lm(shannon_diver_habitat~scale_name, data=analyses_df_meth)
summary(mod_lm_meth_shannon)

em_cat_meth_shannon <- emmeans(mod_lm_meth_shannon,  c("scale_name"))
pairs(em_cat_meth_shannon)
letter_clas_meth_shannon <- multcomp::cld(em_cat_meth_shannon, alpha = 0.05,
                                          Letters = LETTERS)
letter_clas_meth_shannon

# 7. simpson diversity models for eco dataframe ----------------------

mod_lm_eco_simpson <- lm(simpson_diver_habitat~scale_name, data=analyses_df_eco)
summary(mod_lm_eco_simpson)

em_cat_eco_simpson <- emmeans(mod_lm_eco_simpson,  c("scale_name"))
pairs(em_cat_eco_simpson)
letter_clas_eco_simpson <- multcomp::cld(em_cat_eco_simpson, alpha = 0.05,
                                         Letters = LETTERS)
letter_clas_eco_simpson

# 8. simpson diversity models for meth dataframe ---------------------

mod_lm_meth_simpson <- lm(simpson_diver_habitat~scale_name, data=analyses_df_meth)
summary(mod_lm_meth_simpson)

em_cat_meth_simpson <- emmeans(mod_lm_meth_simpson,  c("scale_name"))
pairs(em_cat_meth_simpson)
letter_clas_meth_simpson <- multcomp::cld(em_cat_meth_simpson, alpha = 0.05,
                                          Letters = LETTERS)
letter_clas_meth_simpson



# paper tables ------------------------------------------------------------

#Just grouping by scale terms
# sum_table <- analyses_df %>%
#   group_by(scale_name)%>%
#   summarise(median=median(scale_size_m2, na.rm=T),
#             mean=mean(scale_size_m2, na.rm=T),
#             sd=sd(scale_size_m2, na.rm=T), 
#             count_na = sum(is.na(record)),
#             n_rows=n(),
#             quantile_10=quantile(scale_size_m2,probs=0.1, na.rm=T),
#             quantile_25=quantile(scale_size_m2,probs=0.25, na.rm=T),
#             quantile_75=quantile(scale_size_m2,probs=0.75, na.rm=T),
#             quantile_90=quantile(scale_size_m2,probs=0.9, na.rm=T))%>%
#   mutate(perc_na=(count_na/n_rows)*100)
# sum_table

#Grouping by scale term within each habitat
sum_table_h <- analyses_df %>%
  group_by(habitat, scale_name)%>%
  summarise(median=median(scale_size_m2, na.rm=T),
            mean=mean(scale_size_m2, na.rm=T),
            sd=sd(scale_size_m2, na.rm=T), 
            count_na = sum(is.na(record)),
            n_rows=n(),
            quantile_10=quantile(scale_size_m2,probs=0.1, na.rm=T),
            quantile_25=quantile(scale_size_m2,probs=0.25, na.rm=T),
            quantile_75=quantile(scale_size_m2,probs=0.75, na.rm=T),
            quantile_90=quantile(scale_size_m2,probs=0.9, na.rm=T))%>%
  mutate(perc_na=(count_na/n_rows)*100)

sum_table_h

# paper figures -----------------------------------------------------------

#Fig. 1:


letter_clas_eco <- letter_clas_eco%>%
  mutate(group=.group,
         scale_type1="ecological scale")
letter_clas_meth <- letter_clas_meth%>%
  mutate(group=tolower(.group),
         scale_type1="method scale")

letter_clas_mod1_2 <- rbind(letter_clas_eco,
                            letter_clas_meth)

(gg_facet_m2 <- ggplot(data=analyses_df%>%
                        filter(!is.na(scale_size_m2))%>%
                        mutate(scale_name=as.factor(scale_name))%>%
                        mutate(scale_name = fct_reorder(scale_name,
                                                        scale_size_m2,
                                                        .fun="median"),
                               scale_type1=gsub("_", " ", scale_type)),
                      aes(x=scale_size_m2, y=scale_name, 
                          #fill=region2,
                          color=scale_name))+
  geom_boxplot(outlier.size = 0.15, lwd=0.3)+
  scale_color_viridis_d()+
  scale_x_continuous(trans='log10', limits=c(NA, 1e+17)) +
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab(expression(paste("Size (", m^2, ")"))) +
  geom_text(data = letter_clas_mod1_2, 
            aes(x = 1e+15, y = scale_name, label = group), 
            vjust=0.3, color = "black", size=2.5,
            family="serif")+
  facet_grid(rows=vars(scale_type1), scales = "free", 
             space = "free"
  )+
  theme_bw() +
  theme(
    legend.position="none",
    axis.title = element_text(size =8, color="black"),
    axis.text = element_text(size = 6, color="black"),
    strip.text= element_text(size=8),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
    text = element_text(family = "serif")
  )
)



ggsave(paste0("./figures/figure1.pdf"), 
       plot = gg_facet_m2, width = 8, height = 6, units = "cm")

ggsave(paste0("./figures/figure1.png"), 
       plot = gg_facet_m2, width = 8, height = 6, dpi = 600, units = "cm")

#Fig. 2:

letter_clas_eco_shannon <- letter_clas_eco_shannon%>%
  mutate(group=.group,
         scale_type1="ecological scale")
letter_clas_meth_shannon <- letter_clas_meth_shannon%>%
  mutate(group=tolower(.group),
         scale_type1="method scale")

letter_clas_mod5_6 <- rbind(letter_clas_eco_shannon,
                            letter_clas_meth_shannon)

(gg_facet_shannon <- ggplot(data=analyses_df%>%
                             filter(!is.na(scale_size_m2))%>%
                             filter(!is.na(scale_type))%>%
                             mutate(scale_name=as.factor(scale_name))%>%
                             mutate(scale_name = fct_reorder(scale_name,
                                                             shannon_diver_habitat,
                                                             .fun="median"),
                                    scale_type1=gsub("_", " ", scale_type)),
                           aes(x=shannon_diver_habitat, y=scale_name, 
                               #fill=region2,
                               color=scale_name))+
  geom_boxplot(outlier.size = 0.15, lwd=0.3)+
  scale_color_viridis_d()+
  xlim(0, 11)+
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab("Shannon Index")+
  geom_text(data = letter_clas_mod5_6, 
            aes(x = 9, y = scale_name, label = group), 
            vjust=0.3, color = "black", size=2.5,
            family="serif")+
  facet_grid(rows=vars(scale_type1), scales = "free", 
             space = "free" 
  )+
  ggtitle("a)")+
  theme_bw() +
  theme(
    legend.position="none",
    axis.title = element_text(size =9, color="black"),
    axis.text = element_text(size = 8, color="black"),
    strip.text=element_text(size=7, color="black"),
    plot.title = element_text(size = 10),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
    text=element_text(family="serif")
  )
)





letter_clas_eco_simpson <- letter_clas_eco_simpson%>%
  mutate(group=.group,
         scale_type1="ecological scale")
letter_clas_meth_simpson <- letter_clas_meth_simpson%>%
  mutate(group=tolower(.group),
         scale_type1="method scale")

letter_clas_mod7_8 <- rbind(letter_clas_eco_simpson,
                            letter_clas_meth_simpson)

(gg_facet_simpson <- ggplot(data=analyses_df%>%
                              filter(!is.na(scale_size_m2))%>%
                              filter(!is.na(scale_type))%>%
                              mutate(scale_name=as.factor(scale_name))%>%
                              mutate(scale_name = fct_reorder(scale_name,
                                                              simpson_diver_habitat,
                                                              .fun="median"),
                                     scale_type1=gsub("_", " ", scale_type)),
                            aes(x=simpson_diver_habitat, y=scale_name, 
                                color=scale_name))+
    geom_boxplot(outlier.size = 0.15, lwd=0.3)+
    scale_color_viridis_d()+
    xlim(0, 1.1)+
    ylab("Scale name") +
    xlab("Simpson Index")+
    geom_text(data = letter_clas_mod7_8, 
              aes(x = 1.05, y = scale_name, label = group), 
              vjust=0.3, color = "black", size=2.5,
              family="serif")+
    facet_grid(rows=vars(scale_type1), scales = "free", 
               space = "free" 
    )+
    ggtitle("b)")+
    theme_bw() +
    theme(
      legend.position="none",
      axis.title = element_text(size =9, color="black"),
      axis.text = element_text(size = 8, color="black"),
      strip.text=element_text(size=7, color="black"),
      plot.title = element_text(size = 10),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_line(linewidth = 0.25),
      text=element_text(family="serif")
    )
)



biplot <- gridExtra::grid.arrange(gg_facet_shannon,
                        gg_facet_simpson,
                        nrow=2,
                        heights=unit(c(60,60), "mm"), 
                        widths=unit(95, "mm"))


ggsave(paste0("./figures/figure2.pdf"), 
       plot = biplot, width = 9.6, height = 12, units = "cm")

ggsave(paste0("./figures/figure2.png"), 
       plot = biplot, width = 9.6, height = 12, dpi = 600, units = "cm")

#Fig. S1

#Based on the summary stats available in 'papers_by_year.xlsx'

years <- data.frame(N=c(37, 52, 66, 76, 76, 34), 
                    year=c('2001-2004', '2004-2008', '2009-2012', '2013-2016', '2017-2020', '2021-2023'), 
                    percent_no_info=c(0, 11.5, 4.5, 7.9, 11.8, 5.9),
                    percent_at_least_one=c(16.2, 23.1, 16.7, 25, 25, 14.7))
years$year <- as.factor(years$year)

yearplot <- ggplot(data=years) + 
  geom_line(aes(x=year, y=percent_no_info, group=1)) +
  geom_line(aes(x=year, y=percent_at_least_one, group=1), linetype="dashed") +
  theme_bw() + ylim(0, 50) +
  labs(x="Year", y="Percentage of studies with missing information") +
  geom_text(aes(x=year, y=43, group=1, label=N))
yearplot

#Fig. S2:

letter_clas_km_eco <- letter_clas_km_eco%>%
  mutate(group=.group,
         scale_type1="ecological scale")

letter_clas_km_meth <- letter_clas_km_meth%>%
  mutate(group=tolower(.group),
         scale_type1="method scale")

letter_clas_mod3_4 <- rbind(letter_clas_km_eco,
                            letter_clas_km_meth)

(gg_facet_km <- ggplot(data=analyses_df%>%
                        filter(!is.na(scale_size_km))%>%
                        mutate(scale_name=as.factor(scale_name))%>%
                        mutate(scale_name = fct_reorder(scale_name,
                                                        scale_size_km,
                                                        .fun="median"),
                               scale_type1=gsub("_", " ", scale_type)),
                      aes(x=scale_size_km, y=scale_name, 
                          #fill=region2,
                          color=scale_name))+
  geom_boxplot(outlier.size = 0.15, lwd=0.3)+
  scale_color_viridis_d()+
  scale_x_continuous(trans='log10', limits=c(NA, 1e+05))+
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab("Size (km)")+
  geom_text(data = letter_clas_mod3_4, 
            aes(x = 1e+04, y = scale_name, label = group),
            hjust=0.5, vjust=0.3, color = "black",
            size=2.5, family="serif")+
  facet_grid(rows=vars(scale_type1), scales = "free", 
             space = "free" 
  )+
  theme_bw() +
  theme(
    legend.position="none",
    axis.title = element_text(size =8, color="black"),
    axis.text = element_text(size = 6, color="black"),
    strip.text=element_text(size=8, color="black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
    text=element_text(family="serif")
  )
)
  
ggsave(paste0("./figures/figure_S1.pdf"), 
       plot = gg_facet_km, width = 10, height = 7, units = "cm")

ggsave(paste0("./figures/figure_S1.png"), 
       plot = gg_facet_km, width = 10, height = 7, dpi = 600, units = "cm")


#Fig. S3:


(gg_facet_habitat <- ggplot(data=analyses_df%>%
                             filter(!is.na(scale_size_m2))%>%
                             filter(!is.na(scale_type))%>%
                             mutate(scale_name=as.factor(scale_name))%>%
                             mutate(scale_name = fct_reorder(scale_name,
                                                             scale_size_m2,
                                                             .fun="median"),
                                    scale_type1=gsub("_", " ", scale_type))%>%
                             drop_na(habitat),
                           aes(x=scale_size_m2, y=scale_name, 
                               color=habitat))+
  geom_boxplot(outlier.size = 0.15, lwd=0.3,
               position = position_dodge2(preserve = "single"))+
  scale_color_viridis_d()+
  scale_x_continuous(trans='log10')+
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab(expression(paste("Size (", m^2, ")"))) +
  facet_grid(rows=vars(scale_type1), scales = "free", 
             space = "free" 
  )+
  theme_bw() +
  theme(
    legend.position="bottom",
    axis.title = element_text(size =12, color="black"),
    axis.text = element_text(size = 11, color="black"),
    strip.text=element_text(size=12, color="black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
    text=element_text(family="serif")
  )
)


ggsave(paste0("./figures/figure_S2.pdf"), 
       plot = gg_facet_habitat, width = 14, height = 17, units = "cm")

ggsave(paste0("./figures/figure_S2.png"), 
       plot = gg_facet_habitat, width = 14, height = 17, dpi = 600, units = "cm")




#Fig. S4:


(gg_facet_region <- ggplot(data=analyses_df%>%
                             filter(!is.na(scale_size_m2))%>%
                             filter(!is.na(scale_type))%>%
                             mutate(scale_name=as.factor(scale_name))%>%
                             mutate(scale_name = fct_reorder(scale_name,
                                                             scale_size_m2,
                                                             .fun="median"),
                                    scale_type1=gsub("_", " ", scale_type))%>%
                             drop_na(region),
                           aes(x=scale_size_m2, y=scale_name, 
                               color=region))+
    geom_boxplot(outlier.size = 0.15, lwd=0.3,
                 position = position_dodge2(preserve = "single"))+
    scale_color_viridis_d()+
    scale_x_continuous(trans='log10')+
    ylab("Scale name") +
    xlab(expression(paste("Size (", m^2, ")"))) +
    facet_grid(rows=vars(scale_type1), scales = "free", 
               space = "free" 
    )+
    theme_bw() +
    theme(
      legend.position="bottom",
      axis.title = element_text(size =15, color="black"),
      axis.text = element_text(size = 13, color="black"),
      strip.text=element_text(size=15, color="black"),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_line(linewidth = 0.25),
      text=element_text(family="serif")
    )
)


ggsave(paste0("./figures/figure_S3.pdf"), 
       plot = gg_facet_region, width = 16, height = 20, units = "cm")

ggsave(paste0("./figures/figure_S3.png"), 
       plot = gg_facet_region, width = 16, height = 20, dpi = 600, units = "cm")



#Fig. S5:

(gg_facet_study <- ggplot(data=analyses_df%>%
                           filter(!is.na(scale_size_m2))%>%
                           filter(!is.na(scale_type))%>%
                           mutate(scale_name=as.factor(scale_name))%>%
                           mutate(scale_name = fct_reorder(scale_name,
                                                           scale_size_m2,
                                                           .fun="median"),
                                  scale_type1=gsub("_", " ", scale_type))%>%
                           drop_na(type_study),
                         aes(x=scale_size_m2, y=scale_name, 
                             color=type_study))+
  geom_boxplot(outlier.size = 0.15, lwd=0.3,
               position = position_dodge2(preserve = "single"))+
  scale_color_viridis_d()+
  scale_x_continuous(trans='log10')+
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab(expression(paste("Size (", m^2, ")"))) +
  facet_grid(rows=vars(scale_type1), scales = "free", 
             space = "free" 
  )+
  theme_bw() +
  theme(
    legend.position="bottom",
    axis.title = element_text(size =12, color="black"),
    axis.text = element_text(size = 11, color="black"),
    strip.text=element_text(size=12, color="black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
    text=element_text(family="serif")
  )
)

ggsave(paste0("./figures/figure_S4.pdf"), 
       plot = gg_facet_study, width = 14, height = 17, units = "cm")

ggsave(paste0("./figures/figure_S4.png"), 
       plot = gg_facet_study, width = 14, height = 17, dpi = 600, units = "cm")



# Fig S6. Diversity by habitat type

analyses_df <- read.csv("./data/Data_S1b.csv", sep=",", header=T)


habitat_type <- c("grassland", "shrubland", "forest")


for (i in 1:length(habitat_type)){
  
  analyses_habitat <- analyses_df%>%
    filter(n_rows>=7)%>%
    filter(habitat%in%habitat_type[i])%>%
    drop_na(c(simpson_diver_habitat,
              shannon_diver_habitat))
  
  # create separated dataframes for methodological and ecological scales 
  
  analyses_df_meth <- analyses_habitat%>%
    filter(scale_type=="method_scale")
  
  analyses_df_eco <- analyses_habitat%>%
    filter(scale_type=="ecological_scale")
  
  
  # 5. shannon diversity models for eco dataframe 
  
  mod_lm_eco_shannon <- lm(shannon_diver_habitat~scale_name, data=analyses_df_eco)
  summary(mod_lm_eco_shannon)
  
  em_cat_eco_shannon <- emmeans(mod_lm_eco_shannon,  c("scale_name"))
  pairs(em_cat_eco_shannon)
  letter_clas_eco_shannon <- multcomp::cld(em_cat_eco_shannon, alpha = 0.05,
                                           Letters = LETTERS)
  letter_clas_eco_shannon
  
  # 6. shannon diversity models for meth dataframe 
  
  mod_lm_meth_shannon <- lm(shannon_diver_habitat~scale_name, data=analyses_df_meth)
  summary(mod_lm_meth_shannon)
  
  em_cat_meth_shannon <- emmeans(mod_lm_meth_shannon,  c("scale_name"))
  pairs(em_cat_meth_shannon)
  letter_clas_meth_shannon <- multcomp::cld(em_cat_meth_shannon, alpha = 0.05,
                                            Letters = LETTERS)
  letter_clas_meth_shannon
  
  # 7. simpson diversity models for eco dataframe 
  
  mod_lm_eco_simpson <- lm(simpson_diver_habitat~scale_name, data=analyses_df_eco)
  summary(mod_lm_eco_simpson)
  
  em_cat_eco_simpson <- emmeans(mod_lm_eco_simpson,  c("scale_name"))
  pairs(em_cat_eco_simpson)
  letter_clas_eco_simpson <- multcomp::cld(em_cat_eco_simpson, alpha = 0.05,
                                           Letters = LETTERS)
  letter_clas_eco_simpson
  
  # 8. simpson diversity models for meth dataframe 
  
  mod_lm_meth_simpson <- lm(simpson_diver_habitat~scale_name, data=analyses_df_meth)
  summary(mod_lm_meth_simpson)
  
  em_cat_meth_simpson <- emmeans(mod_lm_meth_simpson,  c("scale_name"))
  pairs(em_cat_meth_simpson)
  letter_clas_meth_simpson <- multcomp::cld(em_cat_meth_simpson, alpha = 0.05,
                                            Letters = LETTERS)
  letter_clas_meth_simpson
  
  
  
  #Fig. S7:
  
  letter_clas_eco_shannon <- letter_clas_eco_shannon%>%
    mutate(group=.group,
           scale_type1="ecological scale")
  letter_clas_meth_shannon <- letter_clas_meth_shannon%>%
    mutate(group=tolower(.group),
           scale_type1="method scale")
  
  letter_clas_mod5_6 <- rbind(letter_clas_eco_shannon,
                              letter_clas_meth_shannon)
  
  (gg_facet_shannon <- ggplot(data=analyses_habitat%>%
                                filter(!is.na(scale_size_m2))%>%
                                filter(!is.na(scale_type))%>%
                                mutate(scale_name=as.factor(scale_name))%>%
                                mutate(scale_name = fct_reorder(scale_name,
                                                                shannon_diver_habitat,
                                                                .fun="median"),
                                       scale_type1=gsub("_", " ", scale_type)),
                              aes(x=shannon_diver_habitat, y=scale_name, 
                                  #fill=region2,
                                  color=scale_name))+
      geom_boxplot(outlier.size = 0.15, lwd=0.3)+
      scale_color_viridis_d()+
      xlim(0, 11)+
      # geom_violin(width=1.1, size=0.2,
      #             color="transparent") +
      ylab("Scale name") +
      xlab("Shannon Index")+
      geom_text(data = letter_clas_mod5_6, 
                aes(x = 9, y = scale_name, label = group), 
                vjust=0.3, color = "black", size=2.5,
                family="serif")+
      facet_grid(rows=vars(scale_type1), scales = "free", 
                 space = "free" 
      )+
      labs(title = str_to_title(habitat_type[i]),
           subtitle = "a)")+
      theme_bw() +
      theme(
        legend.position="none",
        axis.title = element_text(size =9, color="black"),
        axis.text = element_text(size = 8, color="black"),
        strip.text=element_text(size=7, color="black"),
        plot.title = element_text(size = 10),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(linewidth = 0.25),
        text=element_text(family="serif")
      )
  )
  
  
  
  
  
  letter_clas_eco_simpson <- letter_clas_eco_simpson%>%
    mutate(group=.group,
           scale_type1="ecological scale")
  letter_clas_meth_simpson <- letter_clas_meth_simpson%>%
    mutate(group=tolower(.group),
           scale_type1="method scale")
  
  letter_clas_mod7_8 <- rbind(letter_clas_eco_simpson,
                              letter_clas_meth_simpson)
  
  (gg_facet_simpson <- ggplot(data=analyses_habitat%>%
                                filter(!is.na(scale_size_m2))%>%
                                filter(!is.na(scale_type))%>%
                                mutate(scale_name=as.factor(scale_name))%>%
                                mutate(scale_name = fct_reorder(scale_name,
                                                                simpson_diver_habitat,
                                                                .fun="median"),
                                       scale_type1=gsub("_", " ", scale_type)),
                              aes(x=simpson_diver_habitat, y=scale_name, 
                                  color=scale_name))+
      geom_boxplot(outlier.size = 0.15, lwd=0.3)+
      scale_color_viridis_d()+
      xlim(0, 1.1)+
      ylab("Scale name") +
      xlab("Simpson Index")+
      geom_text(data = letter_clas_mod7_8, 
                aes(x = 1.05, y = scale_name, label = group), 
                vjust=0.3, color = "black", size=2.5,
                family="serif")+
      facet_grid(rows=vars(scale_type1), scales = "free", 
                 space = "free" 
      )+
      labs(subtitle = "b)")+
      theme_bw() +
      theme(
        legend.position="none",
        axis.title = element_text(size =9, color="black"),
        axis.text = element_text(size = 8, color="black"),
        strip.text=element_text(size=7, color="black"),
        plot.title = element_text(size = 10),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(linewidth = 0.25),
        text=element_text(family="serif")
      )
  )
  
  
  
  biplot <- gridExtra::grid.arrange(gg_facet_shannon,
                                    gg_facet_simpson,
                                    nrow=2,
                                    heights=unit(c(70,65), "mm"), 
                                    widths=unit(95, "mm"))
  
  
  ggsave(paste0("./figures/figureS5", habitat_type[i], ".pdf"), 
         plot = biplot, width = 9.6, height = 14, units = "cm")
  
  ggsave(paste0("./figures/figureS5", habitat_type[i], ".png"), 
         plot = biplot, width = 9.6, height = 14, dpi = 600, units = "cm")
  
}


