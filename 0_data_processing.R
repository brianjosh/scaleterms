
#**********************************************************
#* Code supporting "Lost in space: When scale terms blur 
#* #actual study size in plant community ecology"
#* Authors: (Hung, Perez-Navarro, Brian)
#* Code: Data cleaning and Management
#* ********************************************************


library(dplyr) 
library(tidyverse) 

## load raw data

# raw_df <- read.csv("./data/raw_data_clean.csv",  sep=",", header=T)

## add number of cases per category

n_cases <- raw_df%>%
  group_by(scale_name)%>%
  summarise(n_rows=n())

raw_df <- raw_df%>%
  left_join(n_cases, by="scale_name")

raw_df$scale_name%>%unique()

## add spatial category name

raw_df <- raw_df%>%
  mutate(scale_type=case_when(
    scale_name%in%c("quadrat", "plot", 
                    "transect", "site")~"method_scale",
    TRUE~"ecological_scale"
  ), 
  record = case_when(scale_size_m2=="NA" & scale_size_km=="NA"  ~ 'NA',
                     scale_size_m2>0 ~ "m2",
                     scale_size_km>0 ~ 'km',
                     scale_size_m2>0&scale_size_km>0~"both")
  )


write.table(raw_df, "./data/clean_scales_df.csv", 
            sep=",", dec=".", col.names = T,
            row.names=F)
size_df <- read.csv("./data/clean_scales_df.csv", sep=",", header=T)


raw_df <- read.csv("./data/clean_scales_df.csv",  sep=",", header=T)


## cleaning the database 

analyses_df <- raw_df%>%
  filter(n_rows>=7)


## creating new dataframe for ecological scale terms 
analyses_df_eco <- analyses_df%>%
  filter(!scale_name%in%c("quadrat", "plot", 
                           "transect", "site"))

## creating new dataframe for methodological scale terms 
analyses_df_meth <- analyses_df%>%
  filter(scale_name%in%c("quadrat", "plot", 
                          "transect", "site"))

sum_table <- raw_df %>%
  filter(n_rows>=7)%>% 
  group_by(scale_name)%>%
  summarise(median=median(scale_size_m2, na.rm=T),
            mean=mean(scale_size_m2, na.rm=T),
            sd=sd(scale_size_m2, na.rm=T), 
            count_na = sum(is.na(scale_size_m2)),
            n_rows=n(),
            quantile_10=quantile(scale_size_m2,probs=0.1, na.rm=T),
            quantile_25=quantile(scale_size_m2,probs=0.25, na.rm=T),
            quantile_75=quantile(scale_size_m2,probs=0.75, na.rm=T),
            quantile_90=quantile(scale_size_m2,probs=0.9, na.rm=T)
  )%>%
  mutate(perc_na=(count_na/n_rows)*100)


# repeat that code with scale_size_km_2

size_df <- size_df%>%
  mutate(scale_type=case_when(
    scale_name2%in%c("quadrat","plot",
                     "transect", "site")~ "method scale",
    scale_name2%in%c("regional", "landscape",
                     "area", "large-scale",
                     "local", "neighbourhood",
                     "patch")~"ecological scale"
  ))


