#Code supporting "Use of individual scale terms in plant ecology varies by 4.7 orders of magnitude"
#(Hung, Perez-Navarro, Brian)

library(tidyverse)
library(emmeans)
library(effectsize)
library(effects)
library(lmerTest)
library(RColorBrewer)
library(MuMIn)
library(multcomp)
library(forcats)

analyses_df <- read.csv("Data_S1.csv", header=T)
#The dataset only has terms that have >6 observations, and includes all variables required
#for the analysis

unique(analyses_df$scale_name)
n_distinct(analyses_df$title)

analyses_df <- analyses_df %>%
  mutate(record = case_when(scale_size_m2=="NA" & scale_size_km=="NA"  ~ 'NA',
                            scale_size_m2>0 | scale_size_km>0  ~ 'YES')) 

##creating new dataframe for ecological scale terms 
analyses_df_eco <- analyses_df%>%
  filter(scale_type=="ecological scale")

##creating new dataframe for methodological scale terms 
analyses_df_meth <- analyses_df%>%
  filter(scale_type=="method scale")

sum_table <- analyses_df %>%
  group_by(scale_name)%>%
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
sum_table

#Function to calculate order of magnitude differences: 

orderofmag <- function(small, big){
  x <- log10(big/small)
  print(x)
}

#Fig. 1:

gg_facet_m2 <- ggplot(data=analyses_df%>%
                     filter(!is.na(scale_size_m2))%>%
                     mutate(scale_name=as.factor(scale_name))%>%
                     mutate(scale_name = fct_reorder(scale_name,
                                                     scale_size_m2,
                                                     .fun="median")),
                   aes(x=scale_size_m2, y=scale_name, 
                       #fill=region2,
                       color=scale_name))+
  geom_boxplot()+
  scale_color_viridis_d()+
  scale_x_continuous(trans='log10', limits=c(NA, 1e+17)) +
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab(expression(paste("Size (", m^2, ")"))) +
  facet_grid(rows=vars(scale_type), scales = "free", 
             space = "free" 
  )+
  theme_bw() +
  theme(
    legend.position="none",
    axis.title = element_text(size =10, color="black"),
    axis.text = element_text(size = 8, color="black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
  )
gg_facet_m2

#Fig. 2:

gg_facet_shannon <- ggplot(data=analyses_df%>%
                             filter(!is.na(scale_size_m2))%>%
                             filter(!is.na(scale_type))%>%
                             mutate(scale_name=as.factor(scale_name))%>%
                             mutate(scale_name = fct_reorder(scale_name,
                                                             shannon_diver,
                                                             .fun="median")),
                           aes(x=shannon_diver, y=scale_name, 
                               #fill=region2,
                               color=scale_name))+
  geom_boxplot()+
  scale_color_viridis_d()+
  xlim(0, 11)+
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab("Shannon Index")+
  facet_grid(rows=vars(scale_type), scales = "free", 
             space = "free" 
  )+
  theme_bw() +
  theme(
    legend.position="none",
    axis.title = element_text(size =10, color="black"),
    axis.text = element_text(size = 8, color="black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
  )
gg_facet_shannon

#Fig. S1:

gg_facet_km <- ggplot(data=analyses_df%>%
                        filter(!is.na(scale_size_km))%>%
                        mutate(scale_name=as.factor(scale_name))%>%
                        mutate(scale_name = fct_reorder(scale_name,
                                                        scale_size_km,
                                                        .fun="median")),
                      aes(x=scale_size_km, y=scale_name, 
                          #fill=region2,
                          color=scale_name))+
  geom_boxplot()+
  scale_color_viridis_d()+
  scale_x_continuous(trans='log10', limits=c(NA, 1e+05))+
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab("Size (km)")+
  facet_grid(rows=vars(scale_type), scales = "free", 
             space = "free" 
  )+
  theme_bw() +
  theme(
    legend.position="none",
    axis.title = element_text(size =10, color="black"),
    axis.text = element_text(size = 8, color="black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
  )
gg_facet_km

#Fig. S2:

gg_facet_region <- ggplot(data=analyses_df%>%
                            filter(!is.na(scale_size_m2))%>%
                            filter(!is.na(scale_type))%>%
                            mutate(scale_name=as.factor(scale_name))%>%
                            mutate(scale_name = fct_reorder(scale_name,
                                                            scale_size_m2,
                                                            .fun="median"))%>%
                            drop_na(region),
                          aes(x=scale_size_m2, y=scale_name, 
                              color=region))+
  geom_boxplot()+
  scale_color_viridis_d()+
  scale_x_continuous(trans='log10')+
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab(expression(paste("Size (", m^2, ")"))) +
  facet_grid(rows=vars(scale_type), scales = "free", 
             space = "free" 
  )+
  theme_bw() +
  theme(
    legend.position="right",
    axis.title = element_text(size =10, color="black"),
    axis.text = element_text(size = 8, color="black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
  )
gg_facet_region

#Fig. S3:

gg_facet_habitat <- ggplot(data=analyses_df%>%
                             filter(!is.na(scale_size_m2))%>%
                             filter(!is.na(scale_type))%>%
                             mutate(scale_name=as.factor(scale_name))%>%
                             mutate(scale_name = fct_reorder(scale_name,
                                                             scale_size_m2,
                                                             .fun="median"))%>%
                             drop_na(habitat),
                           aes(x=scale_size_m2, y=scale_name, 
                               color=habitat))+
  geom_boxplot()+
  scale_color_viridis_d()+
  scale_x_continuous(trans='log10')+
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab(expression(paste("Size (", m^2, ")"))) +
  facet_grid(rows=vars(scale_type), scales = "free", 
             space = "free" 
  )+
  theme_bw() +
  theme(
    legend.position="right",
    axis.title = element_text(size =10, color="black"),
    axis.text = element_text(size = 8, color="black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
  )
gg_facet_habitat

#Fig. S4:

gg_facet_study <- ggplot(data=analyses_df%>%
                           filter(!is.na(scale_size_m2))%>%
                           filter(!is.na(scale_type))%>%
                           mutate(scale_name=as.factor(scale_name))%>%
                           mutate(scale_name = fct_reorder(scale_name,
                                                           scale_size_m2,
                                                           .fun="median"))%>%
                           drop_na(type_study),
                         aes(x=scale_size_m2, y=scale_name, 
                             color=type_study))+
  geom_boxplot()+
  scale_color_viridis_d()+
  scale_x_continuous(trans='log10')+
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab(expression(paste("Size (", m^2, ")"))) +
  facet_grid(rows=vars(scale_type), scales = "free", 
             space = "free" 
  )+
  theme_bw() +
  theme(
    legend.position="right",
    axis.title = element_text(size =10, color="black"),
    axis.text = element_text(size = 8, color="black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
  )
gg_facet_study

#Fig. S5:

gg_facet_simpson <- ggplot(data=analyses_df%>%
                             filter(!is.na(scale_size_m2))%>%
                             filter(!is.na(scale_type))%>%
                             mutate(scale_name=as.factor(scale_name))%>%
                             mutate(scale_name = fct_reorder(scale_name,
                                                              simpson_diver,
                                                              .fun="median")),
                           aes(x=simpson_diver, y=scale_name, 
                               #fill=region2,
                               color=scale_name))+
  geom_boxplot()+
  scale_color_viridis_d()+
  xlim(0, 1.25)+
  # geom_violin(width=1.1, size=0.2,
  #             color="transparent") +
  ylab("Scale name") +
  xlab("Simpson's Index of Diversity")+
  facet_grid(rows=vars(scale_type), scales = "free", 
             space = "free" 
  )+
  theme_bw() +
  theme(
    legend.position="none",
    axis.title = element_text(size =10, color="black"),
    axis.text = element_text(size = 8, color="black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.25),
  )
gg_facet_simpson

#Modelling differences between scale terms

# models for scale name vs scale size m2 for eco dataframe

mod_lm_eco <- lm(log(scale_size_m2)~scale_name, data=analyses_df_eco)
sum_eco <- summary(mod_lm_eco)

em_cat_eco <- emmeans(mod_lm_eco,  c("scale_name"))
pairs(em_cat_eco)

letter_clas_eco <- multcomp::cld(em_cat_eco, alpha = 0.10,
                                 Letters = LETTERS)
letter_clas_eco

mod_lm_eco_all <- lm(log(scale_size_m2)~scale_name + habitat
                      + region + type_study, data=analyses_df_eco)
summary(mod_lm_eco_all)
Anova(mod_lm_eco_all)

# models for scale name vs scale size m2 for meth dataframe 

mod_lm_meth <- lm(log(scale_size_m2)~scale_name, data=analyses_df_meth)
summary(mod_lm_meth)

em_cat_meth <- emmeans(mod_lm_meth,  c("scale_name"))
pairs(em_cat_meth)

letter_clas_meth <- multcomp::cld(em_cat_meth, alpha = 0.10,
                                  Letters = LETTERS)
letter_clas_meth

mod_lm_meth_all <- lm(log(scale_size_m2)~scale_name + habitat
                      + region + type_study, data=analyses_df_meth)
summary(mod_lm_meth_all)
Anova(mod_lm_meth_all)

#emmeans for km and scale name for eco dataframe 
mod_lm_km_eco <- lm(scale_size_km~scale_name, data=analyses_df_eco)
summary(mod_lm_km_eco)

em_cat_km_eco <- emmeans(mod_lm_km_eco,  c("scale_name"))
pairs(em_cat_km_eco)
letter_clas_km_eco <- multcomp::cld(em_cat_km_eco, alpha = 0.10,
                                    Letters = LETTERS)
letter_clas_km_eco

#emmeans for km and scale name for meth dataframe 
mod_lm_km_meth <- lm(scale_size_km~scale_name, data=analyses_df_meth)
summary(mod_lm_km_meth)

em_cat_km_meth <- emmeans(mod_lm_km_meth,  c("scale_name"))
pairs(em_cat_km_meth)
letter_clas_km_meth <- multcomp::cld(em_cat_km_meth, alpha = 0.10,
                                     Letters = LETTERS)
letter_clas_km_meth

## shannon diversity models for eco dataframe
mod_lm_eco_shannon <- lm(shannon_diver~scale_name, data=analyses_df_eco)
summary(mod_lm_eco_shannon)

em_cat_eco_shannon <- emmeans(mod_lm_eco_shannon,  c("scale_name"))
pairs(em_cat_eco_shannon)
letter_clas_eco_shannon <- multcomp::cld(em_cat_eco_shannon, alpha = 0.10,
                                 Letters = LETTERS)
letter_clas_eco_shannon

## shannon diversity models for meth dataframe
mod_lm_meth_shannon <- lm(shannon_diver~scale_name, data=analyses_df_meth)
summary(mod_lm_meth_shannon)

em_cat_meth_shannon <- emmeans(mod_lm_meth_shannon,  c("scale_name"))
pairs(em_cat_meth_shannon)
letter_clas_meth_shannon <- multcomp::cld(em_cat_meth_shannon, alpha = 0.10,
                                         Letters = LETTERS)
letter_clas_meth_shannon

## simpson diversity models for eco dataframe
mod_lm_eco_simpson <- lm(simpson_diver~scale_name, data=analyses_df_eco)
summary(mod_lm_eco_simpson)

em_cat_eco_simpson <- emmeans(mod_lm_eco_simpson,  c("scale_name"))
pairs(em_cat_eco_simpson)
letter_clas_eco_simpson <- multcomp::cld(em_cat_eco_simpson, alpha = 0.10,
                                         Letters = LETTERS)
letter_clas_eco_simpson

## simpson diversity models for meth dataframe
mod_lm_meth_simpson <- lm(simpson_diver~scale_name, data=analyses_df_meth)
summary(mod_lm_meth_simpson)

em_cat_meth_simpson <- emmeans(mod_lm_meth_simpson,  c("scale_name"))
pairs(em_cat_meth_simpson)
letter_clas_meth_simpson <- multcomp::cld(em_cat_meth_simpson, alpha = 0.10,
                                          Letters = LETTERS)
letter_clas_meth_simpson

###################################

#Code for simulating and sampling the virtual community
#(The full data set loaded in line 14 has the output from this procedure already included as
#columns; this code just shows how it was done)
#Code prepared by Maria Angeles Perez-Navarro: maria_angeles.perez-navarro@kcl.ac.uk

devtools::install_github("MoBiodiv/mobsim", build_vignettes = TRUE)
library(mobsim)
library(sars)
library(vegan)

mobsim_files <- list.files("./mobsim-master/R/",pattern="*.R",full.names = T)
sapply(mobsim_files,source,.GlobalEnv)

# 0. Preparing R functions -------------------------------------------------

#* preparing functions for creating random richness (including random variability)
#* -sim_rich- and to produce specific values of richness (not including variability)
#* -sim_pred-

#* preparing function for estimating species-abundance relationship

sim_rich <- function(A, c=NULL, z=NULL){
  
  c <- 3 #species pool per sq unit
  z <- 0.25
  error <- rnorm(7, 0, 3)
  #set.seed(3)
  sp_model <- c*(A)^z+sample(error*A^0.08, 1)
  n_sp_sim <- round(sp_model, 0)
  return(n_sp_sim)
  
} 

sim_pred <- function(A, c=NULL, z=NULL){
  
  c <- 3 #species pool per sq unit
  z <- 0.25
  #error <- rnorm(7, 0, 4)
  #set.seed(3)
  sp_model <- c*(A)^z
  n_sp_sim <- round(sp_model, 0)
  return(n_sp_sim)
  
} 

sim_sad_modif <- function(s_pool, n_sim,
                          sad_type = c("lnorm", "bs", "gamma", "geom", "ls",
                                       "mzsm","nbinom", "pareto", "poilog", "power",
                                       "powbend", "weibull"),
                          sad_coef = list("cv_abund" = 1),
                          fix_s_sim = FALSE,
                          drop_zeros = TRUE)
{
  sad_type <- match.arg(sad_type)
  
  if (!is.numeric(n_sim) || n_sim <= 0)
    stop("n_sim has to be a positive integer number")
  
  n_sim <- round(n_sim, digits = 0)
  
  if (class(sad_coef) != "list" | is.null(names(sad_coef))) stop("coef must be a named list!")
  
  # Handles parameters that give the community size
  if (sad_type %in% c("bs", "ls", "mzsm")) {
    S <- switch(sad_type,
                bs = sad_coef$S,
                ls = sad_coef$alpha * log ( 1 + sad_coef$N / sad_coef$alpha ),
                mzsm = sum(sad_coef$theta / (1:sad_coef$J) *
                             (1 - (1:sad_coef$J)/sad_coef$J)^(sad_coef$theta - 1))
    )
    S <- round(S)
    if (!is.null(s_pool)){
      warning(paste("For the selected SAD model the value of s_pool is ignored.
  s_pool calculated from the SAD model coefficients is", S, "species."))
    }
    s_pool <- S
  } else {
    if (is.null(s_pool) || is.na(s_pool) || !is.numeric(s_pool) || s_pool <= 0)
      stop("The argument s_pool is mandatory for the selected sad and has to be a positive integer number.")
    s_pool <- round(s_pool, digits = 0)
  }
  
  if (s_pool > 1){
    
    #alternative parameterization for lnorm and poilog
    if ((sad_type == "lnorm" || sad_type == "poilog") &&
        names(sad_coef)[1] == "cv_abund"){
      
      mean_abund <- n_sim/s_pool
      sd_abund <-  mean_abund * sad_coef$cv_abund
      sigma1 <- sqrt(log(sd_abund^2/mean_abund^2 + 1))
      mu1 <- log(mean_abund) - sigma1^2/2
      
      # mean1 <- exp(mu1 + sigma1^2/2)
      # sd1 <- exp(mu1 + sigma1^2/2) * sqrt(exp(sigma1^2) - 1)
      # cv1 <- sd1/mean1
      
      if (sad_type == "lnorm" & !is.null(sad_coef$cv_abund) ){
        sad_coef <- list("meanlog" = mu1, "sdlog" = sigma1)
      }
      
      if (sad_type == "lnorm" & !is.null(sad_coef$meanlog) & !is.null(sad_coef$meanlog)){
        mu1 <- sad_coef$meanlog
        sigma1 <- sad_coef$meanlog
        sad_coef <- list("meanlog" = mu1, "sdlog" = sigma1)
      }
      
      if (sad_type == "poilog"){
        sad_coef <- list("mu" = mu1, "sig" = sigma1)
        
      }
      
    }
    
    # Generates the "community"
    if (sad_type %in% c("gamma","geom","lnorm","nbinom","weibull")){
      sadr <- utils::getFromNamespace(paste("r", sad_type, sep=""), ns = "stats")
    } else {
      sadr <- utils::getFromNamespace(paste("r", sad_type, sep=""), ns = "sads")
    }
    
    abund_pool <- do.call(sadr, c(list(n = s_pool), sad_coef))
    
    # abund_pool <- abund_pool[abund_pool > 0]
    rel_abund_pool <- abund_pool/sum(abund_pool)
    rel_abund_pool <- sort(rel_abund_pool, decreasing = T)
    names(rel_abund_pool) <- paste0("species",seq_along(rel_abund_pool))
    
    ab_df <- as.data.frame(abund_pool)
    ab_df$species <- paste0("species",seq_along(abund_pool))
    ab_df$abund_rel <- ab_df$abund_pool/sum(abund_pool)
    ab_df <- ab_df[, c("species", "abund_pool", "abund_rel")]
    
    # sample_vec <- base::sample(x = names(rel_abund_pool),
    #                            size = n_sim, replace = TRUE,
    #                            prob = rel_abund_pool)
    # 
    # sample_vec <- factor(sample_vec, levels = names(rel_abund_pool))
    # 
    abund_local <- table(sample_vec)
    # s_local <- sum(abund_local > 0)
    
    s_local <- ab_df%>%
      filter(abund_pool>0)%>%
      nrow()
    
    if (fix_s_sim == TRUE & s_local < s_pool){
      
      s_diff <- s_pool - s_local
      abund_local[abund_local == 0] <- 1
      ab_df <- ab_df%>%
        mutate(abund_pool1=case_when(
          abund_pool==0~1,
          abund_pool>0~abund_pool
        ))
      
      n <- sum(ab_df$abund_pool1)
      ab_df <- ab_df%>%
        select(-abund_pool)%>%
        rename(abund_pool=abund_pool1)
      
      #randomly remove individuals until target level is reached
      while (n > n_sim){
        rel_abund <- abund_local/sum(abund_local)
        # draw proportional to relative abundance
        irand <- sample(1:s_pool, size = 1, prob = rel_abund)
        if (abund_local[irand] > 1) abund_local[irand] <- abund_local[irand] - 1
        n <- sum(abund_local)
        
        ab_df <- as.data.frame(abund_local)
        ab_df$species <- paste0("species",seq_along(abund_local))
        ab_df$abund_rel <- ab_df$abund_local/sum(abund_local)
        ab_df <- ab_df[, c("species", "abund_pool", "abund_rel")]
    
       
      }
    }
  } else { # end if(s_pool > 1)
    abund_local <- n_sim
  }
  
  #names(abund2) <- paste("species", 1:length(abund2), sep = "")
  if (drop_zeros == T){
    ab_df <- ab_df%>%
      filter(abund_pool>0)
  }
  
  #class(abund_local) <- c("sad","integer")
  return(ab_df)
}
# 1. create SAR model based on random data or existing package dataset --------

#* using sim_rich function to simulate random richness

random_sar <- tibble(
  a=(seq(1,1000,length.out=2000)*rnorm(2000,1,0.3))^2  #rlnorm(200, log(30), log(2.8))
)%>%
  rowwise()%>%
  mutate(s=sim_rich(a))%>%
  filter(s>=0)%>%
  as.data.frame()

ggplot(data=random_sar,
       aes(x=a,y=s))+
  geom_point(color="black")+
  xlab("Area")+
  ylab("Number of species")+
  theme_bw()

#* using datasets from SARS package
#* carefully check SARS package vignette
#* https://cran.r-project.org/web/packages/sars/vignettes/sars-r-package.html

#* names of datasets in the sars package
d <- data(package = "sars")
sars_datasets <- d$results[, "Item"]

#* run logarithmic species-area relationship for example and some plant datasets

# fit_random <- sar_loga(data = random_sar%>%as.data.frame(), 
#                        grid_start = "exhaustive", grid_n = 1000) #grid_start = "partial"
# fit_aegean <- sar_loga(data = aegean2,
#                        grid_start = "exhaustive", grid_n = 1000) 
# fit_galap <- sar_loga(data = galap, 
#                       grid_start = "exhaustive", grid_n = 1000) 

# 2. load size dataset -------------------------------------------------------

size_df <- read.csv("./data/clean_data2_modifiedhabitats.csv", sep=",", header=T)

# 3. add species richness predictions ----------------------------------------

hist(size_df$scale_size_m2_2)
hist(log(size_df$scale_size_m2_2))

size_df <- size_df%>%
  rowwise()%>%
  mutate(sp_cz=sim_pred(A=scale_size_m2_2)-2)#simulate number of species

size_df%>%
  filter(sp_cz<=0)

size_df%>%
  filter(scale_size_m2_2<=1)
#* for sites with less than 1 sq meter it might be normal to
#* obtain zero or negative richness values given the formula

hist(size_df$sp_cz)

ggplot(data=size_df,
       aes(x=scale_size_m2_2,
           y=sp_cz))+
  geom_point(color="black")+
  ylim(-3,100)+
  xlim(0, 1000000)+
  xlab("Area")+
  ylab("Number of species")+
  ggtitle("Own function")+
  theme_bw()

# replace negative values by 0
size_df <- size_df%>%
  mutate(sp_cz=case_when(
    sp_cz<0~0,
    TRUE~sp_cz
  ))

# 4. add species abundance and diversity  prediction -------------------------

#* first we will create a total random number of individual as a linear function
#* of surface

size_df <- size_df%>%
  mutate(n_ind=5*scale_size_m2_2)

#* 5 or any other number should be greater than c
#* I've chosen 5 as the number of individuals per sq meter. This is an overall
#* number to use in all type of communities thought this is not realistic
#* eg. we expect more individuals in grassland communities and less for trees
#* this simulation do not account neither for habitat fragmentation

com_list <- list()

for(i in 1:nrow(size_df)){
  
  size_ex <- size_df[i, ]
  
  if(is.na(size_ex$scale_size_m2_2)){
    
    size_ex$simpson_diver <- NA_integer_
    size_ex$shannon_diver <- NA_integer_
    
  }else if(size_ex$sp_cz==0){
    
    size_ex$simpson_diver <- 0
    size_ex$shannon_diver <- 0
    
  }else if( !is.na(size_ex$scale_size_m2_2)& 
            size_ex$sp_cz!=0){
    
    abun <- tryCatch(
      
      sim_sad_modif(s_pool = size_ex$sp_cz, 
                    n_sim = size_ex$n_ind , 
                    sad_type = "lnorm",
                    sad_coef = list("meanlog" = 5, "sdlog" = 0.85)),
      
      error=function(err) NA
      
    )
    
    head(abun)
    str(abun)
    hist(abun$abund_pool)
    
    #calculate simpson and shannon diversity
    
    if(is.data.frame(abun)==F){#unique(is.na(abun))
      
      size_ex$diver <- NA_integer_
      
    }else if(is.data.frame(abun)==T){#unique(!is.na(abun))
      
      size_ex$simpson_diver <- vegan::diversity(abun$abund_pool, "simpson")
      size_ex$shannon_diver <- vegan::diversity(abun$abund_pool, "shannon")
      # check other diversity indexes as well
      
    }
    
  }#conditions ends
  
  com_list[[i]] <- size_ex
  print(i/nrow(size_df)*100)# to see percentage of rows already processed
  
}# loop ends

size_diver <- bind_rows(com_list)%>%
  ungroup()

size_diver <- size_diver%>%
  rename(richness=sp_cz)

# 5. save table, plot and analyze scale size name vs diversity ---------------

write.table(size_diver, "./data/clean_diversity_df_modifiedhabitats2.csv", 
            sep=",", dec=".", col.names = T,
            row.names=F)

  
  
  
  
  
  