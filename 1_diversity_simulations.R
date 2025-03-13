

#**********************************************************
#* Code supporting "Lost in space: When scale terms blur 
#* #actual study size in plant community ecology"
#* Authors: (Hung, Perez-Navarro, Brian)
#* Code: Diversity simulations
#* ********************************************************

#devtools::install_github("MoBiodiv/mobsim", build_vignettes = TRUE)
#library(mobsim)
library(sars)
library(dplyr) 
library(tidyverse) 
library(ggplot2) # for making plots
library(vegan)

mobsim_files <- list.files("./mobsim-master/R/",pattern="*.R",full.names = T)#in case mobsim doesn't run download it from github into the computer
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

sim_rich_tree <- function(A, c=NULL, z=NULL){
  
  c <- 0.8 #species pool per sq unit
  z <- 0.3
  error <- rnorm(7, 0, 3)
  #set.seed(3)
  sp_model <- c*(A)^z+sample(error*A^0.08, 1)
  n_sp_sim <- round(sp_model, 0)
  return(n_sp_sim)
  
} 

sim_rich_shrub <- function(A, c=NULL, z=NULL){
  
  c <- 2 #species pool per sq unit -here 1 m2-
  z <- 0.25
  error <- rnorm(7, 0, 3)
  #set.seed(3)
  sp_model <- c*(A)^z+sample(error*A^0.08, 1)
  n_sp_sim <- round(sp_model, 0)
  return(n_sp_sim)
  
} 

sim_rich_grass <- function(A, c=NULL, z=NULL){
  
  c <- 6 #species pool per sq unit
  z <- 0.22
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


sim_pred_tree <- function(A, c=NULL, z=NULL){
  
  c <- 0.8 #species pool per sq unit
  z <- 0.3
  #error <- rnorm(7, 0, 4)
  #set.seed(3)
  sp_model <- c*(A)^z
  n_sp_sim <- round(sp_model, 0)
  return(n_sp_sim)
  
} 

sim_pred_shrub <- function(A, c=NULL, z=NULL){
  
  c <- 2 #species pool per sq unit
  z <- 0.25
  #error <- rnorm(7, 0, 4)
  #set.seed(3)
  sp_model <- c*(A)^z
  n_sp_sim <- round(sp_model, 0)
  return(n_sp_sim)
  
} 

sim_pred_grass <- function(A, c=NULL, z=NULL){
  
  c <- 6 #species pool per sq unit
  z <- 0.22
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
    
    sample_vec <- base::sample(x = names(rel_abund_pool),
                               size = n_sim, replace = TRUE,
                               prob = rel_abund_pool)

    sample_vec <- factor(sample_vec, levels = names(rel_abund_pool))
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
  mutate(s=sim_rich(a),
         organism="unreal")%>%
  filter(s>=0)%>%
  as.data.frame()

random_sar_tree <- tibble(
  a=(seq(1,1000,length.out=2000)*rnorm(2000,1,0.3))^2  #rlnorm(200, log(30), log(2.8))
)%>%
  rowwise()%>%
  mutate(s=sim_rich_tree(a),
         organism="tree")%>%
  filter(s>=0)%>%
  as.data.frame()

random_sar_shrub <- tibble(
  a=(seq(1,1000,length.out=2000)*rnorm(2000,1,0.3))^2  #rlnorm(200, log(30), log(2.8))
)%>%
  rowwise()%>%
  mutate(s=sim_rich_shrub(a),
         organism="shrub")%>%
  filter(s>=0)%>%
  as.data.frame()

random_sar_grass <- tibble(
  a=(seq(1,1000,length.out=2000)*rnorm(2000,1,0.3))^2  #rlnorm(200, log(30), log(2.8))
)%>%
  rowwise()%>%
  mutate(s=sim_rich_grass(a),
         organism="grass")%>%
  filter(s>=0)%>%
  as.data.frame()

random_sar <- bind_rows(list(random_sar,
                             random_sar_tree,
                             random_sar_grass,
                             random_sar_shrub))

ggplot(data=random_sar,
       aes(x=a,y=s,
           color=organism))+
  geom_point()+
  xlab("Area")+
  ylab("Number of species")+
  theme_bw()


#* using datasets from SARS package
#* carefully check SARS package vignette
#* https://cran.r-project.org/web/packages/sars/vignettes/sars-r-package.html

#* names of datasets in the sars package
# d <- data(package = "sars")
# sars_datasets <- d$results[, "Item"]
# 

#* run logarithmic species-area relationship for example and some plant datasets

# display_sars_models()# to see all different sar equations 
# fit_random <- sar_power(data = random_sar_shrub%>%as.data.frame(),
#                        grid_start = "exhaustive", grid_n = 1000) #grid_start = "partial"
# plot(fit_random)
# fit_aegean2 <- sar_power(data = aegean2,
#                        grid_start = "exhaustive", grid_n = 1000)
# plot(fit_aegean2)
# fit_galap <- sar_power(data = galap,
#                       grid_start = "exhaustive", grid_n = 1000)# also sar_loga
# plot(fit_galap)

# 2. load size dataset -------------------------------------------------------


size_df <- read.csv("./data/clean_scales_df.csv", sep=",", header=T)

# 3. add species richness predictions ----------------------------------------

hist(size_df$scale_size_m2)
hist(log(size_df$scale_size_m2))
unique(size_df$habitat)

size_df <- size_df%>%
  rowwise()%>%
  mutate(sp_cz=sim_pred(A=scale_size_m2)-2,#simulate number of species
         sp_cz_grass=sim_pred_grass(A=scale_size_m2),
         sp_cz_shrub=sim_pred_shrub(A=scale_size_m2),
         sp_cz_tree= sim_pred_tree(A=scale_size_m2),
         sp_cz_habitat=case_when(
           habitat%in%c("multiple", "palustrine wetland", 
                        "desert and tundra")~ sp_cz,
           is.na(habitat)~sp_cz,
           habitat=="shrubland"~sp_cz_shrub, #also sp_cz as they are the same
           habitat=="grassland"~sp_cz_grass,
           habitat=="forest"~sp_cz_tree
         ))

size_df%>%
  filter(sp_cz<=0)

size_df%>%
  filter(scale_size_m2<=1)
#* for sites with less than 1 sq meter it might be normal to
#* obtain zero or negative richness values given the formula

hist(size_df$sp_cz)

ggplot(data=size_df,
       aes(x=scale_size_m2,
           y=sp_cz))+
  geom_point(color="black")+
  ylim(-3,100)+
  xlim(0, 1000000)+
  xlab("Area")+
  ylab("Number of species")+
  ggtitle("Own function")+
  theme_bw()

ggplot(data=size_df,
       aes(x=scale_size_m2,
           y=sp_cz_habitat,
           color=habitat))+
  geom_point()+
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
    TRUE~sp_cz),
    sp_cz_habitat=case_when(
      sp_cz_habitat<0~0,
      TRUE~sp_cz_habitat)
  )

# 4. add species abundance and diversity  prediction -------------------------

#* first we will create a total random number of individual as a linear function
#* of surface


size_df <- size_df%>%
  mutate(n_ind=5*scale_size_m2,
         n_ind_habitat=case_when(
           habitat%in%c("multiple", "palustrine wetland", 
                        "desert and tundra")~ 5*scale_size_m2,
           is.na(habitat)~5*scale_size_m2,
           habitat=="shrubland"~3*scale_size_m2,
           habitat=="forest"~0.2*scale_size_m2,
           habitat=="grassland"~100*scale_size_m2
         ))

#* 5 or any other number should be greater than c
#* I've chosen 5 as the number of individuals per sq meter. This is an overall
#* number to use in all type of communities thought this is not realistic
#* eg. we expect more individuals in grassland communities and less for trees
#* this simulation do not account neither for habitat fragmentation
#* 
#* #I've considered a conservative estimate for the remaining habitat categories
#* 3 individuals per sq meter for shrubland, 100 per sq meter as average in grassland
#* and 0.2 (2 per 10sq meters) as average tree number in forests.



## estimation based on general SAR and species abundance across habitats 

com_list <- list()

for(i in 1:nrow(size_df)){
  
  size_ex <- size_df[i, ]
  
  if(is.na(size_ex$scale_size_m2)){
    
    size_ex$simpson_diver <- NA_integer_
    size_ex$shannon_diver <- NA_integer_
    
  }else if(size_ex$sp_cz==0){
    
    size_ex$simpson_diver <- 0
    size_ex$shannon_diver <- 0
    
  }else if( !is.na(size_ex$scale_size_m2)& 
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


## estimation based on specific SAR and species abundance per habitat type 

com_list <- list()

for(i in 1:nrow(size_df)){
  
  size_ex <- size_df[i, ]
  
  if(is.na(size_ex$scale_size_m2)){
    
    size_ex$simpson_diver_habitat <- NA_integer_
    size_ex$shannon_diver_habitat <- NA_integer_
    
  } else if(size_ex$sp_cz_habitat == 0){
    
    size_ex$simpson_diver_habitat <- 0
    size_ex$shannon_diver_habitat <- 0
    
  } else {  # Removed unnecessary condition check
    
    abun <- tryCatch(
      {
        sim_sad_modif(
          s_pool = size_ex$sp_cz_habitat, 
          n_sim = size_ex$n_ind_habitat, 
          sad_type = "lnorm",
          sad_coef = list("meanlog" = 5, "sdlog" = 0.85)
        )
      },
      error = function(err) return(NULL)  # Returning NULL instead of NA for better checks
    )
    
    # Only process abun if it's valid
    if(!is.null(abun) && is.data.frame(abun)){
      
      # Diagnostic prints only if abun is valid
      head(abun)
      str(abun)
      
      if("abund_pool" %in% names(abun)){  # Check if column exists before plotting
        hist(abun$abund_pool)
      }
      
      # Calculate Simpson and Shannon diversity indexes
      size_ex$simpson_diver_habitat <- vegan::diversity(abun$abund_pool, "simpson")
      size_ex$shannon_diver_habitat <- vegan::diversity(abun$abund_pool, "shannon")
      
    } else {
      size_ex$simpson_diver_habitat <- NA_integer_
      size_ex$shannon_diver_habitat <- NA_integer_
    }
  }
  
  com_list[[i]] <- size_ex
  print(paste0(round(i / nrow(size_df) * 100, 2), "% completed"))  # Better formatted progress update
  
} #loop ends


size_diver <- bind_rows(com_list)%>%
  ungroup()

size_diver <- size_diver%>%
  rename(richness_habitat=sp_cz_habitat)




# 5. save table --------------------------------------------------------------


write.table(size_diver, "./data/Data_S1b.csv", 
            sep=",", dec=".", col.names = T,
            row.names=F)




