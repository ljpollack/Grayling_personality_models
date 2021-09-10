#data frame
assays <- read.csv("Grayling_behavioraldata_clean.csv")

##load libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(brms)
library(performance)
library(tidybayes)

#### my icc tibble function ####
#function reports icc for both hurdle and negbinomial portions of movel
my_icc_tibble <- function(model, total_re.form, lesser_re.form, new_data = NULL, hurdle_sep = FALSE){
  
  is_multi <- insight:::is_multivariate(model)
  
  PPD <- posterior_predict(model, re.form = total_re.form, newdata = new_data)
  PPD_0 <- posterior_predict(model, re.form = lesser_re.form, newdata = new_data)
  
  if(is_multi){
    PPD_list <- array_branch(PPD, 3)
    PPD_0_list <- array_branch(PPD_0, 3)
    
    PPD_list <- map(.x = PPD_list, function(x)  apply(x, 1, FUN = stats::var, na.rm = T))
    PPD_0_list <- map(.x = PPD_0_list, function(x)  apply(x, 1, FUN = stats::var, na.rm = T))
    
    ICC_list <- map2(.x = PPD_0_list, .y = PPD_list, .f = function(x,y) 1 - (x/y))
    
    map_dfr(ICC_list, .f = function(x) x %>% median_hdi() %>% 
              rename(median_icc = y, ci_2.5 = ymin, ci_97.5 = ymax) %>% 
              select(median_icc, ci_2.5, ci_97.5), .id = "outcome")
  } else{
    
    if(hurdle_sep){
      
      # h is hurdle, d is other distribution
      PPD_h <- PPD
      PPD_h[PPD_h > 0] <- 1
      
      PPD_0_h <- PPD_0
      PPD_0_h[PPD_0_h > 0] <- 1
      
      PPD_d <- PPD
      PPD_d[PPD_d == 0] <- NA
      
      PPD_0_d <- PPD_0
      PPD_0_d[PPD_0_d == 0] <- NA
      
      vars_h <- apply(PPD_h, MARGIN = 1, FUN = var)
      vars_0_h <- apply(PPD_0_h, MARGIN = 1, FUN = var)
      
      vars_d <- apply(PPD_d, MARGIN = 1, FUN = var, na.rm = T)
      vars_0_d <- apply(PPD_0_d, MARGIN = 1, FUN = var, na.rm = T)
      
      icc_draws_h <- tibble(icc_draws = 1 - vars_0_h/vars_h)
      icc_draws_d <- tibble(icc_draws = 1 - vars_0_d/vars_d)
      
      model_name <- deparse(substitute(model))
      
      icc_draws_h <- icc_draws_h %>% 
        median_hdi() %>% 
        rename(median_icc = icc_draws, ci_2.5 = .lower, ci_97.5 = .upper) %>% 
        mutate(model_name = model_name, distribution = "hurdle") %>% 
        select(model_name, distribution, median_icc, ci_2.5, ci_97.5)
      
      icc_draws_d <- icc_draws_d %>% 
        median_hdi() %>% 
        rename(median_icc = icc_draws, ci_2.5 = .lower, ci_97.5 = .upper) %>% 
        mutate(model_name = model_name, distribution = "post-hurdle") %>% 
        select(model_name, distribution, median_icc, ci_2.5, ci_97.5)
      
      bind_rows(icc_draws_h, icc_draws_d)
      
    } else {
      
      vars <- apply(PPD, MARGIN = 1, FUN = var)
      vars_0 <- apply(PPD_0, MARGIN = 1, FUN = var)
      
      icc_draws <- tibble(icc_draws = 1 - (vars_0/vars))
      model_name <- deparse(substitute(model))
      icc_draws %>% 
        median_hdi() %>% 
        rename(median_icc = icc_draws, ci_2.5 = .lower, ci_97.5 = .upper) %>% 
        mutate(model_name = model_name) %>% 
        select(model_name, median_icc, ci_2.5, ci_97.5)
    }
  }
}

#### generating z-scores for important predictor variables ####
assays$ST_length <- scale(assays$ST_length, center=TRUE, scale = TRUE)
assays$K_ST <- scale(assays$K_ST, center=TRUE, scale = TRUE)

assays$X13C_scaled <- scale(assays$X13C, center=TRUE, scale = TRUE)
assays$X15N_scaled <- scale(assays$X15N, center=TRUE, scale = TRUE)

####latency to cross barrier analysis####

#round latency to cross seconds to integer
assays$Latcross_s_r <- as.numeric(round(assays$Latcross_s))

#Run basic model with hurdle
hurd_latcross_1 <- brm(data = assays, family = hurdle_negbinomial(),
                     bf(Latcross_hurdle  ~ 1 + Type + Trial + (1| ID) + (1 | Pool),
                        hu ~ 1 + Type + Trial + (1| ID) + (1 | Pool)),
                     prior = c(set_prior("normal(0,3)", class = "b"),
                               set_prior("cauchy(0,2)", class = "sd")),
                     iter = 5000, warmup = 1000, 
                     chains = 3, cores = future::availableCores(),
                     control = list(adapt_delta = 0.99999999999999, max_treedepth = 27))
  
  ##inspect model
  launch_shinystan(hurd_latcross_1) 
  summary(hurd_latcross_1) ## nothing strongly predicts latency to cross
  icc(hurd_latcross_1, re.form= ~(1| ID), typical = "mean", prob = 0.95, ppd = TRUE) ##time to cross is not repeatable
  ##use function for ICC of hurdle portion of model
  my_icc_tibble(hurd_latcross_1, ~(1|ID) + (1|Pool), ~(1|Pool), hurdle_sep = T) ##neither is repeatable

#Try full model with length and sex
hurd_latcross_2 <- brm(data = assays, family = hurdle_negbinomial(),
                      bf(Latcross_hurdle  ~ 1 + Type + Trial + ST_length + Sex_sim + (1| ID) + (1 | Pool),
                         hu ~ 1 + Type + Trial + ST_length + Sex_sim + (1| ID) + (1 | Pool)),
                      prior = c(set_prior("normal(0,3)", class = "b"),
                                set_prior("cauchy(0,2)", class = "sd")),
                      iter = 5000, warmup = 1000, 
                      chains = 3, cores = future::availableCores(),
                      control = list(adapt_delta = 0.99999999999999, max_treedepth = 27))
  
launch_shinystan(hurd_latcross_2) 

#remove sex
hurd_latcross_3 <- brm(data = assays, family = hurdle_negbinomial(),
                       bf(Latcross_hurdle  ~ 1 + Type + Trial + ST_length + (1| ID) + (1 | Pool),
                          hu ~ 1 + Type + Trial + ST_length + (1| ID) + (1 | Pool)),
                       prior = c(set_prior("normal(0,3)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       iter = 5000, warmup = 1000, 
                       chains = 3, cores = future::availableCores(),
                       control = list(adapt_delta = 0.999999999999999, max_treedepth = 27))

launch_shinystan(hurd_latcross_3)
waic(hurd_latcross_2,hurd_latcross_3) # latcross 3 fits slightly better
  
#try removing angling type
hurd_latcross_4 <- brm(data = assays, family = hurdle_negbinomial(),
                       bf(Latcross_hurdle  ~ 1  + Trial + ST_length + (1| ID) + (1 | Pool),
                          hu ~ 1 + Trial + ST_length + (1| ID) + (1 | Pool)),
                       prior = c(set_prior("normal(0,3)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       iter = 5000, warmup = 1000, 
                       chains = 3, cores = future::availableCores(),
                       control = list(adapt_delta = 0.99999999999999, max_treedepth = 27))  

waic(hurd_latcross_2,hurd_latcross_3, hurd_latcross_4) # hurd_latcross_4 fits slightly better

#remove length
hurd_latcross_5 <- brm(data = assays, family = hurdle_negbinomial(),
                       bf(Latcross_hurdle  ~ 1  + Trial + (1| ID) + (1 | Pool),
                          hu ~ 1 + Trial + (1| ID) + (1 | Pool)),
                       prior = c(set_prior("normal(0,3)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       iter = 5000, warmup = 1000, 
                       chains = 3, cores = future::availableCores(),
                       control = list(adapt_delta = 0.99999999999999, max_treedepth = 27)) 

waic(hurd_latcross_1, hurd_latcross_2,hurd_latcross_3, hurd_latcross_4, hurd_latcross_5) # hurd_latcross_4 fits slightly better

#get summary data from latcross_4 - best fit model
summary(hurd_latcross_4)
launch_shinystan(hurd_latcross_4Kint)

my_icc_tibble(hurd_latcross_4, ~(1|ID) + (1|Pool), ~(1|Pool), hurdle_sep = T) 

#Adding Fulton's K to best fit model
hurd_latcross_4K <- brm(data = assays, family = hurdle_negbinomial(),
                       bf(Latcross_hurdle  ~ 1  + Trial + ST_length + K_ST + (1| ID) + (1 | Pool),
                          hu ~ 1 + Trial + ST_length + K_ST + (1| ID) + (1 | Pool)),
                       prior = c(set_prior("normal(0,3)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       iter = 5000, warmup = 1000, 
                       chains = 3, cores = future::availableCores(),
                       control = list(adapt_delta = 0.999999999999999, max_treedepth = 29))  

summary(hurd_latcross_4K)

##marginal effects plot
#save marginal effects as an object, give the specific effect to include, dpar= hu is for hurdle portion of model
hl4k_me <- marginal_effects(hurd_latcross_4K, effects = "K_ST", dpar = "hu")
#then plot the effect from the marginal effects object using ggplot
hl4k_me$K_ST %>%
  ggplot(aes(x = K_ST, y = estimate__)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("Probability of\ncrossing barrier") +
  xlab("Fulton's K (scaled)")

#add scaled 13C
hurd_latcross_4C <- brm(data = assays, family = hurdle_negbinomial(),
                        bf(Latcross_hurdle  ~ 1  + Trial + ST_length + X13C_scaled + (1| ID) + (1 | Pool),
                           hu ~ 1 + Trial + ST_length + X13C_scaled + (1| ID) + (1 | Pool)),
                        prior = c(set_prior("normal(0,3)", class = "b"),
                                  set_prior("cauchy(0,2)", class = "sd")),
                        iter = 5000, warmup = 1000, 
                        chains = 3, cores = future::availableCores(),
                        control = list(adapt_delta = 0.999999999999999, max_treedepth = 29)) 
summary(hurd_latcross_4C)

#add scaled 15N
hurd_latcross_4N <- brm(data = assays, family = hurdle_negbinomial(),
                        bf(Latcross_hurdle  ~ 1  + Trial + ST_length + X15N_scaled + (1| ID) + (1 | Pool),
                           hu ~ 1 + Trial + ST_length + X15N_scaled + (1| ID) + (1 | Pool)),
                        prior = c(set_prior("normal(0,3)", class = "b"),
                                  set_prior("cauchy(0,2)", class = "sd")),
                        iter = 5000, warmup = 1000, 
                        chains = 3, cores = future::availableCores(),
                        control = list(adapt_delta = 0.999999999999999, max_treedepth = 29)) 
summary(hurd_latcross_4N)

####Number of barrier crosses analysis####

#Run basic model with hurdle
hurd_numcross_1 <- brm(data = assays, family = hurdle_negbinomial(),
                         bf(Num_cross  ~ 1 + Type + Trial + (1| ID) + (1 | Pool),
                            hu ~ 1 + Type + Trial + (1| ID) + (1 | Pool)),
                         iter = 5000, warmup = 1000, 
                         chains = 3, cores = future::availableCores(),
                         control = list(adapt_delta = 0.99999999999999999, max_treedepth = 27))
  
  summary(hurd_numcross_1 )
  icc(hurd_numcross_1, re.form= ~(1| ID), typical = "mean", prob = 0.95, ppd = TRUE) ##not repeatable
  my_icc_tibble(hurd_numcross_1, ~(1|ID) + (1|Pool), ~(1|Pool), hurdle_sep = T)

#try with full model that includes length and sex
hurd_numcross_2 <- brm(data = assays, family = hurdle_negbinomial(),
                       bf(Num_cross  ~ 1 + Type + Trial + ST_length + Sex_sim + (1| ID) + (1 | Pool),
                          hu ~ 1 + Type + Trial + ST_length + Sex_sim + (1| ID) + (1 | Pool)),
                       prior = c(set_prior("normal(0,3)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       iter = 5000, warmup = 1000, 
                       chains = 3, cores = future::availableCores(),
                       control = list(adapt_delta = 0.999999999, max_treedepth = 27))
launch_shinystan(hurd_numcross_2)

#remove sex
hurd_numcross_3 <- brm(data = assays, family = hurdle_negbinomial(),
                       bf(Num_cross  ~ 1 + Type + Trial + ST_length  + (1| ID) + (1 | Pool),
                          hu ~ 1 + Type + Trial + ST_length + (1| ID) + (1 | Pool)),
                       prior = c(set_prior("normal(0,3)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       iter = 5000, warmup = 1000, 
                       chains = 3, cores = future::availableCores(),
                       control = list(adapt_delta = 0.999999999, max_treedepth = 27))


waic(hurd_numcross_2,hurd_numcross_3) ## numcross 3 fits every so slightly better
launch_shinystan(hurd_numcross_3)

#remove type
hurd_numcross_4 <- brm(data = assays, family = hurdle_negbinomial(),
                       bf(Num_cross  ~ 1 + Trial + ST_length  + (1| ID) + (1 | Pool),
                          hu ~ 1 +  Trial + ST_length + (1| ID) + (1 | Pool)),
                       prior = c(set_prior("normal(0,3)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       iter = 5000, warmup = 1000, 
                       chains = 3, cores = future::availableCores(),
                       control = list(adapt_delta = 0.999999999, max_treedepth = 27))


waic(hurd_numcross_2,hurd_numcross_3, hurd_numcross_4) ## numcross 4 fits every so slightly better

#getting summary info from best fit model
launch_shinystan(hurd_numcross_4)
summary(hurd_numcross_4)
my_icc_tibble(hurd_numcross_4, ~(1|ID) + (1|Pool), ~(1|Pool), hurdle_sep = T)

#Add Fulton's K to best fit model
hurd_numcross_4K <- brm(data = assays, family = hurdle_negbinomial(),
                       bf(Num_cross  ~ 1 + Trial + ST_length  + K_ST + (1| ID) + (1 | Pool),
                          hu ~ 1 +  Trial + ST_length + K_ST + (1| ID) + (1 | Pool)),
                       prior = c(set_prior("normal(0,3)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       iter = 5000, warmup = 1000, 
                       chains = 3, cores = future::availableCores(),
                       control = list(adapt_delta = 0.999999999, max_treedepth = 27))
summary (hurd_numcross_4K)

#Add scaled 13C
hurd_numcross_4C <- brm(data = assays, family = hurdle_negbinomial(),
                        bf(Num_cross  ~ 1 + Trial + ST_length  + X13C_scaled + (1| ID) + (1 | Pool),
                           hu ~ 1 +  Trial + ST_length + X13C_scaled + (1| ID) + (1 | Pool)),
                        prior = c(set_prior("normal(0,3)", class = "b"),
                                  set_prior("cauchy(0,2)", class = "sd")),
                        iter = 5000, warmup = 1000, 
                        chains = 3, cores = future::availableCores(),
                        control = list(adapt_delta = 0.999999999, max_treedepth = 27))
summary (hurd_numcross_4C)

#Add scaled 15N
hurd_numcross_4N <- brm(data = assays, family = hurdle_negbinomial(),
                        bf(Num_cross  ~ 1 + Trial + ST_length  + X15N_scaled + (1| ID) + (1 | Pool),
                           hu ~ 1 +  Trial + ST_length + X15N_scaled + (1| ID) + (1 | Pool)),
                        prior = c(set_prior("normal(0,3)", class = "b"),
                                  set_prior("cauchy(0,2)", class = "sd")),
                        iter = 5000, warmup = 1000, 
                        chains = 3, cores = future::availableCores(),
                        control = list(adapt_delta = 0.999999999, max_treedepth = 27))
summary (hurd_numcross_4N)

####Time spent center of arena analysis####

#Negbinomial - run basic model
timecenter_1 <- brm(data = assays, family = negbinomial(),
                    Center_s  ~ 1 + Type + Trial + (1| ID) + (1 | Pool),
                    iter = 5000, warmup = 1000, 
                    chains = 3, cores = future::availableCores(),
                    control = list(adapt_delta = 0.9999999999999999, max_treedepth = 29))
summary(timecenter_1)
icc(timecenter_1, re.form= ~(1| ID), typical = "mean", prob = 0.95, ppd = TRUE) #not repeatable


#full model with length and sex
timecenter_2 <- brm(data = assays, family = negbinomial(),
                    Center_s  ~ 1 + Type + Trial + ST_length + Sex_sim + (1| ID) + (1 | Pool),
                    iter = 5000, warmup = 1000, 
                    prior = c(set_prior("normal(0,3)", class = "b"),
                              set_prior("cauchy(0,2)", class = "sd")),
                    chains = 3, cores = future::availableCores(),
                    control = list(adapt_delta = 0.9999999999999999, max_treedepth = 29))
  launch_shinystan(timecenter_2)
  
#remove sex
  timecenter_3 <- brm(data = assays, family = negbinomial(),
                      Center_s  ~ 1 + Type + Trial + ST_length + (1| ID) + (1 | Pool),
                      iter = 5000, warmup = 1000, 
                      prior = c(set_prior("normal(0,3)", class = "b"),
                                set_prior("cauchy(0,2)", class = "sd")),
                      chains = 3, cores = future::availableCores(),
                      control = list(adapt_delta = 0.9999999999999999, max_treedepth = 29))
  waic(timecenter_2,timecenter_3)
  launch_shinystan(timecenter_3)
  
#remove length
  timecenter_4 <- brm(data = assays, family = negbinomial(),
                      Center_s  ~ 1 + Type + Trial + (1| ID) + (1 | Pool),
                      iter = 5000, warmup = 1000, 
                      prior = c(set_prior("normal(0,3)", class = "b"),
                                set_prior("cauchy(0,2)", class = "sd")),
                      chains = 3, cores = future::availableCores(),
                      control = list(adapt_delta = 0.9999999999999999, max_treedepth = 29))
  waic(timecenter_2,timecenter_3, timecenter_4)
  launch_shinystan(timecenter_4)
  
#remove angling type
  timecenter_5 <- brm(data = assays, family = negbinomial(),
                      Center_s  ~ 1 + Trial + (1| ID) + (1 | Pool),
                      iter = 5000, warmup = 1000, 
                      prior = c(set_prior("normal(0,3)", class = "b"),
                                set_prior("cauchy(0,2)", class = "sd")),
                      chains = 3, cores = future::availableCores(),
                      control = list(adapt_delta = 0.9999999999999999, max_treedepth = 29))
  
#model 4 fits best 
waic(timecenter_2,timecenter_3, timecenter_4, timecenter_5) 
#summary data of best fit model
summary(timecenter_4)
icc(timecenter_4, re.form= ~(1| ID), typical = "median", prob = 0.95, ppd = TRUE)  
  
#add Fulton's K
  timecenter_4K <- brm(data = assays, family = negbinomial(),
                      Center_s  ~ 1 + Type + Trial + K_ST + (1| ID) + (1 | Pool),
                      iter = 5000, warmup = 1000, 
                      prior = c(set_prior("normal(0,3)", class = "b"),
                                set_prior("cauchy(0,2)", class = "sd")),
                      chains = 3, cores = future::availableCores(),
                      control = list(adapt_delta = 0.9999999999999999, max_treedepth = 29))
  summary(timecenter_4K)
  
#add 13C
  timecenter_4C <- brm(data = assays, family = negbinomial(),
                       Center_s  ~ 1 + Type + Trial + X13C_scaled + (1| ID) + (1 | Pool),
                       iter = 5000, warmup = 1000, 
                       prior = c(set_prior("normal(0,3)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       chains = 3, cores = future::availableCores(),
                       control = list(adapt_delta = 0.9999999999999999, max_treedepth = 29))
  summary(timecenter_4C)
  
#add 15N
  timecenter_4N <- brm(data = assays, family = negbinomial(),
                       Center_s  ~ 1 + Type + Trial + X15N_scaled + (1| ID) + (1 | Pool),
                       iter = 5000, warmup = 1000, 
                       prior = c(set_prior("normal(0,3)", class = "b"),
                                 set_prior("cauchy(0,2)", class = "sd")),
                       chains = 3, cores = future::availableCores(),
                       control = list(adapt_delta = 0.9999999999999999, max_treedepth = 29))
  summary(timecenter_4N)
  

####multivariate 13C and 15N####

  #make new data frame that only uses Trial 1
  assays_t1 <- assays %>% filter(Trial < 2)
  
  #mutivariate model
  multi_iso_norand <- brm(
    mvbind(X13C, X15N) ~ 1  + Pop + Type + set_rescor(T), 
    data = assays_t1,
    iter = 10000, warmup = 2000, chains = 3, cores = 3,
    control = list(adapt_delta = 0.9, max_treedepth = 10))

    multi_iso_norand

    brms::pp_check(multi_iso_norand, resp = "X13C")
    brms::pp_check(multi_iso_norand, resp = "X15N")

####posthoc analysis of differences between catch types####

#make sure data set does not include standardized outcome variables
assays <- read.csv("Grayling_behavioraldata_clean.csv")
assays_t1 <- assays %>% filter(Trial < 2)

#model of length with location and technique as predictors
length_1 <- brm(data = assays_t1, family = gaussian(),
                ST_length  ~ 1 + Pop + Type,
                iter = 5000, warmup = 1000, 
                chains = 3, cores = future::availableCores(),
                control = list(adapt_delta = 0.99999999, max_treedepth = 25))

length_1

#model of body condition with location and technique as predictors
condition_1 <- brm(data = assays_t1, family = gaussian(),
                K_ST*1000  ~ 1 + Pop + Type,
                iter = 5000, warmup = 1000, 
                chains = 3, cores = future::availableCores(),
                control = list(adapt_delta = 0.99999999, max_treedepth = 25))
