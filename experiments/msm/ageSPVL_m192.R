library(evonet)

initial_pop = 1000
mean_sqrt_age_diff = 1.2
meandeg = 0.7
min_age = 18
max_age = 75

param_list=list(
  model_name = "m192",
#  nsims = 64, ncores = 16,
  nsims = 5, ncores = 1,

  min_age = min_age,
  max_age = max_age,
  proportion_treated = 0.5,
  prob_tx_droput = 0.05,
  prob_sex_age_19	= 0.8,
  max_age_sex	= 55,
  #relation_dur = 1000,  # see below
  susceptibility_var = 0,
  n_steps = 365*50,
  
  initial_pop = initial_pop,
  initial_infected = 0.1*initial_pop,

  nw_form_terms = "~edges + nodefactor(~age>25) + absdiff('sqrt_age') + offset(nodematch('role',diff=TRUE, keep=1:2))",
  target_stats = c(initial_pop*meandeg/2, initial_pop*meandeg*(max_age-25)/(max_age-min_age),  
                   mean_sqrt_age_diff*initial_pop*meandeg/2),
  nw_coef_form = c(-Inf, -Inf),
  dissolution <- "~offset(edges) + offset(nodefactor(~age>25))",
  relation_dur=c(50, 2000),
  
  age_dist_new_adds = "min_age",
  initial_agedata_male = "stable_age_no_hiv_dist",
  birth_model = "exponential_growth",
  baseline_input_exp_growth = 0.007*initial_pop/100,
  
  prob_sex_by_age	= TRUE,
    
  tx_type = "random",
  mean_trtmnt_delay = 0,
  start_treatment_campaign = 3,
  prob_eligible_ART = 1.0,
  tx_limit = "percent_of_currently_diagnosed",
  vl_full_supp = 10,
  tx_schedule_props = c("F"=0,"V"=1.0,"N"=0,"P"=0),
  yearly_incr_tx = 0.0,
  
  testing_model = "interval",
  mean_test_interval_male  = 1,
  
  popsumm_frequency = 30,
  fast_edgelist = TRUE,
  plot_nw = FALSE,
  save_vl_list = FALSE
)

evoparams <- do.call(evonet_setup,param_list)

evoparams$dissolution <- "~offset(edges) + offset(nodefactor(~age>25))"
evoparams$hpc <- FALSE

nw <- nw_setup(evoparams)

modules <- c(
  "aging",
  "testing",
  "treatment_dropout",
  "targeted_treatment2",
  "viral_update_delayed_rebound",
  "cd4_update2",
  "coital_acts",
  "transmission",
  "evo_departures",
  "evo_arrivals",
  "summary_module")

evomodel <- evorun(modules,evoparams,nw)
ageSPVL_m192 <- evomodel
save(ageSPVL_m192, file="../AgeAndSPVL_oversize/ageSPVL_m192.rda")

##### TEMP

obj <- ageSPVL_m192

popatts[[i]] <- obj$pop
popsumm[[i]] <- get(paste("ageSPVL_m",filler,i,sep=""))$popsumm

agecoef.list[[i]] <- iSPVL.list[[i]] <- ageinf.list[[i]] <- prev.list[[i]] <- 
  numinc.list[[i]] <- agematch.list[[i]] <- meanageinf.list[[i]] <- vector()

for (j in 1:length(popatts[[i]])) {
  agecoef.list[[i]][j] <- lm(log(popatts[[i]][[j]]$SetPoint[popatts[[i]][[j]]$Time_Inf>0],10)~
                               popatts[[i]][[j]]$age_infection[popatts[[i]][[j]]$Time_Inf>0])$coef[[2]]
  iSPVL.list[[i]][j] <- mean(log10(popatts[[i]][[j]]$SetPoint[popatts[[i]][[j]]$Time_Inf>0]),na.rm=TRUE)
  ageinf.list[[i]][j] <- mean(popatts[[i]][[j]]$age_infection,na.rm=TRUE)
  prev.list[[i]][j] <- tail(popsumm[[i]][[j]]$prevalence,1)
  numinc.list[[i]][j] <- sum(popatts[[i]][[j]]$Time_Inf>0, na.rm=TRUE)
  agematch.list[[i]][j] <- obj$nwparam[[1]]$coef.form['absdiff.sqrt_age']
  meanageinf.list[[i]][j] <- mean(popatts[[i]][[j]]$age_infection[popatts[[i]][[j]]$Time_Inf>0], na.rm=TRUE)
}
