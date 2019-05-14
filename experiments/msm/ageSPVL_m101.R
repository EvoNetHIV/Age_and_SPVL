library(evonet)

initial_pop = 100
mean_sqrt_age_diff = 1.2
meandeg = 0.7

param_list=list(
  model_name = "m101",
  nsims = 64, ncores = 16,

  min_age = 18,
  max_age = 75,
  proportion_treated = 0.5,
  prob_tx_droput = 0.05,
  prob_sex_age_19	= 0.8,
  max_age_sex	= 55,
  relation_dur = 1000,
  susceptibility_var = 1,
  n_steps = 365*50,
  
  initial_pop = initial_pop,
  initial_infected = 0.1*initial_pop,

  nw_form_terms = "~edges+absdiff('sqrt_age') + offset(nodematch('role',diff=TRUE, keep=1:2))",
  target_stats = c(initial_pop*meandeg/2, mean_sqrt_age_diff*initial_pop*meandeg/2),
  nw_coef_form = c(-Inf, -Inf),
  
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
ageSPVL_m101 <- evomodel
save(ageSPVL_m101, file="../AgeAndSPVL_oversize/ageSPVL_m101.rda")
