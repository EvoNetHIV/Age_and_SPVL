
library(evonet)
options(error=browser) # go into debug mode on error
#--------------------------------------------------------------

initial_pop = 1000
mean_sqrt_age_diff = 0.6
meandeg = 0.7

param_list=list(
  model_name = "m20",
  nw_form_terms = "~edges+absdiff('sqrt_age') + offset(nodematch('role',diff=TRUE, keep=1:2))",
  target_stats = c(initial_pop*meandeg/2, mean_sqrt_age_diff*initial_pop*meandeg/2),
  #mean_sex_acts_day = 0.2,
  min_age = 18,
  max_age = 55,

  #age_dist = seq(50, 10, -10/9)/1110,
  initial_agedata_male = "linear_decrease",
  nw_coef_form = c(-Inf, -Inf),
  prob_sex_by_age	= TRUE,
  prob_sex_age_19	= 0.4,
  max_age_sex	= 55,
  relation_dur = 200,
  
  tx_type = "random",
  mean_trtmnt_delay = 1,
  start_treatment_campaign = 1,
  proportion_treated = 0.5,
  tx_schedule_props  = c("F"=0,"V"=1.0,"N"=0,"P"=0), # Percentage of people who always (F), sometimes (V), or never (N) take therapy
  
  #testing params
  testing_model            = "interval",
  mean_test_interval_male  = 365,
  
  prob_care                = 1.0,  # Percent population that could get treated given "all out" campaign
  prob_eligible_ART        = 1.0,
  tx_limit = "absolute_num", # Choices: "absolute_num" or percentage"
  vl_full_supp    = 1e-3, # Assume highly suppressiver therapy
  
#  "diagnosis.FUN"      = social_testing_diagnosis_module, #change from PrEP sims
#  "treatment.FUN"      = social_treatment_module_john_v3,

  nsims = 1,
  initial_pop = initial_pop,
  initial_infected = 100,
  n_steps = 365*20,
  popsumm_frequency=10,
  fast_edgelist=T,
  plot_nw=F,
  save_vl_list=TRUE

)

evoparams <- do.call(evonet_setup,param_list)
nw <- nw_setup(evoparams)

evoparams$prob_tx_droput = 0.1    # Yes, the parameter name is missing an 'o'

#evoparams$prob_care                = 1.0  # Percent population that could get treated given "all out" campaign
#evoparams$prob_eligible_ART        = 1.0
#evoparams$tx_limit = "absolute_num" # Choices: "absolute_num" or percentage"
#evoparams$vl_full_supp    = 1e-3 # Assume highly suppressiver therapy

modules <- c(
  "aging",
  "testing",
  "treatment_dropout",
  "targeted_treatment2",
  "treatment",
  "viral_update_delayed_rebound",
  "cd4_update2",
  "coital_acts",
  "transmission",
  "deaths",
  "births",
  "summary_module"
)

evomodel <- evorun(modules,evoparams,nw)
ageSPVL_m20 <- evomodel
save(ageSPVL_m20, file="ageSPVL_m20.rda")

