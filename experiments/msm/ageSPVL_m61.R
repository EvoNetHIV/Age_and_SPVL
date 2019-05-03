
library(evonet)
#options(error=recover) # go into debug mode on error
#--------------------------------------------------------------

initial_pop = 1000
mean_sqrt_age_diff = 1.2
meandeg = 0.7

param_list=list(
  model_name = "m61",
  nw_form_terms = "~edges+absdiff('sqrt_age') + offset(nodematch('role',diff=TRUE, keep=1:2))",
  target_stats = c(initial_pop*meandeg/2, mean_sqrt_age_diff*initial_pop*meandeg/2),
  #mean_sex_acts_day = 0.2,
  min_age = 18,
  max_age = 55,
  age_dist_new_adds = "min_age",
  initial_agedata_male = "stable_age_no_hiv_dist",
  birth_model = "exponential_growth",
  
  nw_coef_form = c(-Inf, -Inf),
  prob_sex_by_age	= TRUE,
  prob_sex_age_19	= 0.8,
  max_age_sex	= 55,
  relation_dur = 1000,
  susceptibility_var = 1,
    
  tx_type = "random",
  mean_trtmnt_delay = 0,
  start_treatment_campaign = 1,
  #proportion_treated = 0.2,
  
  #testing params
  testing_model            = "interval",
  mean_test_interval_male  = 1,
  
#  "diagnosis.FUN"      = social_testing_diagnosis_module, #change from PrEP sims
#  "treatment.FUN"      = social_treatment_module_john_v3,

  nsims = 10, ncores = 10,
  initial_pop = initial_pop,
  initial_infected = 200,
  n_steps = 365*50,
  popsumm_frequency=30,
  fast_edgelist=T,
  plot_nw=F,
  save_vl_list=TRUE

)

evoparams <- do.call(evonet_setup,param_list)
nw <- nw_setup(evoparams)

evoparams$baseline_input_exp_growth = 0.007*(evoparams$initial_pop/100) 
evoparams$start_treatment_campaign      = 3
evoparams$prob_eligible_ART        = 1.0
evoparams$tx_limit = "percent_of_currently_diagnosed"      # New option as of 4/30/2019
evoparams$vl_full_supp    = 10
evoparams$tx_schedule_props  = c("F"=0,"V"=1.0,"N"=0,"P"=0) # Percentage of people who always (F), sometimes (V), or never (N) take therapy
evoparams$prob_tx_droput  = 0.05   # Probability of "V"-type agent discontinuing therapy over the course of a year
#evoparams$proportion_treated_begin = 0.5 # Treated before gradual ramp-up 
#evoparams$prop_tx_before=0.5
evoparams$yearly_incr_tx  = 0.0 #  No annual increase in the number of people being treated after the TasP campaign

#### THIS IS IT
evoparams$proportion_treated  = 0.5


if(F) {
  evoparams$proportion_treated  = 0.5
  evoparams$prob_tx_droput  = 0.10   # Probability of "V"-type agent discontinuing therapy over the course of a year
  evoparams$proportion_treated_begin = 0.5 # Treated before gradual ramp-up 
  evoparams$prop_tx_before=0.5
  evoparams$start_treat_before_big_campaign = 0 # Start time of gradual ramp-up prior to the start of the big TasP campaign (this part random)
  evoparams$start_treatment_campaign      = 365*100   # Start time of big TasP campaign (this part depends on tx_type)
  #  Other parameters affecting treatment
  evoparams$prob_care                = 1.0  # Percent population that could get treated given "all out" campaign
  evoparams$prob_eligible_ART        = 1.0
  evoparams$tx_limit = "absolute_num" # Choices: "absolute_num" or percentage"
  evoparams$vl_full_supp    = 10
  evoparams$tx_schedule_props  = c("F"=0,"V"=1.0,"N"=0,"P"=0) # Percentage of people who always (F), sometimes (V), or never (N) take therapy
}
#evoparams$yearly_incr_tx  = 0.0 #  No annual increase in the number of people being treated after the TasP campaign
#evoparams$mean_trtmnt_delay = 15    # Delay betweent diagnosis and start of treatment

modules <- c(
  "aging",
  "testing",
  #"adherence",
  "treatment_dropout",
  "targeted_treatment2",
  "viral_update_delayed_rebound",
  "cd4_update2", # viral_update_cd4_daily
  "coital_acts",
  "transmission",
  "evo_departures",
  "evo_arrivals",
  "summary_module")

evomodel <- evorun(modules,evoparams,nw)
ageSPVL_m61 <- evomodel
save(ageSPVL_m61, file="../AgeAndSPVL_oversize/ageSPVL_m61.rda")

aa=evomodel
bb=aa$pop[[1]]  
bb$tx_stop_time
#edit(treatment_dropout)
bb$treated
bb$tx_init_time
table(bb$tx_schedule)

if(evoparams$save_vl_list==TRUE){
  plot_vl_trajectories(model=evomodel,sim=1,
                       outpath="experiments/msm",
                       name="ageSPVL_vl_m61")
} 
