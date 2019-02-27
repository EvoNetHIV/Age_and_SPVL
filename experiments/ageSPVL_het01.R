
library(evonet)
#options(error=recover) # go into debug mode on error
#--------------------------------------------------------------

# items used in calcuation of params
initial_pop      = 100
initial_infected = 10
mean_sqrt_age_diff = 1.2              # CHECK
meandeg = 0.7                         # CHECK

param_list <- list(
  
  # Top-level simulation parameters
  model_name = "h01",
  n_steps          = 365,
  nsims            = 5,
  ncores           = 5,
  
  popsumm_frequency=30,
  fast_edgelist=T,
  plot_nw=F,
  save_vl_list=TRUE,
  
  # Initial conditions
  initial_pop      = initial_pop,
  initial_infected = initial_infected,
  
  # Demographic parameters
  model_sex = "hetero",
  asmr_data_male = "south_africa_male",
  asmr_data_female = "south_africa_female",
  initial_agedata_male = "stable_age_no_hiv_dist",
  initial_agedata_female = "stable_age_no_hiv_dist",
  birth_model = "exponential_growth",
  baseline_input_exp_growth = 0.007*(initial_pop/100), # Calibrated by hand; assumes initial_agedata_male = "stable_age_no_hiv_dist"
  pop_growth_rate_annual   = 0.01,
  min_age = 16,
  max_age = 99,

  # Factors affecting transmission probabilities (note that some may be overwritten by calling program)
  trans_lambda  =  0.000247 * 15, # some value times Hughes et al. (b/c they may be biased)  
  trans_RR_receptive_vaginal = 1.5,

  # Condom use 
  age_condom_use_halves = 35,
  condom_use_age = TRUE,
  condom_prob <- 0.80, # When condom_use_age == TRUE this gives probability that 16-year old uses condom, otherwise this is the probability for everyone
  condom_use_rel_dur = FALSE,
  RR_cond_male_concurrent  = 1.0,
  RR_cond_fem_concurrent   = 1.0,
  
  # Coital acts
  prob_sex_by_age          = TRUE,
  prob_sex_age_19          = 0.20, # used when prob_sex_by_age == TRUE
  circum_prob              = 0.4,
  prop_AI                  = 0.0,
  mean_prop_acts_AI        = 0.0,
  max_age_sex	             = 55,  # From Steve's code, not John's : check on this carefully
  mean_sex_acts_day        = 0.2, # From Steve's code, not John's : check on this carefully
  
  # Network : to finish
  #nw_coef_form = c(-Inf, -Inf),
  relation_dur = 1000,
  nw_form_terms = "~edges+absdiff('sqrt_age') + offset(nodematch('role',diff=TRUE, keep=1:2))",
  target_stats = c(initial_pop*meandeg/2, mean_sqrt_age_diff*initial_pop*meandeg/2),
  
  # Testing: to finish
  
  mean_test_interval_male  = 1,

  
  # Treatment: to finish

  start_treatment_campaign      = 3,
  prob_eligible_ART        = 1.0,
  tx_limit = "absolute_num",
  vl_full_supp    = 10,
  tx_schedule_props  = c("F"=0,"V"=1.0,"N"=0,"P"=0), # Percentage of people who always (F), sometimes (V), or never (N) take therapy
  prob_tx_droput  = 0.05,   # Probability of "V"-type agent discontinuing therapy over the course of a year
  yearly_incr_tx  = 0.0, #  No annual increase in the number of people being treated after the TasP campaign
  proportion_treated  = 0.5,
  tx_type = "random",        # From Steve
  mean_trtmnt_delay = 0,     # From Steve
  #start_treat_before_big_campaign = 0, # Start time of gradual ramp-up prior to the start of the big TasP campaign (this part random)
  #proportion_treated_begin = 0.5, # Treated before gradual ramp-up 
  #prop_tx_before=0.5,
  #prob_care= 1.0,
  #yearly_incr_tx  = 0.0, #  increase in # being treated after the TasP campaign
  
  
  #testing params
  testing_model            = "interval"

#  "diagnosis.FUN"      = social_testing_diagnosis_module, #change from PrEP sims
#  "treatment.FUN"      = social_treatment_module_john_v3,

)

evoparams <- do.call(evonet_setup,param_list)
nw <- nw_setup(evoparams)

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
  "deaths",
  "births",
  "summary_module")

evomodel <- evorun(modules,evoparams,nw)
ageSPVL_h01 <- evomodel
save(ageSPVL_h01, file="ageSPVL_h01.rda")

#aa=evomodel
#b=aa$pop[[1]]  
#bb$tx_stop_time
#bb$treated
#bb$tx_init_time
#table(bb$tx_schedule)

if(evoparams$save_vl_list==TRUE){
  plot_vl_trajectories(model=evomodel,sim=1,
                       outpath=evoparams$output_path,
                       name="ageSPVL_vl_m42")
} 
