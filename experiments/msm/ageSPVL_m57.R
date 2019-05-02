
library(evonet)
#options(error=recover) # go into debug mode on error
#--------------------------------------------------------------

initial_pop = 1000
mean_sqrt_age_diff = 0.7
meandeg = 0.7

param_list=list(
  model_name = "m57",
  nw_form_terms = "~edges+absdiff('sqrt_age') + offset(nodematch('role',diff=TRUE, keep=1:2))",
  target_stats = c(initial_pop*meandeg/2, mean_sqrt_age_diff*initial_pop*meandeg/2),
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
  testing_model = "interval",
  mean_test_interval_male = 1,
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

# These could probably all be moved into theparam_list above, but I inhereted 
#   much of this code from various places and don't want to break what works.
evoparams$baseline_input_exp_growth = 0.007*(evoparams$initial_pop/100) 
evoparams$start_treatment_campaign = 3
evoparams$prob_eligible_ART = 1.0
evoparams$tx_limit = "percent_of_currently_diagnosed"      # New option as of 4/30/2019
evoparams$vl_full_supp = 10
evoparams$tx_schedule_props = c("F"=0,"V"=1.0,"N"=0,"P"=0) # Percentage of people who always (F), sometimes (V), or never (N) take therapy
evoparams$prob_tx_droput = 0.05                            # Probability of "V"-type agent discontinuing therapy over the course of a year
evoparams$yearly_incr_tx = 0.0                             # No annual increase in the number of people being treated after the TasP campaign
evoparams$proportion_treated = 0.5

modules <- c(
  "aging",
  "testing",
  "treatment_dropout",
  "targeted_treatment2",
  "viral_update_delayed_rebound",
  "cd4_update2",
  "coital_acts",
  "transmission",
  "deaths",
  "births",
  "summary_module"
  )

evomodel <- evorun(modules,evoparams,nw)
if(evoparams$save_vl_list==TRUE){
  plot_vl_trajectories(model=evomodel,sim=1,
                       outpath="experiments/msm",
                       name="ageSPVL_vl_m57")
} 

ageSPVL_m57 <- evomodel
save(ageSPVL_m57, file="../AgeAndSPVL_oversize/ageSPVL_m57.rda")
