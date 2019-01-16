
library(evonet)
options(error=browser) # go into debug mode on error
#--------------------------------------------------------------

initial_pop = 1000
mean_sqrt_age_diff = 0.6
meandeg = 0.7

param_list=list(
  model_name = "m14",
  nw_form_terms = "~edges+absdiff('sqrt_age') + offset(nodematch('role',diff=TRUE, keep=1:2))",
  target_stats = c(initial_pop*meandeg/2, mean_sqrt_age_diff*initial_pop*meandeg/2),
  #mean_sex_acts_day = 0.2,
  min_age = 18,
  max_age = 55,

  #age_dist = seq(50, 10, -10/9)/1110,
  initial_agedata_male = "linear_decrease",
  nw_coef_form = c(-Inf, -Inf),
  prob_sex_by_age	= TRUE,
  prob_sex_age_19	= 0.3,
  max_age_sex	= 55,
  relation_dur = 50,
  
  nsims = 1,
  initial_pop = initial_pop,
  initial_infected = 100,
  n_steps = 365*20,
  popsumm_frequency=30,
  fast_edgelist=T,
  plot_nw=F
)

evoparams <- do.call(evonet_setup,param_list)
nw <- nw_setup(evoparams)
modules <- c(
  "aging",
  "testing",
  "treatment",
  "viral_update",
  "coital_acts",
  "transmission",
  "deaths",
  "births",
  "summary_module"
)

evomodel <- evorun(modules,evoparams,nw)
ageSPVL_m14 <- evomodel
save(ageSPVL_m14, file="ageSPVL_m14.rda")

