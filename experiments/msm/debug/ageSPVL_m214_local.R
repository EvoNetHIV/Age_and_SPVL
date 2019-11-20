#hyak=T
#hyak_par=T

#if(!isTRUE(hyak) & isTRUE(hyak_par)){stop("hyak flags incorrect")}

#hyak_path= '/gscratch/csde/sestans/age'
local_path=getwd()

#if(hyak)outpath=hyak_path else 
  outpath=local_path

#--------------------------------------------------------------
library(evonet)
#library(EpiModelHPC)
#--------------------------------------------------------------

# values of changing parameter
variable_values <- c(75)
name_range <- c('214_test')
variable_names = paste("m",name_range,"",sep="")
model_names = paste("ageSPVL_",variable_names,"a",sep="")

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

#Run model for each value of specified target stats
for(ii in 1:length(variable_values)){
  initial_pop = 10000
  mean_sqrt_age_diff = 1.2       
  meandeg = 0.8
  
  param_list=list(
    
    model_name = variable_names[ii],
    #nsims = 28, ncores = 28,
    nsims = 10, ncores = 1,
    
    max_age = variable_values[ii],#75,                 
    proportion_treated = 0.5,     
    prob_tx_droput = 0.2,        

    #prob_sex_by_age	= TRUE,
    #prob_sex_age_19	= 0.8,        
    #max_age_sex	= 45,

    prob_sex_func = TRUE,
    coital_coef_base =  0.00,
    coital_coef_mult =  1.55,
    coital_coef_age  = -0.16, 
    coital_coef_dur  = 0, # unused still
    
    #coital_coef_base =   1.4,
    #coital_coef_age  = -0.08, 
    #coital_coef_dur  = -0.08, 

    relation_dur = 500,          
    susceptibility_var = 0.00,
    n_steps = 365*50,
    
    initial_pop = initial_pop,
    initial_infected = 0.2*initial_pop,
    
    nw_form_terms = "~edges+absdiff('sqrt_age') + offset(nodematch('role',diff=TRUE, keep=1:2))",
    target_stats = c(initial_pop*meandeg/2, mean_sqrt_age_diff*initial_pop*meandeg/2),
    nw_coef_form = c(-Inf, -Inf),
    
    age_dist_new_adds = "min_age",
    initial_agedata_male = "stable_age_no_hiv_dist",
    birth_model = "exponential_growth",
    baseline_input_exp_growth = 0.007*initial_pop/100,
    
    
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
  
  model_name = model_names[ii] 
  evoparams <- do.call(evonet_setup,param_list)
  nw <- nw_setup(evoparams)
  evomodel <- evorun(modules,evoparams,nw)
  assign(model_name,evomodel)
  file_name <- paste(model_name,".rda",sep="")
  save(list=model_name,
       file = file.path(outpath,file_name) )
  remove(evomodel)
}  

