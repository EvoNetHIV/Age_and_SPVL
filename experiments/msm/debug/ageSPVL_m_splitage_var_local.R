
#library(devtools)
#remotes::install_github("https://github.com/chad-klumb/network", ref="tergmLiteterms")
#remotes::install_github("https://github.com/statnet/tergmLite", ref="tergmLiteterms")
#remotes::install_github("https://github.com/statnet/EpiModel", ref="master")
#install_github("EvoNetHIV/EvoNet",subdir="pkg")  ## OR BUILD FROM LOCAL GIT CLONE

#hyak=T
#hyak_par=T
#
#if(!isTRUE(hyak) & isTRUE(hyak_par)){stop("hyak flags incorrect")}
#
#hyak_path= '/gscratch/csde/sestans/age'
local_path=getwd()
#
#if(hyak)outpath=hyak_path else 
outpath=local_path

#--------------------------------------------------------------
library(evonet)
#library(EpiModelHPC)
#--------------------------------------------------------------

# values of changing parameter
variable_values <- c(75)
name_range <- c('splitage_var_')
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
#for(ii in 1:length(variable_values)){
ii <- 1

  initial_pop = 10000
  mean_sqrt_age_diff = 0.4
  meandeg = 2
  min_age = 18
  max_age = 75
  
  param_list=list(
    
    model_name = variable_names[ii],
    #    nsims = 28, ncores = 28,
    nsims = 1, ncores = 1,
    
    min_age = 18,
    max_age = variable_values[ii],#75,                 
    proportion_treated = 0.5,     
    prob_tx_droput = 0.05,        
    prob_sex_age_19	= 0.8,        
    max_age_sex	= 55,             
    #relation_dur = 1000,          
    susceptibility_var = 0.00,
    n_steps = 365*1,
    
    initial_pop = initial_pop,
    initial_infected = 0.1*initial_pop,

    ##### %<25 = 0.187, 0.188, 0.184, 0.193, 0.185, 0.185, 0.186, 0.181, 0.186, 0.182 --- mean = 0.186
        
    nw_form_terms = "~edges + nodemix(~age>25) + absdiff('sqrt_age') + offset(nodematch('role',diff=TRUE, keep=1:2))",
    target_stats = c(initial_pop*meandeg/2,
                     #initial_pop*meandeg*(max_age-25)/(max_age-min_age),
                     initial_pop*meandeg*0.186*0.814*0.2, initial_pop*meandeg*0.814*0.814*0.8,
                     mean_sqrt_age_diff*initial_pop*meandeg/2),

        nw_coef_form = c(-Inf, -Inf),
    dissolution = "~offset(edges) + offset(nodemix(~age>25))",
    relation_dur = c(50, 316, 2000),
    
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
    #    fast_edgelist = TRUE,
    fast_edgelist = FALSE,
    plot_nw = FALSE,
    save_vl_list = FALSE,
    #    save_network=FALSE
    save_network=TRUE
  )
  
  model_name = model_names[ii] 
  evoparams <- do.call(evonet_setup,param_list)
  evoparams$dissolution <- "~offset(edges) + offset(nodemix(~age>25))"
  nw <- nw_setup(evoparams)
  evomodel <- evorun(modules,evoparams,nw)
  assign(model_name,evomodel)
  file_name <- paste(model_name,".rda",sep="")
  save(list=model_name,
       file = file.path(outpath,file_name) )
  #  remove(evomodel)
  remove(list=model_name)
#}


dd <- networkDynamic::get.edge.activity(evomodel$nw[[1]],
                                        as.spellList = TRUE,
                                        active.default = F)

#take a look
head(dd)

#assign "pop" object (agent attributes) to "aa" for shorthand
aa=evomodel$pop[[1]]

###################################################
#related example, plotting partner ages,fyi

aa=evomodel$pop[[1]]

ff=cbind(
  aa$age[dd$head],
  aa$age[dd$tail])

plot(ff[,1],ff[,2],xlab="age",ylab="age")
# So yes, density in the cumulative network is as expected: 
#   (<25/<25) is very dense, (<25/>25) is moderate, (>25/>25) is very rare.


# dissolution probabilities by type

ix=which(aa$age<25)
i1=which((dd$tail %in% ix) + (dd$head %in% ix) == 2)
i2=which((dd$tail %in% ix) + (dd$head %in% ix) == 1)
i3=which((dd$tail %in% ix) + (dd$head %in% ix) == 0)

plot(dd$onset[i1], pmin(dd$terminus[i1],n_steps+1))
plot(dd$onset[i2], pmin(dd$terminus[i2],n_steps+1))
plot(dd$onset[i3], pmin(dd$terminus[i3],n_steps+1))

days_no_diss <- pmin(dd$terminus,n_steps+1) - pmax(dd$onset,0) - 1
days_diss <- ifelse(dd$terminus==Inf, 0, 1)

dissprob.yy <- sum(days_diss[i1]) / sum((days_no_diss[i1] + days_diss[i1]))
dissprob.yo <- sum(days_diss[i2]) / sum((days_no_diss[i2] + days_diss[i2]))
dissprob.oo <- sum(days_diss[i3]) / sum((days_no_diss[i3] + days_diss[i3]))

1/dissprob.yy
1/dissprob.yo
1/dissprob.oo
