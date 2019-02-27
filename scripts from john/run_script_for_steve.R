# Notes from John's email:
# Scripts that include age-based risk factors are attached.

#To get this to work, you will obviously need to change file names paths in the run script.

# For your benefit, I reduced the number of treatment strategies and treatment percentages 
# in the run script to 1.   If you decided that you don't want treatment, you could remove 
# those loops completely.

# Good news is that the routine for stable age distribution scale with population 
# size (i.e., no further tweaking of other parameters with population size required)
# You can also remove treatment by manipulating the following in the param file 
# (e.g., making treatment start at year 20000)

# As for the biological/social assumptions, I think you will figure out the statnet terms quickly enough.  
# There are, however, a few commands related to age-related condom use and probability of sex that you 
# might want to review. I also jacked up the baseline transmission rate to get an appreciable epidemic. 

param_set="Base10000"
param_set_name=paste("params",param_set,sep="")
param_file=paste(param_set_name,".R",sep="")

setwd("H://TasP2/Exp1_Base10000_second_try")
source(param_file) # Default parameters for Tas-By-Age study

# Biological parameters specific to this simulation
exp_number = paste("Exp1_",param_set_name,"_",sep="")

#percent_treated = list(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
percent_treated = list(0.0)
nw <- nw_setup(evoparams) # Sets up the initial network
nw_saved_john <- nw
cat("Finished setting up the model\n")

for (ii in 1:1){ # Loop with treatment targeting strategies (where 1 = "random treatment")
  evoparams$tx_type <- TasP_criteria[ii] # See "TasP_Params_March2018.R" for list of criteria

   for(jj in 1:length(percent_treated)) {
    model_name = paste("txt_opto_",exp_number,"_",paste(TasP_criteria[ii][[1]],collapse="_"),"_prop",as.character(percent_treated[jj]),sep="")
    evoparams$proportion_treated  = percent_treated[jj][[1]]
    
    cat("About to run evorun: (evoparams$mean_test_interval_under25 = ",evoparams$mean_test_interval_under25,"\n")
    
    nw <- nw_saved_john # Don't assume that evorun leaves nw alone
     evomodel <- evorun(modules,evoparams,nw) # Runs the simulation (Line above b/c cannot be sure that evorun doesn't change nw!)
    closeAllConnections()
    
    assign(model_name,evomodel)
    file_name <- paste(model_name,".RData",sep="")
    save(evomodel,file = file.path(getwd(),file_name) )
    evoplot(model=evomodel, name = paste(model_name,".pdf",sep=""), path= file.path(getwd())) # Saves a pdf with standard plots
    remove(evomodel); remove(list=model_name)
    
    #remove(evomodel); remove(list=model_name)
  } # Percent treated
} # TasP criterion


