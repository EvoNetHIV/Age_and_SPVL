

################################################################
# T&C Intro
#
# From Tanfer and Cummins 1992
# Table 5, White model 2 and Black model 2
# Coital frequency, single women in current relationship)
# Additjonal covariates set to reference cateogires are:
#   living arrangements (not dorm, not with parents)
#   religion (not Catholic, not Conserv Prot)
#   church attendance (<1 per wk)
#   prev. sex rel (0)
#   cohab w/ partner (no)
#   engaged (no)
#   partner rel (not Catholic, not Conserv Prot)
#   religion homog
#   preg risk (using effective birth control)

coef_constant     = c(17.477, 7.475)
coef_age	        = c(0.071, 0.235)
coef_age_1st_sex	= c(-0.511, -0.375)
coef_rel_dur	    = c(-0.121, -0.024)
coef_p_age	      = c(-0.38, -0.49)
coef_age_heterog  = c(1.616, -0.414)
coef_rel_dur_sq   = c(0.001,	0.0001)


################################################################
# T&C Version 1: age and dur as continuous numerical variables

n_years <- 10

init_ego_age <- c(25, 25, 30, 30)
init_par_age <- c(25, 30, 30, 25)
age_1st_sex  <- c(20, 20, 20, 20)

n_models <- length(init_ego_age)

white_cf <- black_cf <- matrix(NA, n_years, n_models)
for (i in 1:n_years) {
  rel_dur <- i - 1
  
  white_cf[i,] <- coef_constant[1] + 
                  coef_age[1] * (init_ego_age+rel_dur) + 
                  coef_age_1st_sex[1] * age_1st_sex +
                  coef_rel_dur[1] * rel_dur +   
                  coef_p_age[1] * (init_par_age+rel_dur) +   
                  coef_age_heterog[1] + (init_par_age-init_ego_age > 2) + 
                  coef_rel_dur_sq[1] * rel_dur^2

  black_cf[i,] <- coef_constant[2] + 
                  coef_age[2] * (init_ego_age+rel_dur) + 
                  coef_age_1st_sex[2] * age_1st_sex +
                  coef_rel_dur[2] * rel_dur +   
                  coef_p_age[2] * (init_par_age+rel_dur) +   
                  coef_age_heterog[2] + (init_par_age-init_ego_age > 2) + 
                  coef_rel_dur_sq[2] * rel_dur^2
}

matplot(cbind(white_cf,black_cf), type='b')

##################################################################
#
# T&C version 1: age and dur as continuous numerical variables

n_years <- 10

init_ego_age <- c(25, 25, 30, 30)
init_par_age <- c(25, 30, 30, 25)
age_1st_sex  <- c(20, 20, 20, 20)

n_models <- length(init_ego_age)

white_cf <- black_cf <- matrix(NA, n_years, n_models)
for (i in 1:n_years) {
  rel_dur <-    ((i-1)>0.5 & (i-1)<=1) + 
              2*((i-1)>1 & (i-1)<=2) + 
              3*((i-1)>2)
  
  
  white_cf[i,] <- coef_constant[1] + 
    coef_age[1] * ((init_ego_age+rel_dur) > 24) + 
    coef_age_1st_sex[1] * ((age_1st_sex>=16) + (age_1st_sex>=19)) +
    coef_rel_dur[1] * rel_dur +   
    coef_p_age[1] * (((init_par_age+rel_dur)>=25) + ((init_par_age+rel_dur)>=30)) +   
    coef_age_heterog[1] + (init_par_age-init_ego_age > 2) + 
    coef_rel_dur_sq[1] * rel_dur^2
  
  black_cf[i,] <- coef_constant[2] + 
    coef_age[2] * ((init_ego_age+rel_dur) > 24) + 
    coef_age_1st_sex[2] * ((age_1st_sex>=16) + (age_1st_sex>=19)) +
    coef_rel_dur[2] * rel_dur +   
    coef_p_age[2] * (((init_par_age+rel_dur)>=25) + ((init_par_age+rel_dur)>=30)) +   
    coef_age_heterog[2] + (init_par_age-init_ego_age > 2) + 
    coef_rel_dur_sq[2] * rel_dur^2
}

matplot(cbind(white_cf), type='b')
matplot(cbind(black_cf), type='b')

##################################################################
# Constructed example

rm(list=ls())
coef_constant <- 15
coef_age <- 0
coef_p_age <- 0
coef_rel_dur <- -0.4
coef_rel_dur_sq <- 0.03

n_years <- 25

init_ego_age <- c(25, 25, 30, 30)
init_par_age <- c(25, 30, 30, 25)

n_models <- length(init_ego_age)

cf <- matrix(NA, n_years, n_models)
for (i in 0:(n_years)) {
  rel_dur=i-1
  cf[i,] <- coef_constant + 
            coef_age * (init_ego_age+rel_dur) + 
            coef_rel_dur * rel_dur +   
            coef_p_age * (init_par_age+rel_dur) +   
            coef_rel_dur_sq * rel_dur^2
}

matplot(cf, type='b')

###############################################3

plot(-1, -1, xlim=c(0,10), ylim=c(0,10), xlab='reldur', ylab='coital freq per mo')
segments(0, 9.0, 0.5, 9.0)
segments(0.5, 8.2, 1, 8.2)
segments(1, 7.0, 2, 7.0)
segments(2, 5.5, 10, 5.5)



################################################3

coef_dur <- -0.3
coef_age_ego <- -0.02
coef_age_alt <- -0.02
freq_new_18 <- 0.3

init_age_ego <- 18
init_age_alt <- 18

rel_durs <- 0:20
age_ego <- (init_age_ego-18) + rel_durs
age_alt <- (init_age_alt-18) + rel_durs
scalar <- freq_new_18/0.5

cf_logodds <- coef_dur*rel_durs + coef_age_ego*age_ego + coef_age_alt*age_alt
cf_prob <- scalar*exp(cf_logodds)/(1+exp(cf_logodds))

plot(rel_durs, cf_prob)

