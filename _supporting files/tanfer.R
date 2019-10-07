

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


n_years <- 50

coef_constant    = c(17.477, 7.475)
coef_age	        = c(0.071, 0.235)
coef_age_1st_sex	= c(-0.511, -0.375)
coef_rel_dur	    = c(-0.121, -0.024)
coef_p_age	      = c(-0.38, -0.49)
coef_age_heterog = c(1.616, -0.414)
coef_rel_dur_sq  = c(0.001,	0.0001)

init_ego_age <- c(25, 25, 30, 30)
init_par_age <- c(25, 30, 30, 25)
#age_1st_sex  <- c(20, 20, 20, 20)
age_1st_sex  <- c(2, 2, 2, 2)

n_models <- length(init_ego_age)

white_cf <- black_cf <- matrix(NA, n_years, n_models)
for (i in 1:n_years) {
  rel_dur <- ((i-1)>0.5 & (i-1)<=1) + 2*((i-1)>1 & (i-1)<=2) + 3*((i-1)>2)

  white_cf[i,] <- coef_constant[1] + 
                  coef_age[1] * ((init_ego_age+rel_dur) > 24) + 
                  coef_age_1st_sex[1] * ((age_1st_sex>16) + (age_1st_sex>18)) +
                  coef_rel_dur[1] * rel_dur +   
                  coef_p_age[1] * (((init_par_age+rel_dur)>24) + ((init_par_age+rel_dur)>29)) +   
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

matplot(white_cf)


##################################################################
#

coef_age <- -1
coef_p_age <- -1
coef_dur <- -1
coef_dur^2 <- 0.1



