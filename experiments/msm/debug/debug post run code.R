#source("experiments/msm/debug/ageSPVL_m211_local.R"); source("experiments/msm/debug/ageSPVL_m212_local.R"); source("experiments/msm/debug/ageSPVL_m213_local.R"); source("experiments/msm/debug/ageSPVL_m214_local.R");

i <- 219
runname <- paste('ageSPVL_m', i, '_testa', sep='')
obj <- get(runname)


popatts[[i]] <- obj$pop
popsumm[[i]] <- obj$popsumm
agecoef.list[[i]] <- iSPVL.list[[i]] <- ageinf.list[[i]] <- prev.list[[i]] <-
numinc.list[[i]] <- agematch.list[[i]] <- meanageinf.list[[i]] <- vector()
for (j in 1:length(popatts[[i]])) {
  agecoef.list[[i]][j] <- lm(log(popatts[[i]][[j]]$SetPoint[popatts[[i]][[j]]$Time_Inf>0],10)~
  popatts[[i]][[j]]$age_infection[popatts[[i]][[j]]$Time_Inf>0])$coef[[2]]
  iSPVL.list[[i]][j] <- mean(log10(popatts[[i]][[j]]$SetPoint[popatts[[i]][[j]]$Time_Inf>0]),na.rm=TRUE)
  ageinf.list[[i]][j] <- mean(popatts[[i]][[j]]$age_infection,na.rm=TRUE)
  prev.list[[i]][j] <- tail(popsumm[[i]][[j]]$prevalence,1)
  numinc.list[[i]][j] <- sum(popatts[[i]][[j]]$Time_Inf>0, na.rm=TRUE)
  agematch.list[[i]][j] <- obj$nwparam[[1]]$coef.form['absdiff.sqrt_age']
  meanageinf.list[[i]][j] <- mean(popatts[[i]][[j]]$age_infection[popatts[[i]][[j]]$Time_Inf>0], na.rm=TRUE)
  }

boxplot(prev.list,xlim=c(201,220))
abline(h=0)

boxplot(agecoef.list,xlim=c(201,220))
abline(h=0)

################################ In debug mode
#plot(dat$pop$age[dat$discord_coital_df$inf_id], dat$pop$age[dat$discord_coital_df$sus_id])
evoplot(obj, name=runname)

