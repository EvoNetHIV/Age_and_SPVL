
age1 <- dat$attr$age[dat$el[[1]][,1]]
age2 <- dat$attr$age[dat$el[[1]][,2]]
  
plot(age1, age2)
mean(abs(sqrt(age1) - sqrt(age2)))                      # 1.191 (should be ~1.2, yay)
mean(abs(age1-age2))                                    # 15 (wow!)

# to be used in nodecov
mean(c(age1,age2))                                      #  41.75

# to be used with age>25      with ~edges
mean(age1<=25 & age2<=25)                               # 0.0233
mean(age1<=25 & age2>25) + mean(age1>25 & age2<=25)     # 0.2064
mean(age1>25 & age2>25)                                 # 0.7703

# to be used with age>25      with ~edges + absdiff(sqrt.age)
mean(age1<=25 & age2<=25)                               # 0.0406
mean(age1<=25 & age2>25) + mean(age1>25 & age2<=25)     # 0.2746
mean(age1>25 & age2>25)                                 # 0.6848

summary(nw$fit$network ~ edges + nodemix(~age>25, levels2=2:3))
#edges mix.age>25.FALSE.TRUE  mix.age>25.TRUE.TRUE            
#3500                   712                  2714   ~edges
# which equals 0.203 and 0.775 for proportion FT and tt, respectively

#3500                   984                  2364   ~edges + absdiff(sqrt_age)
# which equals 0.281 and 0.675 for proportion FT and tt, respectively


summary(nw$fit$network ~ edges + nodemix(~age>25, levels2=2:3))
summary(nw$fit$network ~ edges + nodecov('age'))

summary(nw$fit$network ~ edges + concurrent)
# 3500, 1526 for run 206

summary(nw$fit$network ~ edges + concurrent(by=~age>25))
sum(nw$fit$network%v% 'age'>25)
# 3500, 266, 1260 for run 206
# 1827 nodes < 25, 8173 nodes > 25
#so concurrent = ~15% for both group 


