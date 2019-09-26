
library(ergm)

n <- 100
aaa <- network.initialize(n, directed=FALSE)
aaa %v% 'age' <- runif(n, 18, 65)
bbb <- ergm(aaa~edges+absdiff('age'), target.stats = c(n*0.3, n*0.3*2))
ccc <- simulate(bbb, nsim = 1)            
