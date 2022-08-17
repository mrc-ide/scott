# `scott`: Structured Coalescent Transmission Tree simulation 

Example usage: 

```

# Simulate epidemic under exponential growth 
set.seed( 1111 ) 

library(phydynR)
births <- c(I = 'parms$beta * I' )
deaths <- c(I = 'parms$gamma * I' )
dm <- build.demographic.process(births=births
  , deaths = deaths
  , parameterNames=c('beta', 'gamma') 
  , rcpp=FALSE
  , sde = TRUE)
tfgy = dm( theta = list( beta = 1.5, gamma = 1 ), x0 = c(I = 1 ) , t0 = 0, t1 = 10  )

# Simulate transmission tree
library(scott)
n <- 75 # sample size
s <- rep(10, n ) # sample times 
X <- matrix( c(rep(1, n), rep(0,n)), byrow=FALSE , ncol = 2) # sample states
min_dt = 10 / 1e4 # controls precision of simulation 
o = scott(
  tfgy
  , s
  , X 
  , min_dt
  , event_tol = .1 
)

# Here is the table of transmission events within the sample 
(w <- o$waifw)

```

# Roadmap 2022-08-17

- Convenient subroutines for faster analysis 
	- scott which takes Ne(t) instead of tfgy 
	- scott which takes independent estimate of incidence and prevalence 
- Improve documentation and examples 
- Fast Rcpp version 
