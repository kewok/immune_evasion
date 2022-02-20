library(ssar) # courtesy of these wonderful people: https://github.com/INSP-RH/ssar/

# Unfortunately the code to generate the Gillespie algorithm implementation of the SEIR is not available, so we use a custome SEIR monte-carlo simulation "under the hood" for this step where the fully vaccinated are subtracted from the susceptibles in each step and added to the recovereds.

# The state change/propensity matrix corresponding to 
# A <- c("b*S0*I0/n","u(t)*S0","a*E0","g*I0")
# should look like:

# S -1 -1 0 0
# E +1 0 -1 0 
# I 0 0  +1 -1
# R 0 0 0  +1
# V 0 +1 0 0

v <- matrix(c(
 -1, -1, 0, 0,
 +1, 0, -1, 0, 
 0, 0,  +1, -1,
 0, +1, 0, +1,
 0, +1, 0, 0
		),nrow=5, byrow=TRUE)

seir_custom <- function(N, I0, E0, R0, V0, b, a, g, Time, runs, vaxRate, complete_trajectory = FALSE)
	{
	X <- matrix(c(S0=N-I0[1]-E0[1]-R0[1]-V0, E0[1], I0[1], R0[1], V0), nrow=1)
	parameters <- c('b'=as.numeric(b), 'n'=as.numeric(N), 'a'=as.numeric(a), 'g'=as.numeric(g))

	u <- function(t) vaxRate[round(t+1)]
	# The time-dependent propensity function now becomes:

	pfun <- function(t, X, parameters)
		{
		cbind(parameters['b'] * X[,1] * X[,3]/parameters['n'], u(t) * X[,1], parameters['a']*X[,2], parameters['g']*X[,3])
		}

	res <- ssa(X, pfun, v, parameters, tmin=0, tmax=Time, nsim=runs, plot.sim=FALSE)
	if (complete_trajectory)
		return(res)
	else
		return(res[(nrow(res)-runs + 1):nrow(res),])
	}

