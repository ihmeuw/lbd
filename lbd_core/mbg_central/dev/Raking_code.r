require(faraway)

G <- function(x){
	val <- logit(x)
	return(val)
}

G_Inv <- function(x){
	val <- ilogit(x)
	return(val)
}

Agg <- function(p_i, N_i){
	return(sum(p_i * N_i)/sum(N_i))
}


EvalK <- function(K, p_i, p_N, N_i){
	if (Mult == TRUE){
		p_tilde <- G_Inv(G(p_i) * K)
	} else {
		p_tilde <- G_Inv(G(p_i) + K)
	}
	LHS <- Agg(p_tilde, N_i)
	RHS <- p_N
	SE <- (LHS-RHS)^2
	return(SE)
}

# p_i = our grid of probabilities
# p_N = national target
# sim_N_i = our grid of sample sizes (population raster)
# a = bounds for raking factor, perhaps return warning if optimal is on edge of limits 
FindK <- function(p_i, p_N, N_i, a) { 
	if (Mult == TRUE){
		Limits <- c(0,a)
	} else {
		Limits <- c(-a,a)
	}
	iter <- 1
	Boundary = TRUE
	while (Boundary & iter < 10){
		Limits <- Limits * iter
		val <- optimize(EvalK, Limits, p_i = p_i, p_N = p_N, N_i = N_i, tol = 1e-20)$min
		Boundary <- (round(abs(val)) == Limits[2])
		iter <- iter + 1
	}
	return(val)
}

NumLocs <- 100
sim_p_i <- runif(NumLocs)
sim_N_i <- round(runif(NumLocs)+runif(NumLocs)*100*runif(NumLocs))

Mult <- FALSE
a <- 2

Try <- 99

p_Ns <- seq(0, 1, length = Try + 2)[-c(1,Try+2)]
Scale <- seq(-10, 10, length = Try + 2)[-c(1,Try+2)]
AggP <- rep(NA, length=Try)

RF <- matrix(0, Try, Try)

for (i in 1:Try){
	tmp_p_i <- ilogit(logit(sim_p_i) + Scale[i])
	AggP[i] <- Agg(tmp_p_i, sim_N_i)
	for (j in 1:Try){
		tmp_p <- p_Ns[j]
		RF[i,j] <- FindK(tmp_p_i, tmp_p, sim_N_i, a)
	}
}

persp(AggP,p_Ns,RF,phi=35,theta=45)


