simulation <- function(n, p, iterate, d, b, sigma)
{
	library(MASS)
	set.seed(12345)
	eps <- .Machine$double.eps

	# build parameter
	coef <- matrix(0, nrow = iterate, ncol = p)
	Cp <- RSS <- matrix(0, nrow = iterate,ncol = 1)
	num <- seq(p)
	model <- list()
	distinmodel <- list()
	mean <- rep(0,p)

	# compute matrix coef,model,Cp,RSS
	t <- 0
	while(t < iterate)
	{
		t <- t + 1
		#e <- 0
		e <- matrix(rnorm(n,mean = 0,sd = d), ncol = 1)
		x <- mvrnorm(n,mean,sigma)
	    # z0 <- matrix(rnorm(p,mean = 0,sd = 1), nrow = 1)
	    # x <- z + rep(1,n) %*% z0
		y <- x %*% b + e
		f <- lasso(x,y,trace = FALSE)
		
		# choose variable
		if(f$CpdisSign)
		{
			lCp <- length(f$Cp)
			rCp <- range(f$Cp)
			newCp <- f$Cp * (lCp/(rCp[2] - rCp[1]))
			point <- cbind(seq(lCp) - 1, newCp)
			distance <- apply(point^2, 1, sum)
			logitemp <- distance == min(distance)
		} else 
		{ logitemp <- f$Cp == min(f$Cp[-1]) }
				
		coef[t, ] <- f$beta[logitemp, ]
		nonzero <- coef[t, ] > eps
		model[[t]] <- num[nonzero]
		Cp[t,] <- f$Cp[logitemp]
		RSS[t,] <- f$RSS[logitemp]
	}

	# compute model accurate proportion   distinmodel, prop
	distinmodel <- unique(model)
	l <- length(distinmodel)
	prop <- no.of0 <- rep(0,l)
	for (i in 1:l)
	{
		li <- length(distinmodel[[i]])
		logitemp <- rep(0,iterate)
		for(j in i:iterate)
		{
			lj <- length(model[[j]])
			if(lj == li)
			{
				if(sum(model[[j]] == distinmodel[[i]]) == li)
				{
					logitemp[j] <- 1
				}			
			}
		}
		prop[i] <- sum(logitemp)/iterate
		no.of0[i] <- p - li
	}
	distinprop <- prop
	q <- seq(l)
	distincoor <- q[prop == max(prop)] # distinmodel
	maxpropmodel <- distinmodel[distincoor]
	aver0 <- no.of0 %*% prop

	y1 <- y2 <- y3 <- y4 <- 0
	topprop <- list()
	for(i in 1:iterate)
	{
		temp <- length(model[[i]])
		if(temp == 0) y1 <- y1 + 1
		else if(temp == 9) y2 <- y2 + 1
		else if(temp == 8) y3 <- y3 + 1
		else if(temp == 7) y4 <- y4 + 1
	}
	topprop[[1]] <- y1/iterate
	topprop[[2]] <- y2/iterate
	topprop[[3]] <- y3/iterate
	topprop[[4]] <- y4/iterate
	object <- list(aver0 = aver0, coef = coef, model = model, 
				   Cp = Cp, RSS = RSS, distinmodel = distinmodel, 
				   distinprop = prop, maxpropmodel = maxpropmodel, 
				   topprop = topprop, no.of0 = no.of0)
	class(object) <- "simulation"
	object
}
