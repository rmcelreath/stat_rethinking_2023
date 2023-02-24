# homeweek week 6 optional problem
# simulate bangladesh fertility in agent-based way

sim_fertility <- function(
    n_districts = 60, # number of districts
    p_urban = runif(n_districts,0.1,0.2), # init prob woman moves to city
    n_init = rep(100,n_districts), # init number of women in each district
    f = c(rep(0,20),rep(0.5,20),rep(0,100)), # fertility schedule by age
    m = c(0.2,rep(0.01,19),seq(from=0.01,to=0.5,len=100)), # mortality schedule
    # causal effects
    bKC = rep(0.01,20), # prob adopt C at each parity starting at zero
    bDC = runif(n_districts,0,0.1), # district offsets for prob adopt C
    bUC = 0.1, # effect of urban on prob use C at all parities
    # sim controls
    t_max = 1e2, # number of years to simulate population
    n_max = 1e4 # maximum population size across all districts
) {
    # init population
    # all women start at age 1
    # so need to run sim until reach stable age distribution
    n <- sum(n_init)
    age <- rep(NA,n_max)
    age[1:n] <- 1 # newborns
    D <- rep(NA,n_max) # districts
    D[1:n] <- rep(1:n_districts,times=n_init)
    K <- rep(NA,n_max)
    K[1:n] <- 0 # no kids yet
    U <- rep(NA,n_max)
    U[1:n] <- rbern(n,p_urban[D[1:n]])
    C <- rep(NA,n_max)
    C[1:n] <- rep(0,n) # no one starts using contraception

    # sim loop
    for ( i in 1:t_max ) {
        # loop over living women
        n_births <- 0
        for ( j in 1:n_max ) {
            if ( !is.na(age[j]) ) {
                # she's alive!
                # survive to next year?
                if ( runif(1) > m[age[j]] ) {
                    age[j] <- age[j] + 1 # get older
                    # adopt C?
                    if ( C[j]==0 ) {
                        pC <- bKC[K[j]+1] + bDC[D[j]] + bUC*U[j]
                        if ( pC > 1 ) pC <- 1
                        C[j] <- rbern(1,pC)
                    }
                    if ( C[j] ==0 ) {
                        # birth?
                        if ( runif(1) < f[age[j]] ) {
                            K[j] <- K[j] + 1
                            # is female?
                            if ( runif(1) < 0.5 ) {
                                # add to tally to update at end of this loop
                                n_births <- n_births + 1
                            }
                        }
                    }
                } else {
                    # death - remove from population
                    age[j] <- NA
                }
            }
        }#j

        # now add new women to population
        open_spots <- which(is.na(age))
        n_births <- min( n_births , length(open_spots) ) #bound population
        if ( n_births > 0 )
            for ( j in 1:n_births ) {
                k <- open_spots[j]
                age[k] <- 1
                C[k] <- 0
                D[k] <- sample(60,size=1)
                U[k] <- rbern(1, p_urban[D[k]] )
                K[k] <- 0
            }#j
    }#i

    out <- data.frame(
        A = age, U = U, C = C, K = K, D = D
    )

    # remove empty slots in population
    x <- which(!is.na(out$A))
    out <- out[x,]

    return(out)

}

sim_dat <- sim_fertility()

plot( sim_dat$A , sim_dat$K , xlab="age" , ylab="children" , lwd=3 , col=2 )

dat <- list(
    C = sim_dat$C,
    D = sim_dat$D,
    nD = max(sim_dat$D),
    U = sim_dat$U
)

m0 <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] + b[D]*U,
        # define effects using other parameters
        # this is the non-centered Cholesky machine
        transpars> vector[nD]:a <<- abar[1] + v[,1],
        transpars> vector[nD]:b <<- abar[2] + v[,2],
        transpars> matrix[nD,2]:v <-
            compose_noncentered( sigma , L_Rho , Z ),
        # priors - note that none have parameters inside them
        # that is what makes them non-centered
        matrix[2,nD]:Z ~ normal( 0 , 1 ),
        vector[2]:abar ~ normal(0,1),
        cholesky_factor_corr[2]:L_Rho ~ lkj_corr_cholesky( 4 ),
        vector[2]:sigma ~ exponential(1),
        # convert Cholesky to Corr matrix
        gq> matrix[2,2]:Rho <<- Chol_to_Corr(L_Rho)
    ) , data=dat , chains=4 , cores=4 )

