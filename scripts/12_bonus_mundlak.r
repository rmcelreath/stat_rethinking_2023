# endogenous group confound example

set.seed(8672)

N_groups <- 30
N_id <- 200
a0 <- (-2)
bZY <- (-0.5)
g <- sample(1:N_groups,size=N_id,replace=TRUE) # sample into groups
Ug <- rnorm(N_groups,1.5) # group confounds
X <- rnorm(N_id, Ug[g] ) # individual varying trait
Z <- rnorm(N_groups) # group varying trait (observed)
Y <- rbern(N_id, p=inv_logit( a0 + X + Ug[g] + bZY*Z[g] ) )

table(g)

# confounded by correlation
precis(glm(Y~X+Z[g],family=binomial),2)

# fixed effects
# X deconfounded, but Z unidentified now!
precis(glm(Y~X+Z[g]+as.factor(g),family=binomial),pars=c("X","Z"),2)

dat <- list(Y=Y,X=X,g=g,Ng=N_groups,Z=Z)

# naive model
m0 <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a + bxy*X + bzy*Z[g],
        a ~ dnorm(0,10),
        c(bxy,bzy) ~ dnorm(0,1)
    ) , data=dat , chains=4 , cores=4 )

# fixed effects
mf <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bxy*X + bzy*Z[g],
        a[g] ~ dnorm(0,10),
        c(bxy,bzy) ~ dnorm(0,1)
    ) , data=dat , chains=4 , cores=4 )

# random effects
mr <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bxy*X + bzy*Z[g],
        transpars> vector[Ng]:a <<- abar + z*tau,
        z[g] ~ dnorm(0,1),
        c(bxy,bzy) ~ dnorm(0,1),
        abar ~ dnorm(0,1),
        tau ~ dexp(1)
    ) , data=dat , chains=4 , cores=4 , sample=TRUE )

# random effects + Xbar
# The Mundlak Machine
xbar <- sapply( 1:N_groups , function(j) mean(X[g==j]) )
dat$Xbar <- xbar
mrx <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bxy*X + bzy*Z[g] + buy*Xbar[g],
        transpars> vector[Ng]:a <<- abar + z*tau,
        z[g] ~ dnorm(0,1),
        c(bxy,buy,bzy) ~ dnorm(0,1),
        abar ~ dnorm(0,1),
        tau ~ dexp(1)
    ) , data=dat , chains=4 , cores=4 , sample=TRUE )

# random effects + latent U
# The Latent Mundlak Machine
mru <- ulam(
    alist(
        # Y model
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bxy*X + bzy*Z[g] + buy*u[g],
        transpars> vector[Ng]:a <<- abar + z*tau,
        # X model
        X ~ normal(mu,sigma),
        mu <- aX + bux*u[g],
        vector[Ng]:u ~ normal(0,1),
        # priors
        z[g] ~ dnorm(0,1),
        c(aX,bxy,buy,bzy) ~ dnorm(0,1),
        bux ~ dexp(1),
        abar ~ dnorm(0,1),
        tau ~ dexp(1),
        sigma ~ dexp(1)
    ) , data=dat , chains=4 , cores=4 , sample=TRUE )

precis(mf)
precis(mr)
precis(mrx)
precis(mru)

# density plots
# bxy
post <- extract.samples(mf)
dens(post$bxy,lwd=3,col=1,xlab="b_XY",ylim=c(0,2))
abline(v=1,lty=2)

post <- extract.samples(m0)
dens(post$bxy,lwd=3,col=grau(),add=TRUE)

post <- extract.samples(mr)
dens(post$bxy,lwd=3,col=2,add=TRUE)

post <- extract.samples(mrx)
dens(post$bxy,lwd=3,col=4,add=TRUE)

post <- extract.samples(mru)
dens(post$bxy,lwd=8,col="white",add=TRUE)
dens(post$bxy,lwd=4,col=3,add=TRUE)


# bzy
post <- extract.samples(mf)
dens(post$bzy,lwd=3,col=1,xlab="b_ZY",ylim=c(0,2))
abline(v=-0.5,lty=2)

post <- extract.samples(m0)
dens(post$bzy,lwd=3,col=grau(),add=TRUE)

post <- extract.samples(mr)
dens(post$bzy,lwd=3,col=2,add=TRUE)

post <- extract.samples(mrx)
dens(post$bzy,lwd=3,col=4,add=TRUE)

post <- extract.samples(mru)
dens(post$bzy,lwd=8,col="white",add=TRUE)
dens(post$bzy,lwd=4,col=3,add=TRUE)


##########
# show better estimates of intercepts

af <- coef(mf)[1:N_groups]
ar <- coef(mr)[1:N_groups]

plot( af , col=4 )
points( 1:N_groups , ar , col=2 )
points( 1:N_groups , a0+Ug , col=1 )


# treatment effect in each group now
# counterfactual increase of X at individual level, stratified by each group

# fixed estimates
pf0 <- link(mf,data=list(g=1:N_groups,X=rep(0,N_groups)))
pf1 <- link(mf,data=list(g=1:N_groups,X=rep(1,N_groups)))
cf <- apply( pf1 - pf0 , 2 , mean )

# random estimates
pr0 <- link(mr,data=list(g=1:N_groups,X=rep(0,N_groups)))
pr1 <- link(mr,data=list(g=1:N_groups,X=rep(1,N_groups)))
cr <- apply( pr1 - pr0 , 2 , mean )

# true
ctrue <- inv_logit( a0 + Ug + 1 ) - inv_logit( a0 + Ug )

plot( ctrue , ylim=c(0,0.3) )
points( 1:N_groups , cf , col=4 )
points( 1:N_groups , cr , col=2 )

plot( cf - ctrue , col=4 , ylim=c(-0.1,0.1) )
points( 1:N_groups , cr - ctrue , col=2 )
abline(h=0,lty=2)

mean((cf - ctrue)^2)
mean((cr - ctrue)^2)
