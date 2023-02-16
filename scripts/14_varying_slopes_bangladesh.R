# week 7
# varying effects, clusters and features, non-centering

library(rethinking)

# simple varying intercepts model
library(rethinking)
data(bangladesh)
d <- bangladesh


# visualize a and b distributions as independent gaussian

blank2()

library(ellipse)
rho <- 0.5
SIGMA <- matrix( c(1,rho,rho,1) , 2 , 2 )
MU <- c( 0 , 0 )

# SIGMA <- rlkjcorr(1,2,eta=4)

plot( NULL , xlim=c(-3,3) , ylim=c(-3,3) , xlab="a" , ylab="b" )
for ( l in seq(from=0.25,to=0.95,len=5) ) {
    el <- ellipse( SIGMA , centre=MU , level=l )
    #lines( (el) , col=2 , lwd=3 )
    polygon( (el) , col=col.alpha(2,0.25) , border=NA )
}

Y <- rmvnorm(6,c(0,0),sigma=SIGMA)
points( Y , lwd=4 , col="white" )
points( Y , lwd=3 , col=2 )

# lkjcorr prior predictive

plot( NULL , xlim=c(-2.5,2.5) , ylim=c(-2.5,2.5) , xlab="a" , ylab="b" )
for ( i in 1:10 ) {
    RHO <- rlkjcorr(1,2,eta=4)
    s <- rexp(1,1)
    tau <- rexp(1,1)
    SIGMA <- diag(c(s,tau)) %*% RHO %*% diag(c(s,tau))
    el <- ellipse( SIGMA , centre=MU , level=0.89 )
    lines( (el) , col=col.alpha(2,0.5) , lwd=3 )
    #polygon( (el) , col=col.alpha(2,0.25) , border=NA )
}

SIGMA <- rlkjcorr(1e4,2,eta=4)
dens(SIGMA[,1,2],lwd=3,col=2,xlab="correlation")

###########
# non-centered varying slopes with and without covariance

dat <- list(
    C = d$use.contraception,
    D = as.integer(d$district),
    U = d$urban,
    A = standardize(d$age.centered),
    K = d$living.children )

# no covariance
mCDUnc <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] + b[D]*U,
        # define effects using other parameters
        save> vector[61]:a <<- abar + za*sigma,
        save> vector[61]:b <<- bbar + zb*tau,
        # z-scored effects
        vector[61]:za ~ normal(0,1),
        vector[61]:zb ~ normal(0,1),
        # ye olde hyper-priors
        c(abar,bbar) ~ normal(0,1),
        c(sigma,tau) ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )

# covariance - centered
mCDUcov <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] + b[D]*U,
        # define effects using other parameters
        transpars> vector[61]:a <<- v[,1],
        transpars> vector[61]:b <<- v[,2],
        # priors - centered correlated varying effects
        matrix[61,2]:v ~ multi_normal(abar,Rho,sigma),
        vector[2]:abar ~ normal(0,1),
        corr_matrix[2]:Rho ~ lkj_corr(4),
        vector[2]:sigma ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )

# covariance - non-centered
mCDUcov_nc <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] + b[D]*U,
        # define effects using other parameters
        # this is the non-centered Cholesky machine
        transpars> vector[61]:a <<- abar[1] + v[,1],
        transpars> vector[61]:b <<- abar[2] + v[,2],
        transpars> matrix[61,2]:v <-
            compose_noncentered( sigma , L_Rho , Z ),
        # priors - note that none have parameters inside them
        # that is what makes them non-centered
        matrix[2,61]:Z ~ normal( 0 , 1 ),
        vector[2]:abar ~ normal(0,1),
        cholesky_factor_corr[2]:L_Rho ~ lkj_corr_cholesky( 4 ),
        vector[2]:sigma ~ exponential(1),
        # convert Cholesky to Corr matrix
        gq> matrix[2,2]:Rho <<- Chol_to_Corr(L_Rho)
    ) , data=dat , chains=4 , cores=4 )

precis(mCDUcov_nc,3,pars=c("Rho","sigma"))
precis(mCDUnc,3,pars=c("sigma","tau"))

# posterior rho
post <- extract.samples(mCDUcov_nc)
dens( post$Rho[,1,2] , xlim=c(-1,1) , lwd=3 , col=2 , xlab="posterior correlation a,b" )
abline(v=0,lty=2,lwd=0.5)
prior_rho <- rlkjcorr(1e4,2,eta=4)
dens( prior_rho[,1,2] , lwd=2 , lty=2 , add=TRUE )

# posterior MVN of a,b
plot( NULL , xlim=c(-2,1) , ylim=c(-1,2) , xlab="a" , ylab="b" )
abline(v=0,lty=2,lwd=0.5)
abline(h=0,lty=2,lwd=0.5)
SIGMA <- cov(cbind( apply(post$a,2,mean) , apply(post$b,2,mean) ) )
MU <- apply(post$abar,2,mean)
for ( l in seq(from=0.25,to=0.95,len=5) ) {
    el <- ellipse( SIGMA , centre=MU , level=l )
    #lines( (el) , col=2 , lwd=3 )
    polygon( (el) , col=col.alpha(2,0.25) , border=NA )
}

#points( apply(post$a,2,mean) , apply(post$b,2,mean) , col="white", lwd=3 )
points( apply(post$a,2,mean) , apply(post$b,2,mean) , col=1, lwd=2 )

post2 <- extract.samples(mCDUnc)
points( apply(post2$a,2,mean) , apply(post2$b,2,mean) , col=1, lwd=2 )

# plot estimates

Uval <- 1
xcol <- ifelse(Uval==0,2,4)
p2 <- link( mCDUnc , data=list(D=1:61,U=rep(Uval,61)) )
#p2 <- link( mCDUcov_nc , data=list(D=1:61,U=rep(Uval,61)) )

# blank2(w=2,h=0.8)
plot( NULL , xlab="district" , lwd=3 , col=2 , xlim=c(1,61), ylim=c(0,1) , ylab="prob use contraception" )
abline(h=0.5,lty=2,lwd=0.5)

#points( 1:61 , apply(p,2,mean) , xlab="district" , lwd=3 , col=grau(0.8) , ylim=c(0,1) )

points( 1:61 , apply(p2,2,mean) , xlab="district" , lwd=3 , col=xcol , ylim=c(0,1) )

#for ( i in 1:61 ) lines( c(i,i) , PI(p2[,i]) , lwd=8 , col=col.alpha(xcol,0.5) )

# show other feature
Uvalx <- 1-Uval
xcolx <- ifelse(Uvalx==0,2,4)
p2x <- link( mCDUcov_nc , data=list(D=1:61,U=rep(Uvalx,61)) )
points( 1:61 , apply(p2x,2,mean) , lwd=3 , col=xcolx )

# show raw proportions - have to skip 54
n <- table(dat$D,dat$U)
Cn <- xtabs(dat$C ~ dat$D + dat$U)
pC <- as.numeric( Cn[,Uval+1]/n[,Uval+1] )
pC <- c( pC[1:53] , NA , pC[54:60] )
#points( pC , lwd=2 )

# only some labels via locator
nn <- as.numeric(n[,Uval+1])
nn <- c( nn[1:53] , 0 , nn[54:60] )
#identify( 1:61 , pC , labels=nn , cex=1 )




# shrinkage plot now
# blank2(w=1)

idx <- 34
idx <- 1:61


plot( NULL , xlab="prob C (rural)" , ylab="prob C (urban)" , xlim=c(0,1), ylim=c(0,1) )

plot( NULL , xlab="prob C (rural)" , ylab="prob C (urban)" , xlim=c(0.1,0.65), ylim=c(0.2,0.75) )

abline(h=0.5,lty=2,lwd=0.5)
abline(v=0.5,lty=2,lwd=0.5)

# point sizes proportional to smaple size in district
n <- table(dat$D)
n <- c( n[1:53] , 0 , n[54:60] )

# uncorrelated model
post <- extract.samples(mCDUnc)
logitp0 <- post$a
logitp1 <- post$a + post$b
p0 <- inv_logit(logitp0)
p1 <- inv_logit(logitp1)
#points( apply(p0,2,mean) , apply(p1,2,mean) , lwd=6 , col="white" )
points( apply(p0,2,mean)[idx] , apply(p1,2,mean)[idx] , lwd=2 , col=1 )

# correlated model
post <- extract.samples(mCDUcov_nc)
logitp0 <- post$a
logitp1 <- post$a + post$b
p0 <- inv_logit(logitp0)
p1 <- inv_logit(logitp1)
points( apply(p0,2,mean)[idx] , apply(p1,2,mean)[idx] , lwd=5 , col="white" )
points( apply(p0,2,mean)[idx] , apply(p1,2,mean)[idx] , lwd=3 , col=2 , cex=1 )



n <- table(dat$D,dat$U)
Cn <- xtabs(dat$C ~ dat$D + dat$U)
pC0 <- as.numeric( Cn[,1]/n[,1] )
pC1 <- as.numeric( Cn[,2]/n[,2] )

points( (pC0)[idx] , (pC1)[idx] , lwd=2 , cex=1 , pch=16 )
#points( (pC0) , (pC1) , lwd=2 , cex=2*n[,1]/100 + 0.5 )
#points( (pC0) , (pC1) , lwd=2 , cex=2*n[,2]/100 + 0.5 , col=4 )

p0x <- apply(p0,2,mean)
p1x <- apply(p1,2,mean)
for ( i in 1:61 ) {
    lines( c(pC0[i], p0x[i] ) , c(pC1[i], p1x[i] ) , col=grau() )
}


# show raw proportions - have to skip 54
n <- table(dat$D)
Cn <- xtabs(dat$C ~ dat$D)
pC <- as.numeric( Cn/n )
pC <- c( pC[1:53] , NA , pC[54:60] )
points( pC , lwd=2 )

# only some labels via locator
n <- table(dat$D)
n <- as.numeric(n)
n <- c( n[1:53] , 0 , n[54:60] )
identify( 1:61 , pC , labels=n , cex=1 )



####
# simple Cholesky factor example

R <- matrix(c(1,0.6,0.6,1),2,2)
L <- chol(R)
