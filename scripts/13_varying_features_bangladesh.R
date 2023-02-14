# week 7
# varying effects, clusters and features, non-centering

library(rethinking)

# simple varying intercepts model
library(rethinking)
data(bangladesh)
d <- bangladesh

dat <- list(
    C = d$use.contraception,
    D = as.integer(d$district) )

mCD <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D],
        vector[61]:a ~ normal(abar,sigma),
        abar ~ normal(0,1),
        sigma ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )


# plot estimates
p <- link( mCD , data=list(D=1:61) )
# blank2(w=2)
plot( NULL , xlab="district" , lwd=3 , col=2 , xlim=c(1,61), ylim=c(0,1) , ylab="prob use contraception" )

points( 1:61 , apply(p,2,mean) , xlab="district" , lwd=3 , col=2 , ylim=c(0,1) , ylab="prob use contraception" )

 for ( i in 1:61 ) lines( c(i,i) , PI(p[,i]) , lwd=8 , col=col.alpha(2,0.5) )

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




#####################
# add urban category

dat <- list(
    C = d$use.contraception,
    D = as.integer(d$district),
    U = ifelse(d$urban==1,1,0) )

# total U
mCDU <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] + b[D]*U,
        vector[61]:a ~ normal(abar,sigma),
        vector[61]:b ~ normal(bbar,tau),
        c(abar,bbar) ~ normal(0,1),
        c(sigma,tau) ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )

traceplot(mCDU,pars="tau",lwd=2,n_cols=1)
trankplot(mCDU,pars="tau",lwd=3,n_cols=1)

# non-centered version
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

# plot estimates

Uval <- 0
xcol <- ifelse(Uval==0,2,4)
p <- link( mCDUnc , data=list(D=1:61,U=rep(Uval,61)) )
# blank2(w=2,h=0.8)
plot( NULL , xlab="district" , lwd=3 , col=2 , xlim=c(1,61), ylim=c(0,1) , ylab="prob use contraception" )
abline(h=0.5,lty=2,lwd=0.5)

points( 1:61 , apply(p,2,mean) , xlab="district" , lwd=3 , col=xcol , ylim=c(0,1) , ylab="prob use contraception" )

 for ( i in 1:61 ) lines( c(i,i) , PI(p[,i]) , lwd=8 , col=col.alpha(xcol,0.5) )

# show raw proportions - have to skip 54
n <- table(dat$D,dat$U)
Cn <- xtabs(dat$C ~ dat$D + dat$U)
pC <- as.numeric( Cn[,Uval+1]/n[,Uval+1] )
pC <- c( pC[1:53] , NA , pC[54:60] )
points( pC , lwd=2 )

# only some labels via locator
nn <- as.numeric(n[,Uval+1])
nn <- c( nn[1:53] , 0 , nn[54:60] )
identify( 1:61 , pC , labels=nn , cex=1 )

# show standard deviations
post <- extract.samples(mCDUnc)
dens(post$sigma,xlab="posterior standard deviation",lwd=3,col=2,xlim=c(0,1.2))
dens(post$tau,lwd=3,col=4,add=TRUE,adj=0.2)
curve(dexp(x,1),from=0,to=1.3,add=TRUE,lwd=2,lty=2)

####
# shrinkage plot now
post <- extract.samples(mCDUnc)
logitp0 <- post$a
logitp1 <- post$a + post$b

# blank2(w=1)
#plot( NULL , xlab="log-odds C (rural)" , ylab="log-odds C (urban)" , xlim=c(-2,1), ylim=c(-1.5,1.5) )

plot( NULL , xlab="prob C (rural)" , ylab="prob C (urban)" , xlim=c(0.1,0.7), ylim=c(0.2,0.75) )
abline(h=0.5,lty=2,lwd=0.5)
abline(v=0.5,lty=2,lwd=0.5)

# plausibility ellipses
library(ellipse)
xxx <- sample(1:61,size=6)
for ( i in xxx ) {
    SIGMA <- cov( cbind( logitp0[,i] , logitp1[,i] ) )
    MU <- c( mean(logitp0[,i]) , mean(logitp1[,i]) )
    el <- ellipse( SIGMA , centre=MU , level=0.5 )
    lines( inv_logit(el) , col=col.alpha(2,0.3) , lwd=2 )
    #polygon( inv_logit(el) , col=col.alpha(2,0.2) , border=NA )
}

# posterior means
p0 <- inv_logit(logitp0)
p1 <- inv_logit(logitp1)
points( apply(p0,2,mean) , apply(p1,2,mean) , lwd=6 , col="white" )
points( apply(p0,2,mean) , apply(p1,2,mean) , lwd=3 , col=2 )

n <- table(dat$D,dat$U)
Cn <- xtabs(dat$C ~ dat$D + dat$U)
pC0 <- as.numeric( Cn[,1]/n[,1] )
pC1 <- as.numeric( Cn[,2]/n[,2] )

points( (pC0) , (pC1) , lwd=2 , cex=2*apply(n,1,sum)/100 + 0.5 )

for ( i in 1:61 ) {
    lines( c(pC0[i],p0x[i])) , c(pC1[i],p1x[i]) )
}

