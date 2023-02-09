# age estimation from family structure

# sim

N <- 20 # num families
M <- 3 # num children per woman
IBI <- c(1,4) # IBI range

# age at first birth for each mother
AFB <- 16 + rpois(N,3)

K <- matrix(NA,N,M)
for ( i in 1:N ) {
    # record mom's age when each child born
    K[i,1] <- AFB[i]
    for ( k in 2:M )
        K[i,k] <- K[i,k-1] + runif(1,IBI[1],IBI[2])
}#i

# now sample mom's age and determine kids' ages from that
A <- matrix(NA,N,M+1)
for ( i in 1:N ) {
    delta <- runif(1,5,40)
    A[i,1] <- K[i,M] + delta
    for ( k in 1:M ) A[i,k+1] <- A[i,1] - K[i,k]
}

# K holds mom's ages at each birth
# A holds true ages of moms [,1] and kids [,2:M] in birth order

# now observation error
# error proportional to age
Aobs <- round( A + rnorm(N*(M+1),0,A/10) )

plot(A,Aobs)
abline(a=0,b=1)

dat <- list( N=N , M=M , A=Aobs )

m <- cstan( file="17_measurement_error.stan" , data=dat , chains=1 )

Aobs[,1] - Aobs[,2]
precis(m,2)