// age estimation model
data{
    int N; // num mothers
    int M; // num kids per mother
    matrix[N,M+1] A; // observed ages
}
parameters{
    matrix[N,M+1] Ae;
    real log_sigma;
}
model{
    log_sigma ~ normal(0,0.5);
    to_vector(Ae) ~ normal(0,)

    for ( i in 1:N ) {
        // mom's age
        (A[i,1] - A[i,2]) ~ normal( 
            age_at_first_birth[i] , 
            exp(log_sigma + age_at_first_birth[i]*0.1) );


    }//i
}
