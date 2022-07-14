data {
  int N;
  int N_samples;
  vector[N] y_test;
  vector[N] x_test;
  vector[N_samples] alpha;
  vector[N_samples] beta;
  vector[N_samples] sigma;
}
parameters {
}
model {
}
generated quantities {
  matrix[N_samples,N] log_lik;
  for(n in 1:N){
     for(i in 1:N_samples) {
      log_lik[i,n] = normal_lpdf(y_test[n]|alpha[i] + beta[i]*x_test[n], sigma[i]);
    }  
  }
}