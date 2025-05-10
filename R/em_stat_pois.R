# library(iZID)
# library(lamW)

pmf_zero_inflated_poisson2 <- function(x, lambda, theta = 0){
  index_0 = which(x == 0)
  ret = rep(0, length(x))
  ret[index_0] = theta + (1-theta)*dpois(0, lambda)
  ret[-index_0] = (1-theta)*dpois(x[-index_0], lambda)
  return(ret)
}

estimate_zi_poisson_homo2 <- function(data) {
  x.bar = mean(data)
  return <- x.bar
  return(return)
}

em.stat.poisson <- function(x,alpha.ini,k0,C,labels,group.num){ #,prior.weight = 0.05,earlystop = NA, is_ecdf = F){
  if (length(unique(labels))!=group.num){
    print("The number of different labels does not match the total number of groups!")
    return(NULL)
  }
  
  M.stat <- matrix(nrow = 1+k0,ncol = nrow(alpha.ini))
  # em_zero_inflated_poisson <- function(x, lambda_vec, pi_vec, theta_vec, max_iter = 10){
  n = length(x)
  k = group.num
  max_iter = k0
  
  mu <- rep(0,group.num)
  for (i in 1:group.num){
    mu[i] <- mean(x[labels==i])
  }
  w0 = matrix(0.001,nrow = length(x),ncol = group.num)
  for (i in 1:group.num){
    w0[labels==i,i] <- 1-0.001*(group.num-1)
  }
  
  for (ll in 1:nrow(alpha.ini)){
    pi_vec = c(alpha.ini[ll,],1-sum(alpha.ini[ll,]))
    lambda_vec = mu
    w = w0
    w_old = w0
    for (i in 1:max_iter){
      for (j in 1:k){
        w[,j] = pi_vec[j] * pmf_zero_inflated_poisson2(x, lambda_vec[j])
      }
      w = w / rowSums(w)
      pi_vec = colSums(w) / n
      lambda_vec = colSums(w * x) / colSums(w)
      w_old = w
      likelihood = rep(0, n)
      for (j in 1:k){
        likelihood = likelihood + (pi_vec[j] * (pmf_zero_inflated_poisson2(x, lambda_vec[j])))
      }
      M.stat[i+1,ll] = sum(log(likelihood))
    }
  }
  return(c("pvalue" = 1-pchisq(max(M.stat, na.rm = T),1),"stat" = max(M.stat, na.rm = T)))
}

likelihood_zero_inflated_poisson <- function(x){
  n = length(x)
  para = mean(x)
  ret = pmf_zero_inflated_poisson(x, para)
  return(sum(log(ret)))
}

rzipoisson <- function(lambda, n) {
  data <- rpois(n, lambda)
  return(data)
} 