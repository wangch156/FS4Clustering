pmf_zero_inflated_poisson <- function(x, lambda, theta){
  index_0 = which(x == 0)
  ret = rep(0, length(x))
  ret[index_0] = theta + (1-theta)*dpois(0, lambda)
  ret[-index_0] = (1-theta)*dpois(x[-index_0], lambda)
  # if (x == 0){
  #   return(theta + (1-theta)*dpois(x, lambda))
  # } else {
  #   return((1-theta)*dpois(x, lambda))
  # }
  return(ret)
}
estimate_zi_poisson <- function(data, weights) {
  index_zero = which(data == 0)
  r0 <- sum(weights[index_zero]) / sum(weights)

  x.bar = sum(data * weights) / sum(weights)

  gamma <- x.bar / (1 - r0)

  lambda.hat <- lamW::lambertW0(-gamma * exp(-gamma)) + gamma

  pi.hat <- 1 - x.bar / lambda.hat


  return <- c(lambda.hat, pi.hat)
  return(return)
}

estimate_zi_poisson_homo <- function(data) {
  num.zeros <- sum(data == 0)
  r0 <- 1 / length(data) * num.zeros

  x.bar = mean(data)

  gamma <- x.bar / (1 - r0)

  lambda.hat <- lamW::lambertW0(-gamma * exp(-gamma)) + gamma

  pi.hat <- 1 - x.bar / lambda.hat


  return <- c(lambda.hat, pi.hat)
  return(return)
}

em.stat.zero_inflated_poisson <- function(x,alpha.ini,k0,C,labels,group.num){ #,prior.weight = 0.05,earlystop = NA, is_ecdf = F){
  if (length(unique(labels))!=group.num){
    print("The number of different labels does not match the total number of groups!")
    return(NULL)
  }
  
  M.stat <- matrix(nrow = 1+k0,ncol = nrow(alpha.ini))
  # em_zero_inflated_poisson <- function(x, lambda_vec, pi_vec, theta_vec, max_iter = 10){
  n = length(x)
  k = group.num
  max_iter = k0
# em_zero_inflated_poisson <- function(x, lambda_vec, pi_vec, theta_vec, max_iter = 10){

  mu <- rep(0,group.num)
  for (i in 1:group.num){
    mu[i] <- mean(x[labels==i])
  }
  w0 = matrix(0.001,nrow = length(x),ncol = group.num)
  for (i in 1:group.num){
    w0[labels==i,i] <- 1-0.001*(group.num-1)
  }
  # theta0 <- rep(sum(x==0)/length(x),group.num)
  theta0 <- rep(0.5,group.num)
  for (ll in 1:nrow(alpha.ini)){
    pi_vec = c(alpha.ini[ll,],1-sum(alpha.ini[ll,]))
    lambda_vec = mu
    theta_vec = theta0
    w = w0
    for (i in 1:max_iter){
      for (j in 1:k){
        w[,j] = pi_vec[j] * pmf_zero_inflated_poisson(x, lambda_vec[j], theta_vec[j])
      }
      w = w / rowSums(w)
      pi_vec = colSums(w) / n
      for (j in 1:k){
        est = estimate_zi_poisson(x, w[,j])
        lambda_vec[j] = est[1]
        theta_vec[j] = est[2]
      }
      w_old = w
      likelihood = rep(0, n)
      for (j in 1:k){
        likelihood = likelihood + (pi_vec[j] * (pmf_zero_inflated_poisson(x, lambda_vec[j], theta_vec[j])))
      }
      M.stat[i+1,ll] = sum(log(likelihood))
    }
  }
  return(c("pvalue" = 1-pchisq(max(M.stat, na.rm = T),3),"stat" = max(M.stat, na.rm = T)))
}

likelihood_zero_inflated_poisson <- function(x){
  n = length(x)
  para = estimate_zi_poisson_homo(x)
  ret = pmf_zero_inflated_poisson(x, para[1], para[2])
  return(sum(log(ret)))
}
