my.normal <- function(theta,x){
  # The first component of theta is mean, while the second is r
  mu = theta[1]
  sigma2 = theta[2]
  return((1/sqrt(sigma2))*exp(-(x-mu)^2/(2*sigma2)))
}

pl <- function(theta,x,alpha,C,group.num){
  # alpha is a vector of length **n-1**, where n is the number of groups.
  # The penalty term should reach its maximum at alpha = (1/n,...,1/n) 
  # and tend to -Inf when one component reaches 0
  theta.tmp <- matrix(theta,nrow = group.num)
  tmp <- apply(theta.tmp, 1, my.normal,x = x)
  alpha <- c(alpha,1-sum(alpha))
  sum(log(tmp %*% alpha))+C*log(prod(alpha))
}

l <- function(theta,x,alpha,C,group.num){
  # alpha is a vector of length **n-1**, where n is the number of groups.
  # The penalty term should reach its maximum at alpha = (1/n,...,1/n) 
  # and tend to -Inf when one component reaches 0
  theta.tmp <- matrix(theta,nrow = group.num)
  tmp <- apply(theta.tmp, 1, my.normal,x = x)
  alpha <- c(alpha,1-sum(alpha))
  return((tmp %*% alpha))
}

pl.homo <- function(theta0,x,alpha,C,group.num){
  return(pl(c(rep(theta0[1],group.num),rep(theta0[2],group.num)),x,alpha,C,group.num))
}

l.homo <- function(theta0,x,alpha,C,group.num){
  return(l(c(rep(theta0[1],group.num),rep(theta0[2],group.num)),x,alpha,C,group.num))
}

# qalpha <- function(alpha,w,C){
#   alpha <- c(alpha,1-sum(alpha))
#   tmp <- sum(apply(w, 2, sum)*log(alpha))+C*log(prod(alpha))
# }

em.stat.normal <- function(x,alpha.ini,k0,C,labels,group.num,prior.weight = 0.05,earlystop = NA, is_ecdf = F){
  # Calculate EM statistics
  # k0 is the number of iterations, and alpha.ini are the initial alpha values
  # C is the coefficient of the penalty term
  # The first n components of `theta` are means, while the last n components are dispersion parameters (r)
  # Labels should be in the form of 1,2,...,n
  # Each row of alpha.ini contains a set of initial values
  if (length(unique(labels))!=group.num){
    print("The number of different labels does not match the total number of groups!")
    return(NULL)
  }
  
  # EM.time <- 0
  # opt.time <- 0
  # weight.time <- 0
  
  M.stat <- matrix(nrow = 1+k0,ncol = nrow(alpha.ini))
  
  # tmp <- proc.time()[3]
  theta.0 <- c(mean(x), var(x))
  pl0 <- pl.homo(theta.0,x,rep(1/group.num,group.num-1),C,group.num)
  mu <- numeric(group.num)
  sigma <- numeric(group.num)
  for (i in 1:group.num){
    mu[i] <- mean(x[labels==i])
    # mu[i] <- 20
    
    ##### hetero-var #####
    # sigma[i] <- var(x[labels==i])
    ##### homo-var ######
    sigma[i] <- var(x)
    
  }
  # print(pl0)
  # mu[mu<=1e-9] <- 1e-9*1.01
  sigma[sigma<=1e-2] <- 1e-2*1.01
  # mu[mu>=1e5] <- 1e5
  sigma[sigma>=1e5] <- 1e5
  # cat("Preparation Time: ",proc.time()[3]-tmp,"\n")
  w0 = matrix(0.001,nrow = length(x),ncol = group.num)
  for (i in 1:group.num){
    w0[labels==i,i] <- 1-0.001*(group.num-1)
  }
  #########
  for (j in 1:nrow(alpha.ini)){
    # Calculate the initial M values
    
    # tmp <- proc.time()[3]
    k <- 0
    alpha.old <- alpha.ini[j,]
    theta.old <- c(mu,sigma)
    pl.value <- pl(theta.old,x,alpha.old,C,group.num)
    # print(pl.value)
    M.stat[k+1,j] <- 2 * (pl.value-pl0)
    # cat("Time for Homogenuous model: ",proc.time()[3]-tmp,"\n")
    for (k in 1:k0){
      #cat("Iter ",k,"\n")
      # tmp <- proc.time()[3]
      theta.new <- theta.old
      w = matrix(0, nrow = length(x) ,ncol = group.num)
      for (i in 1:(group.num-1)){
        w[,i] = (alpha.old[i] * l.homo(c(theta.old[i],theta.old[i+group.num] ),x,alpha.old,C,group.num)) / 
          l(theta.old,x,alpha.old,C,group.num) 
      } 
      w[,group.num] =((1-sum(alpha.old)) * l.homo(c(theta.old[group.num],theta.old[group.num+group.num] ),x,alpha.old,C,group.num)) / 
        l(theta.old,x,alpha.old,C,group.num) 
      # print(w)
      w =  (1-prior.weight)*w +prior.weight*w0
      # var_sum = 0
      for (i in 1:group.num){
        w_i = w[, i]
        mean = sum(w_i * x) / sum(w_i)
        theta.new[i] = mean  
        #### homo-var #######
        # var_sum = var_sum + sum(w_i * (x -mean)^2)
        ##### hetero-var ######
        theta.new[i+group.num] = (sum(w_i * (x -mean)^2)  + 0.5*var(x))/ (sum(w_i)+0.5)
        
        # if ((sum(w_i * (x -mean)^2) / sum(w_i)) < 1e-2) {
        #   theta.new[i+group.num] = 1e-2
        # } else {
        #   theta.new[i+group.num] = sum(w_i * (x -mean)^2) / sum(w_i)
        # }
      }
      #### homo-var #######
      
      # var_homo = (var_sum +  0.25*var(x)) / (n+0.25)
      # # # print(theta.new)
      # theta.new[(group.num+1):(2*group.num)] = rep(var_homo, group.num)
      # print(theta.new)
      #### homo-var #######
      
      # EM.time <- EM.time+proc.time()[3]-tmp
      # theta.new[theta.new<=1e-10] <- 1e-10
      # theta.new[theta.new>=1e5] <- 1e5
      
      # tmp <- proc.time()[3]
      alpha.new = (colSums(w) + C) / (length(x) + group.num * C)
      alpha.new = alpha.new[1:(group.num - 1)]
      # opt.time <- opt.time+proc.time()[3]-tmp
      
      M.stat[k+1,j] <- 2 * (pl(theta.new,x,alpha.new,C,group.num)-pl0)
      
      ## Early Stopping
      if (!is.na(earlystop)){
        if (abs(M.stat[k+1,j]-M.stat[k,j])<=earlystop){
          #print("Early Stopped!")
          M.stat[is.na(M.stat[,j]),j] <- M.stat[k+1,j]
          break
        }
      }
      alpha.old <- alpha.new
      theta.old <- theta.new
      # tmp <- proc.time()[3]
    }
  }
  
  # cat("EM update: ",EM.time,"\n")
  # cat("Optimization for alpha: ",opt.time,"\n")
  # cat("Weight matrix: ",weight.time,"\n")
  # apply(M.stat,1,max)
  # M.stat
  # return(max(apply(M.stat,1,max)))
  # M.stat
  
  # Calculate EM statistics
  return(c("pvalue" = 1-pchisq(max(M.stat),3),"stat" = max(M.stat)))
}
