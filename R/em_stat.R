em.stat <- function(x,alpha.ini,k0,C,labels,group.num,prior.weight = 0.05,earlystop = NA, is_ecdf = F, 
                    cut_max = 50){
  # require(nloptr)
  # C_time = Sys.time()
  
  # Calculate EM statistics
  # k0 is the number of iterations, and alpha.ini are the initial alpha values
  # C is the coefficient of the penalty term
  # The first n components of `theta` are means, while the last n components are dispersion parameters (r)
  # Labels should be in the form of 1,2,...,n
  # Each row of alpha.ini contains a set of initial values
  
  x_table = as.matrix(table(x))
  x_reduce = as.numeric(rownames(x_table))
  reduce_num = as.vector(x_table[,1])
  if (length(unique(labels))!=group.num){
    print("The number of different labels does not match the total number of groups!")
    return(NULL)
  }
  M.stat <- matrix(nrow = 1+k0,ncol = nrow(alpha.ini))
  w_homo = rep(1/length(x),length(reduce_num))
  m_homo = mean(x)
  var_homo = var(x)
  fit.nb <- function(x){
    # this function returns the MLE of mean and r
    obj.f <- function(theta){-sum(dnbinom(x,theta[2],mu = theta[1],log = T))}
    deriv.pl <- function(theta){nl.grad(x0 = theta,fn = obj.f)}
    nloptr(x0 = c(min(max(mean(x),1e-9*1.01),2999),5),eval_f = obj.f,lb=rep(1e-9,2),ub = rep(3000,2),opts = list("algorithm"="NLOPT_LN_NELDERMEAD","ftol_rel" = 1e-9,"maxeval" = 5000,"maxtime" = 200,"xtol_rel" = 1e-4))$solution
  }
  if(cut_max != max(x)){
    theta.0 = homo_opt(w_homo, x_reduce, reduce_num,  newtown_step = 100, r = 1)
    theta_hat_0 = c(rep(theta.0[1],group.num),rep(theta.0[2],group.num))
  } else {
    theta.0 <- fit.nb(x)
    theta_hat_0 = c(rep(theta.0[1],group.num),rep(theta.0[2],group.num))
  }
  pl0 = compute_pl(theta_hat_0, as.integer(x_reduce), reduce_num, rep(1/group.num,group.num), group.num, C)

  mu <- numeric(group.num)
  r <- numeric(group.num)
  if(cut_max != max(x)){
    for (i in 1:group.num) {
    x_i = x[labels==i]
    reduce_num_i = rep(1, length(x_i))
    w_homo_i = rep(1/length(x_i), length(x_i))
    m_homo = mean(x_i)
    var_homo = var(x_i)
    theta.0 = homo_opt(w_homo_i, x_i, reduce_num_i,  newtown_step = 100, r = 1)
    mu[i] <- theta.0[1]
    r[i] <- theta.0[2]
    }
  } else {
    for (i in 1:group.num){
      mu[i] <- mean(x[labels==i])
      r[i] <- fit.nb(x[labels==i])[2]
    }
  }

  mu[mu<=1e-9] <- 1e-9*1.01
  r[r<=1e-9] <- 1e-9*1.01
  x_order = order(x, decreasing = F)
  x = x[x_order]
  labels = labels[x_order]

  
  for (j in 1:nrow(alpha.ini)){
    # Calculate the initial M values
    
    k <- 0
    alpha.old <- alpha.ini[j,]

    theta.old <- c(mu,r)
    pl.value = compute_pl(theta.old, as.integer(x_reduce), reduce_num, c(alpha.old, 1-sum(alpha.old)), group.num, C)
    M.stat[1,j] <- 2 * (pl.value-pl0)
    w0 = matrix(0.001,nrow = length(x),ncol = group.num)
    for (i in 1:group.num){
      w0[labels==i,i] <- 1-0.001*(group.num-1)
    }
    w <- w0
    w0_reduce = matrix(0, ncol = group.num, nrow = length(reduce_num))
    for (g in 1:group.num) {
      w_g = w0[,g]
      w0_reduce[,g] = reduce_w(w_g, reduce_num)
    }
    w_orign = w0
    w_reduce = w0_reduce
    
    ret = FSEM(group.num, prior.weight, C , w_reduce, w0_reduce, x_reduce,as.integer(x_reduce), reduce_num,theta.old, k0,50,cut_max)
    # print(Sys.time() - C_time)
    M.stat[2:(k0+1),j] = 2*(ret- pl0)
  
  }
  # Calculate EM statistics
  return(c("pvalue" = 1-pchisq(max(M.stat),3),"stat" = max(M.stat)))
}
