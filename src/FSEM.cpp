// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <math.h>
#include <algorithm>


using namespace Rcpp;
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

double sum_partial_log_rr(Eigen::VectorXd x, Eigen::VectorXd w,
                          Eigen::VectorXd reduce_num,
                          double m, double r, int max_x, int cut_max = 50){
  double ret = 0;
  int n = x.size();
  int cut = min(cut_max, max_x) + 1;
  Eigen::VectorXd sum_Process;
  sum_Process.setZero(cut);
  sum_Process(0) = -1/(r*r);
  for(int i=1; i< cut; i++){
    sum_Process(i) = sum_Process(i-1) + (-1/((r+i)*(r+i)));
  }
  // cout << sum_Process << endl;
  for(int i=0; i<n; i++){
    int tmp = int(x[i]);
    if(tmp == 0)  ret += 0;
    else if(tmp <= cut) ret +=  w(i) * reduce_num(i) *sum_Process(tmp-1);
    else {
      ret += w(i) * reduce_num(i) * (sum_Process(cut-1) + (1/(r+tmp-1) - 1/(r+cut-1)) +
        0.5 * (- 1/((r+tmp-1)*(r+tmp-1)) + 1/((r+cut-1)*(r+cut-1))) 
                                       + 2/6 * (1/((r+tmp-1)*(r+tmp-1)*(r+tmp-1)) 
                                                  - 1/((r+cut-1)*(r+cut-1)*(r+cut-1))));
                                                  
    }
    ret += w(i) * reduce_num(i)  * (1/r - 1/(r+m));
  }
  return ret;
}
double sum_partial_r(Eigen::VectorXd x, Eigen::VectorXd w,
                     Eigen::VectorXd reduce_num,
                     double m, double r, int max_x, int cut_max = 50){
  double ret = 0;
  int n = x.size();
  int cut = min(cut_max, max_x) + 1;
  Eigen::VectorXd sum_Process;
  sum_Process.setZero(cut);
  sum_Process(0) = 1/r;
  for(int i=1; i< cut; i++){
    sum_Process(i) = sum_Process(i-1) + 1/(r+i);
  }
  // cout << sum_Process << endl;
  for(int i=0; i<n; i++){
    int tmp = int(x[i]);
    if(tmp == 0)  ret += 0;
    else if(tmp <= cut) ret +=  w(i) * reduce_num(i) * sum_Process(tmp-1);
    else {
      ret += w(i) * reduce_num(i) * (sum_Process(cut-1) + (log(r+tmp-1) - log(r+cut-1)) +
        0.5 * (1/(r+tmp-1) - 1/(r+cut-1)) + 1/6 * (-1/((r+tmp-1)*(r+tmp-1))
                                                     + 1/((r+cut-1)*(r+cut-1))));
      
      
    }
    ret += w(i) * reduce_num(i) * (log(r / (r+m)));
    
  }
  return ret;
}

Eigen::VectorXd EM_update_theta(double r, Eigen::VectorXd w_reduce, Eigen::VectorXd w0_reduce,
                                     Eigen::VectorXd x_reduce, Eigen::VectorXd reduce_num, 
                                     int newtown_step = 50, double th = 1e-6, int cut_max = 50){
  int max_x = x_reduce.maxCoeff() + 2;
  w_reduce += w0_reduce;
  
  double m = (x_reduce.array() * w_reduce.array() * reduce_num.array()).matrix().sum() /
    (w_reduce.array() * reduce_num.array()).matrix().sum();
  Eigen::VectorXd ret;
  ret.setZero(2);
  ret(0) = m;
  double r_cur = r;
  while(sum_partial_r(x_reduce,w_reduce,reduce_num,m,r_cur,max_x,cut_max) < 0){
    r_cur = r_cur/2;
  }
  while (newtown_step > 0) {
    newtown_step --;
    double deri = sum_partial_r(x_reduce,w_reduce,reduce_num,m,r_cur,max_x,cut_max);
    if(abs(deri) < th) break;
    else  r_cur = r_cur - ((1.0/sum_partial_log_rr(x_reduce,w_reduce,reduce_num,m,r_cur,max_x,cut_max)) * deri);
    if(r_cur > 1e3) r_cur = 1e3;
  }
  ret(1) = r_cur;
  return ret;
}

Eigen::VectorXd compute_alpha(Eigen::MatrixXd w, Eigen::VectorXd reduce_num, 
                              double len, double C, int group_num) {
  
  Eigen::VectorXd alpha;
  alpha.setZero(group_num);
  for(int g = 0; g<group_num; g++){
    Eigen::VectorXd w_g = w.col(g);
    
    alpha(g) = ((w_g.array() * reduce_num.array()).matrix().sum() + C) / 
      (len + double(group_num) * C);
  }
  return alpha;
}

Eigen::VectorXd compute_f(double m, double r, Eigen::VectorXi x){
  int max_x = x.maxCoeff() + 5;
  int n = x.size();
  Eigen::VectorXd Pre_pocess;
  Pre_pocess.setZero(max_x);
  double p1 = m/(r+m);
  double p2 = r/(r+m);
  Pre_pocess(0) = pow(p2, r);
  for(int i = 1; i<max_x ; i++){
    Pre_pocess(i) = Pre_pocess(i-1) * (double)((r+i-1) / (double)(i)) * p1;
  }
  Eigen::VectorXd density;
  density.setZero(n);
  for(int i = 0; i<n ; i++){
    int tmp = x(i);
    density(i) = Pre_pocess(tmp);
  }
  return(density);
}

double compute_R1n (Eigen::VectorXd theta, Eigen::VectorXi x, Eigen::VectorXd reduce_num,
                    Eigen::VectorXd alpha, int group_num, double C){
  double R1n = 0.0;
  int n = x.size();
  Eigen::MatrixXd fg;
  fg.setZero(group_num, n);
  Eigen::VectorXd f;
  f.setZero(n);
  Eigen::VectorXd alpha_new;
  alpha_new = alpha;
  for(int g=0; g<group_num ; g++){
    fg.row(g) = compute_f(theta(g), theta(g+group_num), x);
    f +=  (alpha_new(g) * fg.row(g).array()).matrix();
  }
  R1n =  double((f.array().log() * reduce_num.array()).matrix().sum() + 
    C*alpha_new.array().log().matrix().sum());
  return (R1n);
}

Eigen::MatrixXd compute_w(Eigen::VectorXd theta, Eigen::VectorXi x, Eigen::VectorXd alpha, 
                          int group_num, double prior_weight){
  int n = x.size();
  Eigen::MatrixXd w;
  w.setZero(n,group_num);
  Eigen::MatrixXd fg;
  fg.setZero(group_num, n);
  Eigen::VectorXd f;
  f.setZero(n);
  Eigen::VectorXd alpha_new;
  alpha_new = alpha;
  for(int g=0; g<group_num ; g++){
    fg.row(g) = compute_f(theta(g), theta(g+group_num), x);
    f +=  (alpha_new(g) * fg.row(g).array()).matrix();
  }
  for(int g=0; g<group_num; g++){
    Eigen::VectorXd tmp_g = alpha_new(g) * fg.row(g);
    w.col(g) = (1-prior_weight) * (tmp_g.array() / f.array()).matrix();
  } 
  return(w);
}


// [[Rcpp::export]]
Eigen::VectorXd FSEM(int group_num, double prior_wight, double C,  Eigen::MatrixXd w_reduce, 
                     Eigen::MatrixXd w0_reduce,  Eigen::VectorXd x_reduce, Eigen::VectorXi x_reduce_int,
                     Eigen::VectorXd reduce_num, Eigen::VectorXd theta_old,
                     int k_step = 30, int newtown_step = 50, int cut_max = 50) {
  VectorXd R1n;
  R1n.setZero(k_step);
  VectorXd theta_new; 
  theta_new.setZero(2*group_num);
  VectorXd alpha; 
  alpha.setZero(group_num);
  for(int k=0; k<k_step; k++){
    for(int i=0; i<group_num; i++){
      VectorXd w_reduce_i=w_reduce.col(i);
      VectorXd w0_reduce_i=w0_reduce.col(i);
      if(k == 0){
        VectorXd ret;
        ret = EM_update_theta(theta_old(group_num+i),0*w_reduce_i,w0_reduce_i,x_reduce,reduce_num, 50,1e-6,cut_max);
        theta_new(i) = ret(0);
        theta_new(i+group_num) = ret(1);
        // cout << theta_old<< "theta_old\n";
      } else { 
        VectorXd ret;
        ret = EM_update_theta(theta_old(group_num+i),w_reduce_i,prior_wight * w0_reduce_i,x_reduce,reduce_num, 50,1e-6,cut_max);
        theta_new(i) = ret(0);
        theta_new(i+group_num) = ret(1);
      }
    }
    
    // regularization
    for(int i=0; i<(2*group_num); i++){
      if(theta_new(i) < 1e-10){
        theta_new(i) = 1e-10;
      }
    }
    
    // update alpha
    MatrixXd w_reduce_all; 
    if(k == 0){
      w_reduce_all = w0_reduce;
    } else {
      w_reduce_all = w_reduce + prior_wight * w0_reduce;
    }
    double len = reduce_num.sum();
    alpha = compute_alpha(w_reduce_all,reduce_num,len,C,group_num);
    
    
    // update w_reduce
    w_reduce = compute_w(theta_new, x_reduce_int, alpha, group_num, prior_wight);
    theta_old = theta_new;
    
    // traj
    R1n(k) = compute_R1n(theta_new, x_reduce_int, reduce_num,
         alpha, group_num,C);
    // cout << R1n << "R1nvec\n"; 
  }
  return R1n;
}
