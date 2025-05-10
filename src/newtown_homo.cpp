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

double sum_partial_log_rr_2(Eigen::VectorXd x, Eigen::VectorXd w,
                          Eigen::VectorXd reduce_num,
                          double m, double r, int max_x){
  double ret = 0;
  int n = x.size();
  int cut = max_x + 1;
  Eigen::VectorXd sum_Process;
  sum_Process.setZero(cut);
  sum_Process(0) = -1/(r*r);
  for(int i=1; i< cut; i++){
    sum_Process(i) = sum_Process(i-1) + (-1/((r+i)*(r+i)));
  }
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
double sum_partial_r_2(Eigen::VectorXd x, Eigen::VectorXd w,
                     Eigen::VectorXd reduce_num,
                     double m, double r, int max_x){
  double ret = 0;
  int n = x.size();
  int cut = max_x + 1;
  Eigen::VectorXd sum_Process;
  sum_Process.setZero(cut);
  sum_Process(0) = 1/r;
  for(int i=1; i< cut; i++){
    sum_Process(i) = sum_Process(i-1) + 1/(r+i);
  }
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

// [[Rcpp::export]]
Eigen::VectorXd homo_opt(Eigen::VectorXd w_reduce,
                                Eigen::VectorXd x_reduce, Eigen::VectorXd reduce_num, 
                                int newtown_step = 50, double th = 1e-10, double r=5){
  int max_x = x_reduce.maxCoeff() + 2;
  
  double m = (x_reduce.array() * w_reduce.array() * reduce_num.array()).matrix().sum() /
    (w_reduce.array() * reduce_num.array()).matrix().sum();

  
  Eigen::VectorXd ret;
  ret.setZero(2);
  ret(0) = m;
  double r_cur = r;
  while(sum_partial_r_2(x_reduce,w_reduce,reduce_num,m,r_cur,max_x) < 0){
    r_cur = r_cur/2;
  }
  while (newtown_step > 0) {
    newtown_step --;
    double deri = sum_partial_r_2(x_reduce,w_reduce,reduce_num,m,r_cur,max_x);
    if(abs(deri) < th) break;
    else  r_cur = r_cur - ((1.0/sum_partial_log_rr_2(x_reduce,w_reduce,reduce_num,m,r_cur,max_x)) * deri);
    if(r_cur > 1e3) r_cur = 1e3;
    if(r_cur < 0) r_cur = 1e-50;
    
  }
  ret(1) = r_cur;

  return ret;
}
