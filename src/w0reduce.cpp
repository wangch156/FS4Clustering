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


// [[Rcpp::export]]
Eigen::VectorXd  reduce_w(Eigen::VectorXd  w, Eigen::VectorXd  reduce_num) {
  int n = reduce_num.size();
  Eigen::VectorXd w_reduce;
  w_reduce.setZero(n);
  int sum_w = 0;
  
  for(int i=0; i<n; i++){
    int num = int(reduce_num(i));
    for(int j =sum_w; j < (sum_w + num); j++){
      w_reduce(i) += w(j);
    }
    double sum_i = w_reduce(i);
    w_reduce(i) = sum_i / num;
    sum_w += num;
  }
  return w_reduce;
}

