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

Eigen::VectorXd compute_ff(double m, double r, Eigen::VectorXi x){
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

// [[Rcpp::export]]
double compute_pl (Eigen::VectorXd theta, Eigen::VectorXi x, Eigen::VectorXd reduce_num,
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
    fg.row(g) = compute_ff(theta(g), theta(g+group_num), x);
    f +=  (alpha_new(g) * fg.row(g).array()).matrix();
  }
  R1n =  double((f.array().log() * reduce_num.array()).matrix().sum() + 
    C*alpha_new.array().log().matrix().sum());
  return (R1n);
}
