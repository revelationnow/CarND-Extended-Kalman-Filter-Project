#include "kalman_filter.h"
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = VectorXd(2);
  y = z - (H_ * x_);

  MatrixXd S = MatrixXd(2,2);
  S = (H_ * P_ * H_.transpose()) + R_;

  MatrixXd K = MatrixXd(4,2);
  K = P_ * H_.transpose() * S.inverse();

  x_ = x_ + (K * y);
  P_ = (MatrixXd::Identity(4,4) - (K * H_)) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  //z(1) = ((z(1) + 6.28318) % 6.28318) - 3.14159;
  double phi = fmod(z(1) + 6.28318, 6.28318) - 3.14159;

  VectorXd h_func_x = VectorXd(3);
  h_func_x(0) = sqrt(x_(0) * x_(0) + x_(1)*x_(1));
  h_func_x(1) = atan(x_(1)/x_(0));
  h_func_x(2) = (x_(0)*x_(2) + x_(1)*x_(3))/h_func_x(0);

  VectorXd y = VectorXd(3);
  y = z - h_func_x;
  y(1) = phi - h_func_x(1);
  std::cout<<"UpdateEKF : Calculated y"<<std::endl;

  MatrixXd S = MatrixXd(3,3);
  S = H_ * P_ * H_.transpose();
  std::cout<<"UpdateEKF : Calculated S"<<std::endl;

  MatrixXd K = MatrixXd(4,3);
  K = P_ * H_.transpose() * S.inverse();
  std::cout<<"UpdateEKF : Calculated K"<<std::endl;

  x_ = x_ + (K*y);
  P_ = (MatrixXd::Identity(4,4) - (K * H_)) * P_;
  std::cout<<"UpdateEKF : Calculated P"<<std::endl;

  
}








