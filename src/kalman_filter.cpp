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
  std::cout << "\n z = \n" << z << std::endl << std::endl;
  VectorXd y = VectorXd(2);
  y = z - (H_ * x_);

  std::cout << "\n y = \n" << y << std::endl << std::endl;

  MatrixXd S = MatrixXd(2,2);
  S = (H_ * P_ * H_.transpose()) + R_;

  std::cout << "\n S =\n " << S << std::endl << std::endl;
  MatrixXd K = MatrixXd(4,2);
  K = P_ * H_.transpose() * S.inverse();

  std::cout << "\n K = \n" << K << std::endl << std::endl;
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
  std::cout << "\n z = " << z << std::endl << std::endl;

  VectorXd h_func_x = VectorXd(3);
  h_func_x(0) = sqrt(x_(0) * x_(0) + x_(1)*x_(1));
  h_func_x(1) = atan2(x_(1),x_(0));
  h_func_x(2) = (x_(0)*x_(2) + x_(1)*x_(3))/h_func_x(0);

  VectorXd y = VectorXd(3);
  y = z - h_func_x;
  y(1) = fmod(phi - h_func_x(1) + 6.28318, 6.28318) - 3.14159;

  std::cout<< "\ny = \n" << y << std::endl<< std::endl;

  MatrixXd S = MatrixXd(3,3);
  S = H_ * P_ * H_.transpose();

  std::cout<< "\nS = \n" << S << std::endl<< std::endl;

  MatrixXd K = MatrixXd(4,3);
  K = P_ * H_.transpose() * S.inverse();

  std::cout<< "\nK = \n" << K << std::endl<< std::endl;
  std::cout<< "\nH transpose = \n" << H_.transpose() << std::endl<< std::endl;
  std::cout<< "\nS inverse = \n" << S.inverse() << std::endl<< std::endl;


  x_ = x_ + (K*y);
  P_ = (MatrixXd::Identity(4,4) - (K * H_)) * P_;

  
}








