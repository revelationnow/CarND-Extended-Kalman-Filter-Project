#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
   * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() == 0)
	{
	    std::cout<<"Error"<<std::endl;
	    return rmse;
	}
	else if(estimations.size() != ground_truth.size())
	{
	    std::cout<<"Error"<<std::endl;
	    return rmse;
	}


	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i)
  {
    // ... your code here
    VectorXd diff = (estimations[i] - ground_truth[i]);
		diff = diff.array() * diff.array();
		rmse += diff;
    std::cout<<"Estimations : "<<estimations[i]<<std::endl;
    std::cout<<"Ground Truth : "<<ground_truth[i]<<std::endl;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  double px = x_state(0) ;
  double py = x_state(1) ;
  double vx = x_state(2) ;
  double vy = x_state(3) ;
/*
  if(fabs(px) < 0.00001)
    px = 0.00001;

  if(fabs(py) < 0.00001)
    py = 0.00001;
*/
  double c1 = px*px + py*py;
  double c2 = sqrt(c1);
  double c3 = (c1*c2);

  if(fabs(c1) < 0.0001)
  {
    std::cout<<"Divide by zero during Jacobian"<<std::endl;
    return Hj;
  }

  Hj << (px/c2)               , (py/c2)               , 0     , 0    ,
       -1.0*(py/c1)           , (px/c1)               , 0     , 0    ,
        py*(vx*py - vy*px)/c3 , px*(px*vy - py*vx)/c3 , px/c2 , py/c2;

  std::cout<<"\nJacobian matrix : \n" << Hj << std::endl << std::endl;
  return Hj;
}











