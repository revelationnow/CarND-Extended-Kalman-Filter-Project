#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
*/FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_ << 1, 1, 1, 1,
         1, 1, 1, 1,
         1, 1, 1, 1;

  ekf_.noise_ax_ = 9;
  ekf_.noise_ay_ = 9;

  ekf_.Q_ = MatrixXd(4,4);
  ekf_.F_ = MatrixXd(4,4);
  ekf_.P_ = MatrixXd(4,4);
  ekf_.x_ = MatrixXd(4,1);
  ekf_.H_ = MatrixXd(2,4);
  ekf_.R_ = MatrixXd(2,2);

  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  ekf_.Q_ << 0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: First measurement" << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      cout << "EKF: First measurement was RADAR" << endl;
      double rho = measurement_pack.raw_measurements_(0);
      double phi = measurement_pack.raw_measurements_(1);
      double rho_dot = measurement_pack.raw_measurements_(2);
      ekf_.x_(0) = rho * cos(phi);
      ekf_.x_(1) = rho * sin(phi);
      ekf_.x_(2) = rho_dot * cos(phi);
      ekf_.x_(3) = rho_dot * sin(phi);


    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "EKF: First measurement was LASER" << endl;
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
    }
    ekf_.P_ << 100000000, 1.00000000, 1.00000000, 1.00000000,
               1.00000000, 100000000, 1.00000000, 1.00000000,
               1.00000000, 1.00000000, 100000000, 1.00000000,
               1.00000000, 1.00000000, 1.00000000, 100000000;
    last_time_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double elapsed_time = ((double)measurement_pack.timestamp_ - (double)last_time_)/1000000.0;
  last_time_ = measurement_pack.timestamp_;
  cout << "Elapsed time = " << elapsed_time << endl;
  ekf_.F_ << 1, 0, elapsed_time, 0,
             0, 1, 0, elapsed_time,
             0, 0, 1, 0,
             0, 0, 0, 1;
  ekf_.F_(0,2) = elapsed_time;
  ekf_.F_(1,3) = elapsed_time;

  ekf_.Q_ << pow(elapsed_time,4) * ekf_.noise_ax_ / 4, 0, pow(elapsed_time,3) * ekf_.noise_ax_ / 2, 0,
             0, pow(elapsed_time,4) * ekf_.noise_ay_ / 4, 0, pow(elapsed_time,3) * ekf_.noise_ay_ / 2,
             pow(elapsed_time,3) * ekf_.noise_ax_ / 2, 0, pow(elapsed_time,2) * ekf_.noise_ax_ / 1, 0,
             0, pow(elapsed_time,3) * ekf_.noise_ay_ / 2, 0, pow(elapsed_time,2) * ekf_.noise_ay_ / 1;
  ekf_.Q_(0,0) = pow(elapsed_time,4) * ekf_.noise_ax_ / 4;
  ekf_.Q_(0,2) = pow(elapsed_time,3) * ekf_.noise_ax_ / 2;
  ekf_.Q_(1,1) = pow(elapsed_time,4) * ekf_.noise_ay_ / 4;
  ekf_.Q_(1,3) = pow(elapsed_time,3) * ekf_.noise_ay_ / 2;

  ekf_.Q_(2,0) = pow(elapsed_time,3) * ekf_.noise_ax_ / 2;
  ekf_.Q_(2,2) = pow(elapsed_time,2) * ekf_.noise_ax_ / 1;
  ekf_.Q_(3,1) = pow(elapsed_time,3) * ekf_.noise_ay_ / 2;
  ekf_.Q_(3,3) = pow(elapsed_time,2) * ekf_.noise_ay_ / 1;

  cout << "Before Predict : " << endl << endl;
  cout << "x_ = \n" << ekf_.x_ << endl << endl;
  cout << "P_ = \n" << ekf_.P_ << endl << endl;
  cout << "F_ = \n" << ekf_.F_ << endl << endl;
  cout << "Q_ = \n" << ekf_.Q_ << endl << endl;
  ekf_.Predict();
  cout << "-----------------------------------------"<< endl << endl;
  cout << "After Predict : " << endl << endl;
  cout << "x_ = \n" << ekf_.x_ << endl << endl;
  cout << "P_ = \n" << ekf_.P_ << endl << endl;
  cout << "-----------------------------------------"<< endl << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    cout << "Doing Radar Update" << endl;
    Hj_ = tools.CalculateJacobian(measurement_pack.raw_measurements_); 
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    cout << "Doing Laser Update" << endl;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "After measurement : " << endl;
  cout << "x_ = \n" << ekf_.x_ << endl;
  cout << "P_ = \n" << ekf_.P_ << endl ;
  cout << "-----------------------------------------";
  cout << "-----------------------------------------";
  cout << "-----------------------------------------";
  cout << endl << endl << endl;
}
