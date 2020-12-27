#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
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
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1;

  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << 0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0;

  // measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_ << 0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  Tools tools; 
  
  cout << "processing measurement" << endl;
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0, 
               0, 0, 1000, 0,
               0, 0, 0, 1000;
    
    float px, py;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      float rho = measurement_pack.raw_measurements_[0];
      float psi = measurement_pack.raw_measurements_[1];
      float rhodot = measurement_pack.raw_measurements_[2]; 

      cout << "rho: " << rho << ", and psi: " << psi << endl;

      px = rho*cos(psi);
      py = rho*sin(psi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];  
    }
    
    ekf_.x_ << px, py, 1, 1;

    previous_timestamp_ = measurement_pack.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    cout << "initialized, x_ is: " <<  ekf_.x_ << endl;
    return;
  }

  /**
   * Prediction
   */

  // Calculate the time delta from the previous timestamp
  float timeUnit = 1000000.0;
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/timeUnit;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // Set the time dependent values of state transition matrix F
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
  
  // noise for ax and ay are given as 9
  float noise_ax = 9.0;
  float noise_ay = 9.0;
  
  // Calculate variances and powers of dt for updating of the 
  // process noise covariance matrix Q
  float var_ax = noise_ax*noise_ax; 
  float var_ay = noise_ay*noise_ay;

  float dt2 = dt*dt;
  float dt3 = dt2*dt;
  float dt4 = dt3*dt; 

  // Set the values of the process noise covariance matrix Q
  ekf_.Q_(0,0) = (dt4/4)*var_ax;
  ekf_.Q_(0,2) = (dt3/2)*var_ax;
  ekf_.Q_(1,1) = (dt4/4)*var_ay;
  ekf_.Q_(1,3) = (dt3/2)*var_ay;
  ekf_.Q_(2,0) = (dt3/2)*var_ax;
  ekf_.Q_(2,2) = dt2*var_ax;
  ekf_.Q_(3,1) = (dt3/2)*var_ay;
  ekf_.Q_(3,3) = dt2*var_ay;

  cout << ekf_.Q_ << endl;
  // Call prediction step
  ekf_.Predict();

  cout << "x_ after predict = " << ekf_.x_ << endl;
  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  cout << "meas: " << measurement_pack.raw_measurements_ << endl;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_; 
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_; 
    ekf_.R_ = R_laser_; 
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ after update = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
