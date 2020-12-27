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

  // measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // jacobian matrix
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
  
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    
    VectorXd x(4);
    x << 1, 1, 0, 0;
    MatrixXd P(4,4);
    P << 1, 0, 0, 0,
         0, 1, 0, 0, 
         0, 0, 1000, 0,
         0, 0, 0, 1000;
    
    MatrixXd F, H, R, Q;
    
    F = MatrixXd(4,4);
    F << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1;

    Q = MatrixXd(4,4);
    Q << 0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // get the radar measurement outputs, rho, psi and rhodot. 
      float rho = measurement_pack.raw_measurements_[0];
      float psi = measurement_pack.raw_measurements_[1];
      float rhodot = measurement_pack.raw_measurements_[2]; 

      // assign state vector
      x[0] = rho*cos(psi);
      x[1] = rho*sin(psi);
      // For initialization, although rhodot doesn't equate to the
      // actual speed of the vehicle, it may still provide a better estimate. 
      x[2] = rhodot*cos(psi);
      x[3] = rhodot*sin(psi);

      // assign R and H for Radar
      R = R_radar_;
      //H = tools.CalculateJacobian(x);
      H = Hj_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Assign state vector
      x[0] = measurement_pack.raw_measurements_[0];
      x[1] = measurement_pack.raw_measurements_[1];  
      
      // assign R and H for Laser
      R = R_laser_;
      H = H_laser_; 
    }

    // Initialize ekf with the set vector and matrices
    ekf_.Init(x, P, F, H, R, Q);
    
    // set the first timestamp to first measurement
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
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
  // assuming noise, and noise covariances don't change throughout measurement. 
  static float noise_ax = 9.0;
  static float noise_ay = 9.0;
  
  // Calculate variances and powers of dt for updating of the 
  // process noise covariance matrix Q
  static float var_ax = noise_ax*noise_ax; 
  static float var_ay = noise_ay*noise_ay;

  // Calculate powers of time delta for Q matrix
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

  // Call prediction step
  ekf_.Predict();

  /**
   * Update
   */

  // Update will be done based on sensor type
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // If sensor type is radar, then calculate the Jacobian matrix
    // according to the latest state, and set R to measurement 
    // covariance matrix of radar. Then call updateEKF function. 

    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_; 
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // If sensor type is lidar, then set H and R to measurement matrix
    // and measurement covariance matrix of lidar. then call update 
    // function. 
    ekf_.H_ = H_laser_; 
    ekf_.R_ = R_laser_; 
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ after update = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
