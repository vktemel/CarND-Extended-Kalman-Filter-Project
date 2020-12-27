#include "kalman_filter.h"
#include <iostream>
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
  // To predict the state, we need to use the following equations. 
  // X' = F X; where v is random acceleration vector
  // P' = F P F + Q; where Q is the process noise covariance matrix
  x_ = F_ * x_; 
  
  MatrixXd Ftranspose = F_.transpose();
  P_ = F_ * P_ * Ftranspose + Q_; 
}

void KalmanFilter::Update(const VectorXd &z) {
  // Calculate the error
  VectorXd y = z - (H_ * x_); 
  
  // Call UpdateStates with the error as input
  UpdateStates(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // Calculates the error according to the following
  // y = z - h(x')

  // calculate h(x') given the x_ vector
  VectorXd hx(3); 
  hx << sqrt((x_[0]*x_[0])+(x_[1]*x_[1])),
        atan2(x_[1], x_[0]),
        ((x_[0]*x_[2])+(x_[1]*x_[3]))/sqrt((x_[0]*x_[0])+(x_[1]*x_[1]));

  VectorXd y = z - hx;
  // values ot y[1] can be over pi or below -pi; however, EKF expects
  // this value to be between -pi and pi. Therefore, if it's above pi, 
  // 2*pi is subtracted, and if it's below -pi, 2*pi is added. 
  if(y[1] > M_PI)
  {
    y(1) -= 2*M_PI;
  }
  else if(y[1] < -1*M_PI)
  {
    y[1] += 2*M_PI;
  }

  // Call UpdateStates with the error as input
  UpdateStates(y);
}

void KalmanFilter::UpdateStates(const VectorXd &y){
  // This is the Kalman Filter Update which applies to Matrix
  // equations for the Kalman Filter. 
  MatrixXd Htranspose = H_.transpose();
  MatrixXd S = H_ * P_ * Htranspose + R_;
  MatrixXd Sinverse = S.inverse(); 
  MatrixXd K = P_ * Htranspose * Sinverse; 

  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
