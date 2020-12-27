#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // Create the rmse structure, initialized to 0s. 
  VectorXd rmse(4); 
  rmse << 0, 0, 0, 0;

  // Calculate RMSE over the whole measurement
  for(unsigned int i=0; i < estimations.size(); i++)
  {
     // Calculate residual which is difference of estimation to 
     // ground truth
     VectorXd residual = estimations[i]-ground_truth[i];

     // calculate residual squared
     residual = residual.array()*residual.array();

     // sum all values over the whole measurement so far
     rmse += residual;    
  }
  
  // divide to total number of sample points
  rmse = rmse/estimations.size(); 

  // take the square root to calculate root mean square
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  // Create the Jacobian matrix
  MatrixXd Hj(3,4);

  // Extract the states from state vector
  float px = x_state[0];
  float py = x_state[1];
  float vx = x_state[2];
  float vy = x_state[3]; 
  float pxySqSum = px*px + py*py;

  // Set Jacobian Matrix values according to the given structure
  Hj(0,0) = px/sqrt(pxySqSum);
  Hj(0,1) = py/sqrt(pxySqSum);
  Hj(0,2) = 0;
  Hj(0,3) = 0;
  Hj(1,0) = -1*py/pxySqSum;
  Hj(1,1) = px/pxySqSum;
  Hj(1,2) = 0;
  Hj(1,3) = 0;
  Hj(2,0) = (py*(vx*py-vy*px))/pow(pxySqSum, 1.5);
  Hj(2,1) = (px*(vy*px-vx*py))/pow(pxySqSum, 1.5);
  Hj(2,2) = px/sqrt(pxySqSum);
  Hj(2,3) = py/sqrt(pxySqSum); 

  return Hj; 
}
