#include "kalman_filter.h"
#include <iostream>
#include "tools.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
						 // state vector
  x_ = x_in;

  // state covariance matrix
  P_ = P_in;

  // state transition matrix
  F_ = F_in;

  // process covariance matrix
  Q_ = Q_in;

  // measurement matrix
  H_ = H_in;

  // measurement covariance matrix
  R_ = R_in;
						
						}
KalmanFilter::KalmanFilter(){}
KalmanFilter::~KalmanFilter() {}


void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
        x_ = F_ * x_;
		MatrixXd Ft = F_.transpose();
		P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  Tools tools;
  VectorXd zmeasurement = z;
  VectorXd polarMeasurement = VectorXd(3);
  float theta;
  float measurementTheta = z[1];
  while(measurementTheta > M_PI)
  {
  measurementTheta -= 2*M_PI;
  }
  while(measurementTheta < -M_PI)
  {
  measurementTheta += 2*M_PI;
  }
  theta = atan2(x_[1],x_[0]);
  zmeasurement[1] = measurementTheta;
  while(theta > M_PI)
  {
  theta -= 2*M_PI;
  }
  while(theta < -M_PI)
  {
  theta += 2*M_PI;
  }
  float p = sqrt( pow(x_[0], 2)+ pow(x_[1], 2));
  
  
  polarMeasurement << p,theta, ((x_[0]*x_[2]+x_[1]*x_[3])/p);

  MatrixXd Hj = tools.CalculateJacobian(x_);
 	VectorXd y = zmeasurement - polarMeasurement;
	while(y[1] > M_PI)
	{
	y[1] -= 2*M_PI;
	}
	while(y[1] < -M_PI)
	{
	y[1] += 2*M_PI;
	}
	MatrixXd Ht = Hj.transpose();
	MatrixXd S = Hj* P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj) * P_;
}


