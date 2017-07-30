#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
	radarmeasurements = 0;
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
  H_laser_ << 1, 0, 0, 0,
				0, 1, 0, 0;
    VectorXd xin = VectorXd(4);
      MatrixXd	P = MatrixXd(4, 4);
    P << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;
      MatrixXd  F = MatrixXd(4, 4);
	F << 1, 0, 1, 0,
		 0, 1, 0, 1,
		 0, 0, 1, 0,
		 0, 0, 0, 1;

	MatrixXd Q = MatrixXd(4, 4);
	Q << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
	ekf_ = KalmanFilter(xin, P, F, H_laser_, R_radar_, Q);
  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
    noise_ax = 9;
    noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

cout << is_initialized_ << endl;
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



    cout << "Before sensor check." << endl;
    // first measurement
    cout << "EKF: " << endl;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      float ro,theta, ro_dot;

      ro = measurement_pack.raw_measurements_[0];
      theta = measurement_pack.raw_measurements_[1];
      ro_dot = measurement_pack.raw_measurements_[2];
      float px,py;
      px = ro*cos(theta);
      py = ro*sin(theta);

      previous_timestamp_ = measurement_pack.timestamp_;
		ekf_.x_  = VectorXd(4);
      ekf_.x_  << px, py, 0, 0;
	  radarmeasurements++;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
		
		previous_timestamp_ = measurement_pack.timestamp_;
		

    }
    cout << "Done initializing." << endl;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  cout << "time stamp " << measurement_pack.timestamp_ << endl;
  //There is a chance of truncation/ lost information here but is negligable
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;

	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ <<  dt_4/4*pow(noise_ax,2), 0, dt_3/2*pow(noise_ax,2), 0,
			   0, dt_4/4*pow(noise_ay,2), 0, dt_3/2*pow(noise_ay,2),
			   dt_3/2*pow(noise_ax,2), 0, dt_2*pow(noise_ax,2), 0,
			   0, dt_3/2*pow(noise_ay,2), 0, dt_2*pow(noise_ay,2);
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
    ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;
	if(radarmeasurements == 1 && measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
	{
	 float theta = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];
	  
		cout << "first Vx " << (ekf_.x_[2] = ro_dot*cos(theta)) << endl;
		cout << "first Vy " << (ekf_.x_[3] = ro_dot*sin(theta)) << endl;
		radarmeasurements ++;
	}
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

   if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
     // radar updates
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
   } else {
     // laser updates
     //ekf_.R_ = R_laser_;
     //ekf_.Update(measurement_pack.raw_measurements_);
   }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
 
 
}
