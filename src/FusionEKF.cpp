#include "Eigen/Dense"
#include "FusionEKF.h"
#include "tools.h"
#include <iostream>

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

  // initialize matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  // initialize laser measurement matrix
  H_laser_ << 1, 0, 0, 0, 0, 1, 0, 0;

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0, 0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;

  // set the process noise
  noise_ax_ = 9.0;
  noise_ay_ = 9.0;

  // create x, F, P and Q
  ekf_.x_ = Eigen::VectorXd(4);
  ekf_.F_ = Eigen::MatrixXd(4, 4);
  ekf_.P_ = Eigen::MatrixXd(4, 4);
  ekf_.Q_ = Eigen::MatrixXd(4, 4);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {

    // initialize the state covariance matrix (for the velocity component we set
    // high covariance because they're not reliably calculated from first
    // measurement)
    ekf_.P_ << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1000, 0, 0, 0, 0, 1000;

    // initialize state matrix (in the prediction step we will modify it based
    // on delta time)
    ekf_.F_ << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      // meas_package.raw_measurements_ from RADAR contains: rho, theta, rho_dot
      // v components are set to 0 because they cannot be reliably  calculated
      // from one single measurement
      float x = measurement_pack.raw_measurements_[0] *
                cos(measurement_pack.raw_measurements_[1]);
      float y = measurement_pack.raw_measurements_[0] *
                sin(measurement_pack.raw_measurements_[1]);
      float vx = 0;
      float vy = 0;
      ekf_.x_ << x, y, vx, vy;
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state.
      // meas_package.raw_measurements_ from LASER contains px, py;
      ekf_.x_ << measurement_pack.raw_measurements_[0],
          measurement_pack.raw_measurements_[1], 0, 0;
    }

    // set previous timestamp to latest time stamp from measurement pack
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  // calculate delta time in seconds between previous timestamp and new
  // measurement
  float delta_t =
      (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  float delta_t_2 = delta_t * delta_t;
  float delta_t_3 = delta_t_2 * delta_t;
  float delta_t_4 = delta_t_3 * delta_t;

  // set previous timestamp to latest time stamp from measurement pack
  previous_timestamp_ = measurement_pack.timestamp_;

  // set state transition matrix F based on delta time
  ekf_.F_(0, 2) = delta_t;
  ekf_.F_(1, 3) = delta_t;

  // set state covariance matrix Q based on delta time
  ekf_.Q_ << (delta_t_4 / 4) * noise_ax_, 0, (delta_t_3 / 2) * noise_ax_, 0, 0,
      (delta_t_4 / 4) * noise_ay_, 0, (delta_t_3 / 2) * noise_ay_,
      (delta_t_3 / 2) * noise_ax_, 0, delta_t_2 * noise_ax_, 0, 0,
      (delta_t_3 / 2) * noise_ay_, 0, delta_t_2 * noise_ay_;

  ekf_.Predict();

  /**
   * Update
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // for a radar measurement model, the h is a non linear function, therefore
    // it has to be linearized around the state x
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

    ekf_.R_ = R_radar_;

    // update EKF with H and R for a radar measurement
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }

  else {
    // for a laser measurement model, H is a matrix (linear function between
    // state x and measurement y)
    ekf_.H_ = H_laser_;

    ekf_.R_ = R_laser_;

    // update EKF with H and R for a radar measurement
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  std::cout << "x_ = " << ekf_.x_ << std::endl;
  std::cout << "P_ = " << ekf_.P_ << std::endl;
}
