#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

  // calculate predicted output z_pred from output matrix H and state x
  VectorXd z_pred = H_ * x_;

  // calculate y error between measurement z and predicted output z_pred
  VectorXd y = z - z_pred;

  UpdateStateAndCovarianceMatrix(y);
}

void KalmanFilter::UpdateStateAndCovarianceMatrix(const Eigen::VectorXd &y) {

  // calculate Kalman gain K from output matrix H, state covariance matrix P and
  // measurement noise matrix R
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // update state estimate x using Kalman gain K and y error
  x_ = x_ + (K * y);

  // update state covariance matrix P from Kalman gain K and measurement
  // matrix H
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd z_pred(3);

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  // calculate components of z_pred from the state and the non linear
  // measurement function h
  float rho = sqrt(px * px + py * py);
  // correction if calculated rho is too small
  rho = max(rho, 0.001F);

  float theta = atan2(py, px);
  float rho_dot = (px * vx + py * vy) / rho;

  z_pred << rho, theta, rho_dot;

  // calculate y error between measurement z and predicted output z_pred
  VectorXd y = z - z_pred;

  // set calculated angle component of y error between -PI and +PI
  float rho_corrected = y(1);
  while (rho_corrected > M_PI) {
    rho_corrected -= M_PI;
  }
  while (rho_corrected < -M_PI) {
    rho_corrected += M_PI;
  }
  y(1) = rho_corrected;

  UpdateStateAndCovarianceMatrix(y);
}
