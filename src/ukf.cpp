#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
double NormalizeAngle(double &angle)
{
  while (angle > M_PI)
    angle -= 2. * M_PI;
  while (angle < -M_PI)
    angle += 2. * M_PI;
  return angle;
}
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values 
   */

  is_initialized_ = false;
    // state vector dimension
  n_x_ = 5;
  // set augmented dimension
  n_aug_ = 7;
  // define spreading parameter
  lambda_ = 3 - n_aug_;
  // Initialize weights_
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // sigma point prediction
   Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  if (!is_initialized_)
  {
    switch (meas_package.sensor_type_)
    {
      //for radar, we need to convert to cartesian coordinates
    case MeasurementPackage::RADAR:
    {
      double rho = meas_package.raw_measurements_[0];     // range
      double phi = meas_package.raw_measurements_[1];     // angle
      double rho_dot = meas_package.raw_measurements_[2]; //velocity in radial direction
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      x_ << x, y, rho_dot, 0, 0; //rho_dot is only an estimation for vel_abs, since the direction is different

      P_ << 0.5, 0, 0, 0, 0,
          0, 0.5, 0, 0, 0,
          0, 0, 10, 0, 0,
          0, 0, 0, 10, 0,
          0, 0, 0, 0, 1;
      break;
    }
    //for laser we can directly take over the measurements
    case MeasurementPackage::LASER:
    {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      P_ << 0.02, 0, 0, 0, 0,
          0, 0.02, 0, 0, 0,
          0, 0, 50, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
      break;
      break;
    }
    default:
      std::cout << "Sensortype not supported" << std::endl;
      return;
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // calculate dt
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Prediction step
  Prediction(delta_t);

  // Update step
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }

  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    UpdateLidar(meas_package);
  }
}

void UKF::Prediction(double delta_t)
{
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_).setZero();

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_).setZero();

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create augmented mean state
  x_aug.head(5) = x_;

  // create augmented covariance matrix
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; ++i)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  // predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > std::numeric_limits<double>::epsilon())
    {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else
    {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;
    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
  // predicted state mean
  x_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // iterate over sigma points
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2. * M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  // set measurement dimension, lidar can measure px, py
  int n_z = 2;
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  // calculate mean predicted measurement
  z_pred.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    z_pred += weights_(i) * Zsig.col(i);
  }
  // calculate innovation covariance matrix S
  S.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  MatrixXd R = MatrixXd(n_z, n_z).setZero();
  R(0, 0) = pow(std_laspx_, 2);
  R(1, 1) = pow(std_laspy_, 2);

  S += R;

  // Update state and covariance

  // create matrix for cross correlation Tc
  MatrixXd Tc(n_x_, n_z);
  // calculate cross correlation matrix
  Tc.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  // update state mean and covariance matrix
  x_ += K * (z_diff);
  P_ -= K * S * K.transpose();

std::ofstream log("NIS_LASER.txt", std::ios_base::app | std::ios_base::out);//creating a new text file
double NIS_laser = z_diff.transpose() * S.inverse() * z_diff;
log << NIS_laser; //add NIS calculated value
log << "\n";
}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    double r = sqrt(px * px + py * py);
    Zsig(0, i) = r;
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px * v1 + py * v2) / r;
  }

  // mean predicted measurement
  z_pred.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    z_pred += weights_(i) * Zsig.col(i);
  }
  // innovation covariance matrix S
  S.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    z_diff(1) = NormalizeAngle(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z).setZero();
  R(0, 0) = std_radr_ * std_radr_;
  R(1, 1) = std_radphi_ * std_radphi_;
  R(2, 2) = std_radrd_ * std_radrd_;

  S += R;

  // Update state and covariance

  // create matrix for cross correlation Tc
  MatrixXd Tc(n_x_, n_z);
  // calculate cross correlation matrix
  Tc.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    z_diff(1) = NormalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    x_diff(3) = NormalizeAngle(x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  // angle normalization
  z_diff(1) = NormalizeAngle(z_diff(1));

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

std::ofstream log("NIS_RADAR.txt", std::ios_base::app | std::ios_base::out);//creating a new text file
double NIS_RADAR = z_diff.transpose() * S.inverse() * z_diff;
log << NIS_RADAR; //add NIS calculated value
log << "\n";
}