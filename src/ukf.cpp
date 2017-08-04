#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0);
  P_(0,0) = 1.0;
  P_(1,1) = 1.0;
  P_(2,2) = 1.0;
  P_(3,3) = 1.0;
  P_(4,4) = 1.0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

  // Lidar measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Lidar measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  n_x_ = 5;
  n_aug_ = 7;           // 2 more points for process noise
  lambda_ = 3 - n_aug_;
  time_us_ = 0;

  Xsig_pred_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  weights_ = VectorXd(2 * n_aug_ + 1);

  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.rows(); i++)
    weights_(i) = 1.0 / (2.0 * (lambda_ + n_aug_));
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const Radar &radar) {
  if (!is_initialized_) {
    Initialize(radar);
  }
  else {
    double delta_t = radar.timestamp - time_us_;

    Prediction(delta_t);
	TransformSigmaToRadar();
    Update(radar);

    time_us_ = radar.timestamp;
  }
}


void UKF::ProcessMeasurement(const Lidar &lidar) {
  if (!is_initialized_) {
    Initialize(lidar);
  }
  else {
    double delta_t = lidar.timestamp - time_us_;

    Prediction(delta_t);
    Update(lidar);

    time_us_ = lidar.timestamp;
  }
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;              // 0 mean noise
  x_aug(n_x_ + 1) = 0;          // 0 mean noise

  // create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  //cout << "P_aug" << endl;
  //cout << P_aug << endl;

  //create square root matrix
  MatrixXd p_root = P_aug.llt().matrixL();

  //std::cout << "p_root = " << std::endl << p_root << std::endl;

  // create augmented sigma points
  MatrixXd Xsig_aug(n_aug_, 2 * n_aug_ + 1 );
  Xsig_aug.col(0) = x_aug;

  //cout << "X" << endl << x_aug;
  //cout << "Xsig_aug" << endl << Xsig_aug.col(0) << endl;

  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * p_root.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * p_root.col(i);
  }

  const double delta_t2 = delta_t * delta_t; //time diff in sec

  // predict sigma points

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    const double px = Xsig_aug(0, i);
    const double py = Xsig_aug(1, i);
    const double v = Xsig_aug(2, i);
    const double psi = Xsig_aug(3, i);
    const double psi_dot = Xsig_aug(4, i);
    const double mu_a = Xsig_aug(5, i);
    const double mu_yaw = Xsig_aug(6, i);

    if (psi_dot != 0) {
      const double v_psi_k = v / psi_dot;
      Xsig_pred_(0, i) = Xsig_aug(0, i) + v_psi_k * (sin(psi + psi_dot * delta_t) - sin(psi)) + (0.5 * delta_t2 * cos(psi) * mu_a);
      Xsig_pred_(1, i) = Xsig_aug(1, i) + v_psi_k * (-cos(psi + psi_dot * delta_t) + cos(psi)) + (0.5 * delta_t2 * sin(psi) * mu_a);
    }
    else {
      Xsig_pred_(0, i) = Xsig_aug(0, i) + v * cos(psi) * delta_t + (0.5 * delta_t2 * cos(psi) * mu_a);
      Xsig_pred_(1, i) = Xsig_aug(1, i) + v * sin(psi) * delta_t + (0.5 * delta_t2 * sin(psi) * mu_a);
    }

    Xsig_pred_(2,i) = Xsig_aug(2, i) + 0 + delta_t * mu_a;
    Xsig_pred_(3,i) = Xsig_aug(3, i) + psi_dot * delta_t + 0.5 * delta_t2 * mu_yaw;
    Xsig_pred_(4,i) = Xsig_aug(4, i) + 0 + delta_t * mu_yaw;
  }

  // predicted mean
  // x(k+1|k) = Σ(i=1 to n_sgima) wi * X_sigma
  x_ = Xsig_pred_ * weights_;    // this will multiply wi with Xsig_pred and sum it as well

  //cout << "x" << endl << x << endl;

  // predict state covariance matrix
  // P =∑ wi (X − x)(X − x)T
  // The following is equivalent in Eigen. Rather than doing X-Y we can do -X + Y
  MatrixXd diff = (-Xsig_pred_).colwise() + x_;
  P_ = diff * weights_.asDiagonal() * diff.transpose();
}

TrackableObjectState UKF::GetState() {
  TrackableObjectState state;

  state.p_x = x_(PX_INDEX);
  state.p_y = x_(PY_INDEX);
  state.v  = x_(V_INDEX);
  state.yaw = x_(YAW_INDEX);

  state.v_x = cos(state.yaw) * state.v;
  state.v_y = sin(state.yaw) * state.v;

  return state;
}

void UKF::Update(const Radar &radar) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

void UKF::Update(const Lidar &lidar) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

void UKF::Initialize(const Radar &radar) {
  const double x = radar.rho * cos(radar.theta);
  const double y = radar.rho * sin(radar.theta);
  x_.fill(0);
  x_(PX_INDEX) = x;
  x_(PY_INDEX) = y;
  time_us_ = radar.timestamp;
}

void UKF::Initialize(const Lidar &lidar) {
  x_.fill(0);
  x_(PX_INDEX) = lidar.x;
  x_(PY_INDEX) = lidar.y;
  time_us_ = lidar.timestamp;
}

void UKF::TransformSigmaToRadar()
{
	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	//transform sigma points into measurement space
	for (int i = 0; i < Xsig_pred.cols(); i++) {
		const double px = Xsig_pred(0, i);
		const double py = Xsig_pred(1, i);
		const double v = Xsig_pred(2, i);
		const double psi = Xsig_pred(3, i);

		Zsig(0, i) = sqrt(px * px + py * py);
		Zsig(1, i) = atan2(py, px);
		Zsig(2, i) = (px * cos(psi) * v + py * sin(psi) * v) / sqrt(px * px + py * py);
	}

	//cout << "Zsig" << endl << Zsig << endl;

	//calculate mean predicted measurement
	z_pred = Zsig * weights;
	//cout << "Z pred" << endl << z_pred << endl;

	//calculate measurement covariance matrix S
	MatrixXd R(n_z, n_z);
	R.fill(0);
	R(0, 0) = std_radr * std_radr;
	R(1, 1) = std_radphi * std_radphi;
	R(2, 2) = std_radrd * std_radrd;

	MatrixXd diff = (-Zsig).colwise() + z_pred;
	// TODO: normalize angles here
	S = diff * weights.asDiagonal() * diff.transpose() + R;

}
