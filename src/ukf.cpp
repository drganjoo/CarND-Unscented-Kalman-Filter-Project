#include "ukf.h"
#include "Eigen/Dense"
#include "tools.h"
#include <iostream>
#include <string>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() :
    std_a_(5),
    std_yawdd_(1),
    std_laspx_(0.15),
    std_laspy_(0.15),
    std_radr_(0.3),
    std_radphi_(0.03),
    std_radrd_(0.3),
    n_x_(5),
    n_aug_(7)
{
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
    
    time_us_ = 0;
    
    Xsig_pred_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    weights_ = VectorXd(2 * n_aug_ + 1);
    
    const int lambda = GetLambda(n_aug_);
    
    weights_(0) = lambda / (lambda + n_aug_);
    for (int i = 1; i < weights_.rows(); i++)
        weights_(i) = 1.0 / (2.0 * (lambda + n_aug_));
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const Radar &radar)
{
    if (!is_initialized_) {
        Initialize(radar);
    }
    else {
        double delta_t = radar.timestamp - time_us_;
        
        MatrixXd sigma_radar_space, covaraince_radar_space;
        VectorXd pred_radar_space;
        
        Prediction(delta_t);
        TransformSigmaToRadar(&sigma_radar_space, &pred_radar_space, &covaraince_radar_space);
        
        VectorXd z(pred_radar_space.cols());
        z(0) = radar.rho;
        z(1)  = radar.theta;
        z(2) = radar.rhodot;
        
        Update(z, sigma_radar_space, pred_radar_space, covaraince_radar_space);
        
        time_us_ = radar.timestamp;
    }
}


void UKF::ProcessMeasurement(const Lidar &lidar)
{
    if (!is_initialized_) {
        Initialize(lidar);
    }
    else {
        double delta_t = lidar.timestamp - time_us_;
        
        MatrixXd sigma_lidar_space, covaraince_lidar_space;
        VectorXd pred_lidar_space;

        Prediction(delta_t);
        TransformSigmaToLidar(&sigma_lidar_space, &pred_lidar_space, &covaraince_lidar_space);
        
        VectorXd z(pred_lidar_space.cols());
        z(0) = lidar.x;
        z(1)  = lidar.y;
        
        Update(z, sigma_lidar_space, pred_lidar_space, covaraince_lidar_space);
        
        time_us_ = lidar.timestamp;
    }
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
    VectorXd x_aug = VectorXd(7);
    MatrixXd P_aug = MatrixXd(7, 7);

    GenerateAugStateAndCovariance(&x_aug, &P_aug);
    
    //cout << "P_aug" << endl;
    //cout << P_aug << endl;
    
    // create augmented sigma points
    MatrixXd Xsig_aug(n_aug_, 2 * n_aug_ + 1 );
    GenerateSigmaPoints(x_aug, P_aug, &Xsig_aug);
    
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
            Xsig_pred_(0, i) = px + v_psi_k * (sin(psi + psi_dot * delta_t) - sin(psi)) + (0.5 * delta_t2 * cos(psi) * mu_a);
            Xsig_pred_(1, i) = py + v_psi_k * (-cos(psi + psi_dot * delta_t) + cos(psi)) + (0.5 * delta_t2 * sin(psi) * mu_a);
        }
        else {
            Xsig_pred_(0, i) = px + v * cos(psi) * delta_t + (0.5 * delta_t2 * cos(psi) * mu_a);
            Xsig_pred_(1, i) = py + v * sin(psi) * delta_t + (0.5 * delta_t2 * sin(psi) * mu_a);
        }
        
        Xsig_pred_(2,i) = v + 0 + delta_t * mu_a;
        Xsig_pred_(3,i) = psi + psi_dot * delta_t + 0.5 * delta_t2 * mu_yaw;
        Xsig_pred_(4,i) = psi_dot + 0 + delta_t * mu_yaw;
    }
    
    // predicted mean
    // x(k+1|k) = Σ(i=1 to n_sgima) wi * X_sigma
    x_ = Xsig_pred_ * weights_;    // this will multiply wi with Xsig_pred and sum it as well
    
    //cout << "x" << endl << x << endl;
    
    // predict state covariance matrix
    // P =∑ wi (X − x)(X − x)T
    // The following is equivalent in Eigen. Rather than doing X-Y we can do -X + Y
    // and for column wise multiplication of two matrices (X * Y) has to be done
    // X * Y.asDiagonal
    
    MatrixXd diff = (-Xsig_pred_).colwise() + x_;
    P_ = diff * weights_.asDiagonal() * diff.transpose();
}

TrackableObjectState UKF::GetState()
{
    TrackableObjectState state;
    
    state.p_x = x_(PX_INDEX);
    state.p_y = x_(PY_INDEX);
    state.v  = x_(V_INDEX);
    state.yaw = x_(YAW_INDEX);
    
    state.v_x = cos(state.yaw) * state.v;
    state.v_y = sin(state.yaw) * state.v;
    
    return state;
}

void UKF::Initialize(const Radar &radar)
{
    const double x = radar.rho * cos(radar.theta);
    const double y = radar.rho * sin(radar.theta);
    x_.fill(0);
    x_(PX_INDEX) = x;
    x_(PY_INDEX) = y;
    time_us_ = radar.timestamp;
}

void UKF::Initialize(const Lidar &lidar)
{
    x_.fill(0);
    x_(PX_INDEX) = lidar.x;
    x_(PY_INDEX) = lidar.y;
    time_us_ = lidar.timestamp;
}

void UKF::TransformSigmaToRadar(MatrixXd *sigma_radar_space, VectorXd *pred_radar_space,
                                MatrixXd *covariance_radar_space)
{
    const int n_z = 3;
    
    // create matrix for sigma points in measurement space
    *sigma_radar_space = MatrixXd(n_z, 2 * n_aug_ + 1);
    *pred_radar_space = VectorXd(n_z);
    *covariance_radar_space = MatrixXd(n_z, n_z);
    
    //transform sigma points into measurement space
    for (int i = 0; i < Xsig_pred_.cols(); i++) {
        const double px = Xsig_pred_(0, i);
        const double py = Xsig_pred_(1, i);
        const double v = Xsig_pred_(2, i);
        const double psi = Xsig_pred_(3, i);
        
        (*sigma_radar_space)(0, i) = sqrt(px * px + py * py);
        (*sigma_radar_space)(1, i) = atan2(py, px);
        (*sigma_radar_space)(2, i) = (px * cos(psi) * v + py * sin(psi) * v) / sqrt(px * px + py * py);
    }
    
    //cout << "Zsig" << endl << (*sigma_radar_space) << endl;
    
    // Mean predicted measurement would be the weighted sum of each row
    // In Eigen, instead of writing a loop to do this, we can just
    // multiply Zsig * weights to do a row wise sum
    
    *pred_radar_space = *sigma_radar_space * weights_;
    
    //cout << "Z pred" << endl << (*pred_radar_space) << endl;
    
    // Measurement covariance matrix S:
    // S = weights * (Z - z) * (Z - z)transpose + R
    MatrixXd R(n_z, n_z);
    R.fill(0);
    R(0, 0) = std_radr_ * std_radr_;
    R(1, 1) = std_radphi_ * std_radphi_;
    R(2, 2) = std_radrd_ * std_radrd_;
    
    MatrixXd diff = (-(*sigma_radar_space)).colwise() + *pred_radar_space;
    
    // normalize angles
    for (int i = 0; i < diff.cols(); i++)
        diff(1,i) = Tools::NormalizeAngleRad(diff(1,i));
    
    // S = weights * (Z - z) * (Z - z)transpose + R
    *covariance_radar_space = diff * weights_.asDiagonal() * diff.transpose() + R;
}


void UKF::TransformSigmaToLidar(MatrixXd *sigma_lidar_space, VectorXd *pred_lidar_space,
                                MatrixXd *covariance_lidar_space)
{
    const int n_z = 2;
    
    // create matrix for sigma points in measurement space
    *sigma_lidar_space = MatrixXd(n_z, 2 * n_aug_ + 1);
    *pred_lidar_space = VectorXd(n_z);
    *covariance_lidar_space = MatrixXd(n_z, n_z);
    
    // copy over sigma points for px and py already computed
    (*sigma_lidar_space) = Xsig_pred_.topLeftCorner(n_z, sigma_lidar_space->cols());
    
    //cout << "Zsig" << endl << (*sigma_lidar_space) << endl;
    
    // Mean predicted measurement would simply be the px and py already computed in x_
    *pred_lidar_space = x_.topRows(n_z);
    
    //cout << "Z pred" << endl << (*pred_lidar_space) << endl;
    
    // Measurement covariance matrix S:
    // S = weights * (Z - z) * (Z - z)transpose + R
    MatrixXd R(n_z, n_z);
    R.fill(0);
    R(0, 0) = std_laspx_ * std_laspx_;
    R(1, 1) = std_laspy_ * std_laspy_;
    
    MatrixXd diff = (-(*sigma_lidar_space)).colwise() + *pred_lidar_space;
    
    // S = weights * (Z - z) * (Z - z)transpose + R
    *covariance_lidar_space = diff * weights_.asDiagonal() * diff.transpose() + R;
}


void UKF::Update(const VectorXd &z, const MatrixXd &Zsig, const VectorXd &z_pred, const MatrixXd &S)
{
    // T is the Cross Corelation Matrix that maps from the
    // 3 readings of radar to the 7 augmented state matrix that we have
    MatrixXd Tc = MatrixXd(n_x_, z_pred.cols());

    // T = w * (X - x)(Z - z)^T
    MatrixXd diff_x = (-Xsig_pred_).colwise() + x_;
    MatrixXd diff_z = (-Zsig).colwise() + z_pred;   // Zk - z

    // normalize angles
    for (int i = 0; i < diff_z.cols(); i++) {
        diff_z(1, i) = Tools::NormalizeAngleRad(diff_z(1, i));
    }
    
    Tc = diff_x * weights_.asDiagonal() * diff_z.transpose();
    //cout << "Tc" << endl << Tc << endl;
    
    //calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //cout << "K" << endl << K << endl;
    
    //update state mean and covariance matrix
    
    MatrixXd diff_measure = z - z_pred;
    // normalize angles
    for (int i = 0; i < diff_z.cols(); i++) {
        diff_z(1, i) = Tools::NormalizeAngleRad(diff_z(1, i));
    }
    
    x_ = x_ + K * diff_measure;
    P_ = P_ - K * S * K.transpose();
}

void UKF::GenerateSigmaPoints(const VectorXd &x, const MatrixXd &P, MatrixXd *sigma_points)
{
    MatrixXd p_root = P.llt().matrixL();
    const auto n_pts = x.rows();
    
    cout << "n_x" << n_pts << endl;
    cout << p_root << endl;

    const int lambda = GetLambda(n_pts);

    sigma_points->col(0) = x;
    
    for (int i = 0; i < n_pts; i++) {
        const auto sqroot = sqrt(lambda + n_pts)  * p_root.col(i);
        
        sigma_points->col(i + 1) = x + sqroot;
        sigma_points->col(i + 1 + n_pts) = x - sqroot;
    }
}


void UKF::GenerateAugStateAndCovariance(VectorXd *x_aug, MatrixXd *P_aug)
{
    //create augmented mean state
    x_aug->head(n_x_) = x_;
    (*x_aug)(n_x_) = 0;              // 0 mean noise
    (*x_aug)(n_x_ + 1) = 0;          // 0 mean noise
    
    // create augmented covariance matrix
    P_aug->topLeftCorner(n_x_, n_x_) = P_;
    (*P_aug)(n_x_, n_x_) = std_a_ * std_a_;
    (*P_aug)(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;
}
