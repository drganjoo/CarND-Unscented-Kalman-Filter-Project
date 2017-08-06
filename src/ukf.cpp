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
    use_laser_ = true;
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
    
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_.fill(0);
    weights_ = VectorXd(2 * n_aug_ + 1);
    
    const int lambda = GetLambda(n_aug_);
    
    weights_(0) = (double)lambda / (double)(lambda + n_aug_);
    for (int i = 1; i < weights_.rows(); i++)
        weights_(i) = 1.0 / (2.0 * (lambda + n_aug_));
    
    // sum of weights has to be 1
    assert(weights_.sum() - 1 < 0.001);
    
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << std_laspx_ * std_laspx_, 0,
                0, std_laspy_ * std_laspy_;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(const Radar &radar)
{
    if (!is_initialized_) {
        Initialize(radar);
    }
    else {
        if (!use_radar_)
            return;

        PredictState(&radar);
        
        MatrixXd radar_sig, radar_S;
        VectorXd radar_pred;
        
        TransformSigmaToRadar(&radar_sig, &radar_pred, &radar_S);
        
        VectorXd z = radar.GetVector();
        Update(z, radar_sig, radar_pred, radar_S);
    }
}


void UKF::ProcessMeasurement(const Lidar &lidar)
{
    if (!is_initialized_) {
        Initialize(lidar);
    }
    else {
        if (!use_laser_)
            return;

        PredictState(&lidar);

        VectorXd z(2);
        z << lidar.x, lidar.y;

        MatrixXd H = MatrixXd(2, n_x_);
        H << 1, 0, 0, 0, 0,
                    0, 1, 0, 0, 0;

        VectorXd y = z - H * x_;
        MatrixXd Ht = H.transpose();
        MatrixXd S = H * P_ * Ht + R_laser_;
        MatrixXd Si = S.inverse();
        MatrixXd K = P_ * Ht * Si;
        MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
        
        x_ = x_ + (K * y);
        P_ = (I - K * H) * P_;
        
//        MatrixXd lidar_sig, lidar_S;
//        VectorXd lidar_pred;
//
//        TransformSigmaToLidar(&lidar_sig, &lidar_pred, &lidar_S);
//        
//
//        Update(z, lidar_sig, lidar_pred, lidar_S, false);
    }
}


void UKF::PredictState(const Measurement *measurement)
{
    double delta_t_secs = (measurement->timestamp - time_us_) / 1000000.0;
    time_us_ = measurement->timestamp;
    
    assert(fabs(delta_t_secs - 0.05) < 0.001);

    VectorXd x_aug;
    MatrixXd P_aug;
    
    GenerateAugStateAndCovariance(&x_aug, &P_aug);
    
    //cout << "P_aug" << endl;
    //cout << P_aug << endl;
    
    MatrixXd Xsig_aug = GenerateSigmaPoints(x_aug, P_aug);
    PredictSigmaPoints(Xsig_aug, delta_t_secs);
    PredictMeanAndCovariance();
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
    is_initialized_ = true;
}

void UKF::Initialize(const Lidar &lidar)
{
    x_.fill(0);
    x_(PX_INDEX) = lidar.x;
    x_(PY_INDEX) = lidar.y;
    time_us_ = lidar.timestamp;
    is_initialized_ = true;
}

void UKF::TransformSigmaToRadar(MatrixXd *Zsig, VectorXd *z_pred,
                                MatrixXd *S)
{
    const int n_z = 3;
    
    // create matrix for sigma points in measurement space
    *Zsig = MatrixXd(n_z, Xsig_pred_.cols());
    *z_pred = VectorXd(n_z);
    *S = MatrixXd(n_z, n_z);
    
    //transform sigma points into measurement space
    for (int i = 0; i < Xsig_pred_.cols(); i++) {
        const double px = Xsig_pred_(0, i);
        const double py = Xsig_pred_(1, i);
        const double v = Xsig_pred_(2, i);
        const double psi = Xsig_pred_(3, i);
        
        (*Zsig)(0, i) = sqrt(px * px + py * py);
        (*Zsig)(1, i) = atan2(py, px);
        (*Zsig)(2, i) = (px * cos(psi) * v + py * sin(psi) * v) / sqrt(px * px + py * py);
    }
    
    //cout << "Zsig" << endl << (*sigma_radar_space) << endl;
    
    // Mean predicted measurement would be the weighted sum of each row
    // In Eigen, instead of writing a loop to do this, we can just
    // multiply Zsig * weights to do a row wise sum
    
    *z_pred = *Zsig * weights_;
    
    //cout << "Z pred" << endl << (*pred_radar_space) << endl;
    
    // Measurement covariance matrix S:
    // S = weights * (Z - z) * (Z - z)transpose + R
    MatrixXd R(n_z, n_z);
    R.fill(0);
    R(0, 0) = std_radr_ * std_radr_;
    R(1, 1) = std_radphi_ * std_radphi_;
    R(2, 2) = std_radrd_ * std_radrd_;
    
    MatrixXd diff = (-(*Zsig)).colwise() + *z_pred;
    
    // normalize angles
    for (int i = 0; i < diff.cols(); i++)
        diff(1,i) = Tools::NormalizeAngleRad(diff(1,i));
    
    // S = weights * (Z - z) * (Z - z)transpose + R
    *S = diff * weights_.asDiagonal() * diff.transpose() + R;
}


void UKF::TransformSigmaToLidar(MatrixXd *Zsig, VectorXd *z_pred,
                                MatrixXd *S)
{
    const int n_z = 2;
    
    // create matrix for sigma points in measurement space
    *Zsig = MatrixXd(n_z, Xsig_pred_.cols());
    *z_pred = VectorXd(n_z);
    *S = MatrixXd(n_z, n_z);
    
    // copy over sigma points for px and py already computed
    (*Zsig) = Xsig_pred_.topLeftCorner(n_z, Zsig->cols());
    
    // Mean predicted measurement would simply be the px and py already computed in x_
    *z_pred = x_.topRows(n_z);
    
    // Measurement covariance matrix S:
    // S = weights * (Z - z) * (Z - z)transpose + R
    MatrixXd R(n_z, n_z);
    R.fill(0);
    R(0, 0) = std_laspx_ * std_laspx_;
    R(1, 1) = std_laspy_ * std_laspy_;
    
    MatrixXd diff = (-(*Zsig)).colwise() + *z_pred;
    
    // S = weights * (Z - z) * (Z - z)transpose + R
    *S = diff * weights_.asDiagonal() * diff.transpose() + R;
}


void UKF::Update(const VectorXd &z, const MatrixXd &Zsig, const VectorXd &z_pred, const MatrixXd &S, bool normalize /*= true*/)
{
    MatrixXd diff_x = (-Xsig_pred_).colwise() + x_;
    MatrixXd diff_z = (-Zsig).colwise() + z_pred;   // Zk - z
    
    //cout << "diff_x" << endl << diff_x << endl;
    //cout << "diff_z" << endl << diff_z << endl;
    
    if (normalize) {
        for (int i = 0; i < diff_z.cols(); i++) {
            diff_z(1, i) = Tools::NormalizeAngleRad(diff_z(1, i));
        }
    }
    
    // T is the Cross Corelation Matrix that maps from the
    // 3 readings of radar to the 7 augmented state matrix that we have
    // T = w * (X - x)(Z - z)^T
    MatrixXd Tc = diff_x * weights_.asDiagonal() * diff_z.transpose();
    //cout << "Tc" << endl << Tc << endl;
    
    //calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //cout << "K" << endl << K << endl;
    
    //update state mean and covariance matrix
    VectorXd diff_measure = z - z_pred;
    
    if (normalize) {
        for (int i = 0; i < diff_z.cols(); i++) {
            diff_z(1, i) = Tools::NormalizeAngleRad(diff_z(1, i));
        }
    }
    
    x_ = x_ + K * diff_measure;
    P_ = P_ - K * S * K.transpose();
    
    ComputeNIS(z, z_pred, S);
}

MatrixXd UKF::GenerateSigmaPoints(const VectorXd &x, const MatrixXd &P)
{
    const auto n_pts = x.rows();
    MatrixXd sigma_points(n_pts, 2 * n_pts + 1);
    MatrixXd p_root = P.llt().matrixL();
    
    //cout << "n_x" << n_pts << endl;
    //cout << p_root << endl;
    
    const int lambda = GetLambda(n_pts);
    
    sigma_points.fill(0);
    sigma_points.col(0) = x;
    
    for (int i = 0; i < n_pts; i++) {
        const auto sqroot = sqrt(lambda + n_pts)  * p_root.col(i);
        
        sigma_points.col(i + 1) = x + sqroot;
        sigma_points.col(i + 1 + n_pts) = x - sqroot;
    }
    
    return sigma_points;
}


void UKF::GenerateAugStateAndCovariance(VectorXd *x_aug, MatrixXd *P_aug)
{
    *x_aug = VectorXd(n_aug_);
    *P_aug = MatrixXd(n_aug_, n_aug_);
    
    x_aug->fill(0);
    P_aug->fill(0);
    
    //create augmented mean state
    x_aug->head(n_x_) = x_;
    (*x_aug)(n_x_) = 0;              // 0 mean noise
    (*x_aug)(n_x_ + 1) = 0;          // 0 mean noise
    
    // create augmented covariance matrix
    P_aug->topLeftCorner(n_x_, n_x_) = P_;
    (*P_aug)(n_x_, n_x_) = std_a_ * std_a_;
    (*P_aug)(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;
}

void UKF::PredictSigmaPoints(const MatrixXd &Xsig_aug, const double delta_t)
{
    const double delta_t2 = delta_t * delta_t; // time diff in micro sec
    
    assert(Xsig_pred_.rows() == n_x_);
    assert(Xsig_pred_.cols() == 2 * n_aug_ + 1);
    
    // predict sigma points
    Xsig_pred_.fill(0);
    
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
}

void UKF::PredictMeanAndCovariance()
{
    //cout << "Weights: " << endl << weights_ << endl;
    // predicted mean / covariance
    // x = Σ wi * X_sigma
    // P =∑ wi (X − x)(X − x)T

    x_ = Xsig_pred_ * weights_;    // this will multiply wi with Xsig_pred and sum it as well
    
    MatrixXd diff = (-Xsig_pred_).colwise() + x_;
    P_ = diff * weights_.asDiagonal() * diff.transpose();
}

void UKF::ComputeNIS(const VectorXd &z, const VectorXd &z_pred, const MatrixXd &S)
{
    VectorXd diff = z - z_pred;
    nis_ = diff.transpose() * S.inverse() * diff;
}
