#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
    
    UKF();
    virtual ~UKF();
    
    void ProcessMeasurement(const Radar &radar);
    void ProcessMeasurement(const Lidar &lidar);
    
    
    TrackableObjectState GetState();
    
    friend class UKFTest;

private:
    inline double GetLambda(const int pts){ return 3 - pts; }
    
private:
    void PredictState(const Measurement *measurement);
    
    void Initialize(const Radar &radar);
    void Initialize(const Lidar &lidar);

    MatrixXd GenerateSigmaPoints(const VectorXd &x, const MatrixXd &P);
    void GenerateAugStateAndCovariance(VectorXd *x_aug, MatrixXd *P_aug);
    void PredictSigmaPoints(const MatrixXd &Xsig_aug, const double delta_t);
    void PredictMeanAndCovariance();
    
    void TransformSigmaToRadar(MatrixXd *sigma_radar_space, VectorXd *pred_radar_space,
                               MatrixXd *covariance_radar_space);
    void TransformSigmaToLidar(MatrixXd *sigma_lidar_space, VectorXd *pred_lidar_space,
                               MatrixXd *covariance_lidar_space);
    void Update(const VectorXd &z, const MatrixXd &Zsig, const VectorXd &z_pred, const MatrixXd &S);
    
    //void Update(const Lidar &lidar);
    
private:
    const int PX_INDEX = 0;
    const int PY_INDEX = 1;
    const int V_INDEX = 2;
    const int YAW_INDEX = 3;
    
    bool is_initialized_;
    bool use_laser_;     // ignore lidar if set
    bool use_radar_;     // ignore radar if set
    
    VectorXd x_;         // state: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    MatrixXd P_;         // state covariance matrix
    MatrixXd Xsig_pred_; // predicted sigma points matrix
    VectorXd weights_;   // Weights of sigma points
    
    long long time_us_;  // last reading time in microseconds
    
    double std_a_;       // Process noise standard deviation longitudinal acceleration in m/s^2
    double std_yawdd_;   // Process noise standard deviation yaw acceleration in rad/s^2
    const double std_laspx_;   // Laser measurement noise standard deviation position1 in m
    const double std_laspy_;   // Laser measurement noise standard deviation position2 in m
    double std_radr_;    // Radar measurement noise standard deviation radius in m
    double std_radphi_;  // Radar measurement noise standard deviation angle in rad
    double std_radrd_ ;  // Radar measurement noise standard deviation radius change in m/s
    
    const int n_x_;      // State dimension
    const int n_aug_;          // Augmented state dimension
};

#endif /* UKF_H */
