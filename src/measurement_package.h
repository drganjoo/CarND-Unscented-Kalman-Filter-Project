#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"

//class MeasurementPackage {
//public:
//  long timestamp_;
//
//  enum SensorType{
//    LASER,
//    RADAR
//  } sensor_type_;
//
//  Eigen::VectorXd raw_measurements_;
//
//};


#ifndef EKF_MEASUREMENTS_H
#define EKF_MEASUREMENTS_H

struct TrackableObjectState
{
    double p_x;
    double p_y;
    double v_x;
    double v_y;
    double yaw;
    double v;
};

struct Measurement
{
    long long timestamp;
};


struct Radar : Measurement
{
    double rho;
    double theta;
    double rhodot;
    
    Eigen::VectorXd GetVector() const {
        Eigen::VectorXd z(3);
        z << rho, theta, rhodot;
        return z;
    }
};

struct Lidar : Measurement
{
    double x;
    double y;
};


#endif //EKF_MEASUREMENTS_H

#endif /* MEASUREMENT_PACKAGE_H_ */
