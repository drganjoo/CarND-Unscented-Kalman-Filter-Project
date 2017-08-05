//
//  ukf_test.cpp
//  UnscentedKF
//
//  Created by Fahad Zubair on 8/5/17.
//
//

#include "ukf_test.hpp"
#include <iostream>
#include <string>
#include "Eigen/Dense"

using namespace std;

void UKFTest::Test2()
{
    ukf.std_a_ = 0.2;
    ukf.std_yawdd_ = 0.2;
    
    ukf.x_ <<   5.7441,
    1.3800,
    2.2049,
    0.5015,
    0.3528;
    
    //create example covariance matrix
    ukf.P_ = MatrixXd(ukf.n_x_, ukf.n_x_);
    ukf.P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
    -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
    0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
    -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
    -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

    VectorXd x_aug = VectorXd(7);
    MatrixXd P_aug = MatrixXd(7, 7);
    
    ukf.GenerateAugStateAndCovariance(&x_aug, &P_aug);
    
    MatrixXd Xsig;
    ukf.GenerateSigmaPoints(x_aug, P_aug, &Xsig);

    istringstream expected(
"5.7441 5.85768   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441 " \
"1.38  1.34566  1.52806     1.38     1.38     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38     1.38     1.38" \
"2.2049  2.28414  2.24557  2.29582   2.2049   2.2049   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049   2.2049   2.2049" \
"0.5015  0.44339 0.631886 0.516923 0.595227   0.5015   0.5015   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015   0.5015   0.5015" \
"0.3528 0.299973 0.462123 0.376339  0.48417 0.418721   0.3528   0.3528 0.405627 0.243477 0.329261  0.22143 0.286879   0.3528   0.3528" \
"0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641        0" \
"0        0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641");
                           
    PrintAndTest(Xsig, expected);
}

void UKFTest::Test1()
{
    ukf.x_ <<
    5.7441,
    1.3800,
    2.2049,
    0.5015,
    0.3528;
    
    //set example covariance matrix
    ukf.P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
    -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
    0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
    -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
    -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;
    
    MatrixXd sigma(ukf.n_x_, 2 * ukf.n_x_ + 1);
    ukf.GenerateSigmaPoints(ukf.x_, ukf.P_, &sigma);
    
    istringstream expected(
                           "5.7441  5.85768   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441 " \
                           "1.38  1.34566  1.52806     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38 " \
                           "2.2049  2.28414  2.24557  2.29582   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049 " \
                           "0.5015  0.44339 0.631886 0.516923 0.595227   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015 " \
                           "0.3528 0.299973 0.462123 0.376339  0.48417 0.418721 0.405627 0.243477 0.329261  0.22143 0.286879 ");
    PrintAndTest(sigma, expected);
}

void UKFTest::PrintAndTest(const Eigen::MatrixXd &sigma, std::istringstream &expected)
{
    ostringstream output;
    output << sigma << endl;
    
    cout << output.str() << endl;
    
    for (int r = 0; r < ukf.n_x_; r++)
        for (int c = 0; c < 2 * ukf.n_x_ + 1; c++) {
            double value;
            expected >> value;
            
            if (fabs(sigma(r,c) - value) > 0.001 ){
                cout << "mismatch at " << r << "," << c << " expected:" << value << " found:" << sigma(r,c) << endl;
            }
        }
}
