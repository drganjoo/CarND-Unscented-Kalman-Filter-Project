//
//  ukf_test.hpp
//  UnscentedKF
//
//  Created by Fahad Zubair on 8/5/17.
//
//

#ifndef ukf_test_hpp
#define ukf_test_hpp

#include <iostream>
#include <string>
#include "Eigen/Dense"
#include "ukf.h"

class UKFTest
{
public:
    void Test1();
    void Test2();

private:
    void PrintAndTest(const Eigen::MatrixXd &sigma, std::istringstream &expected);
    
    UKF ukf;
};

#endif /* ukf_test_hpp */
