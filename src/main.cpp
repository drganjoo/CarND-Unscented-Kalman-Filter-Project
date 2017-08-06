#ifndef _WIN32

#include <unistd.h>
#include <uWS/uWS.h>

#else

// disable warnings in uWS/uWS code
#pragma warning(push)
#pragma warning(disable:4251)
#pragma warning(disable:4800)
#pragma warning(disable:4996)

#include <uWS/uWS.h>

#pragma warning(pop)

#endif

#include <math.h>
#include <iostream>
#include "json.hpp"
#include "ukf.h"
#include "tools.h"
#include "ukf_test.hpp"
#include <fstream>

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

vector<VectorXd> estimations;
vector<VectorXd> ground_truth;
UKF ukf;
ofstream nis_file;

const char *kVersion = "1.0";

Radar ParseRadar(istream &tokens) {
    Radar r;
    
    // sensor_type, rho_measured, theta_measured, rhodot_measured, timestamp, x_groundtruth, y_groundtruth, vx_groundtruth, vy_groundtruth, yaw_groundtruth, yawrate_groundtruth.
    tokens >> r.rho >> r.theta >> r.rhodot >> r.timestamp;
    return r;
}

Lidar ParseLaser(istream &tokens) {
    Lidar l;
    
    // sensor_type, x_measured, y_measured, timestamp, x_groundtruth, y_groundtruth, vx_groundtruth, vy_groundtruth, yaw_groundtruth, yawrate_groundtruth.
    tokens >> l.x >> l.y >> l.timestamp;
    return l;
}


bool GetMeasurementLine(uWS::WebSocket<uWS::SERVER> &ws, char *data, size_t length, string *measurement_line) {
    bool found = false;
    
    // messages coming in from the simulator should have 42 in the beginning
    if (length > 2 && data[0] == '4' && data[1] == '2') {
        data[length] = 0;
        
        if (strstr(data, "null")) {
            // send response back to the simulator as it does not have a valid JSON object in it
            
            const char msg[] = "42[\"manual\",{}]";
            constexpr size_t length = sizeof(msg) - 1;
            
            ws.send(msg, length, uWS::OpCode::TEXT);
        }
        else {
            // look for json message inside [ and ]
            char *start_bracket = strchr(data, '[');
            if (start_bracket) {
                auto end_bracket = strrchr(data, ']');
                if (end_bracket) {
                    *(end_bracket + 1) = 0;
                    
                    // convert  message into a JSon object, get the measurement line in case
                    // it is a telemetry message
                    auto j = json::parse(start_bracket);
                    auto wsEvent = j[0].get<std::string>();
                    
                    if (wsEvent == "telemetry") {
                        *measurement_line = j[1]["sensor_measurement"];
                        found = true;
                    }
                }
                else {
                    // looks like a bad message, it doesn't have a ]
                    assert(0);
                }
            }
        }
    }
    
    return found;
}


void SendEstimates(uWS::WebSocket<uWS::SERVER> &ws, TrackableObjectState &state, VectorXd &rmse) {
    json msgJson;
    msgJson["estimate_x"] = state.p_x;
    msgJson["estimate_y"] = state.p_y;
    msgJson["rmse_x"] = rmse(0);
    msgJson["rmse_y"] = rmse(1);
    msgJson["rmse_vx"] = rmse(2);
    msgJson["rmse_vy"] = rmse(3);
    
    auto msg = "42[\"estimate_marker\"," + msgJson.dump() + "]";
    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
}


void ProcessMeasurement(uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    string measurement;
    
    if (GetMeasurementLine(ws, data, length, &measurement)) {
        //cout << measurement << endl;
        istringstream iss(measurement);
        
        string type;
        iss >> type;
        
        if (type[0] == 'R') {
            Radar r = ParseRadar(iss);
            ukf.ProcessMeasurement(r);
        }
        else if (type[0] == 'L') {
            Lidar l = ParseLaser(iss);
            ukf.ProcessMeasurement(l);
        }
        else {
            assert(0);
        }
        
        // get ground truth out of the incoming stream, add it to the vector of all ground truths
        // and then send the RMSE back to the simulator
        double x, y, v_x, v_y;
        iss >> x >> y >> v_x >> v_y;
        
        ground_truth.push_back(Vector4d(x, y, v_x, v_y));
        
        TrackableObjectState state = ukf.GetState();
        Vector4d estimates;
        estimates << state.p_x , state.p_y, state.v_x, state.v_y;
        
        estimations.push_back(estimates);
        
        // send estimates back to the simulator
        VectorXd rmse = Tools::CalculateRMSE(estimations, ground_truth);
        SendEstimates(ws, state, rmse);
        
        nis_file << ukf.GetNis() << endl;
    }
}

int main()
{
//    UKFTest test;
//    test.Test4();

    cout << "Version: " << kVersion << endl;
    
    nis_file.open("nis_file.txt");
    if (!nis_file.is_open()) {
        cout << "could not open nis_file.txt";
    }
    
    uWS::Hub h;
    
    h.onMessage(ProcessMeasurement);
    
    h.onConnection([](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });
    
    h.onDisconnection([](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
        std::cout << "Disconnected" << std::endl;
    });
    
    int port = 4567;
    if (h.listen("127.0.0.1", port)) {
        cout << "Listening on: " << port << endl;
        h.run();
        return 0;
    }
    else {
        cout << "Could not start listening";
        return -1;
    }
}
