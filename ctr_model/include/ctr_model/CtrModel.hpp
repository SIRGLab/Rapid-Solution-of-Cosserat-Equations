#pragma once

#include "ctr_model/Calculation.hpp"
#include <ctr_model/MyParamsConfig.h>
#include <dynamic_reconfigure/server.h>

#include <ros/ros.h>
//#include <string>
#include "sensor_msgs/JointState.h"
#include "std_msgs/String.h"
#include <iostream>
#include <math.h>
#include <std_msgs/Bool.h>
#include <std_msgs/Float32.h>

#include "alglib/optimization.h"
#include "alglib/stdafx.h"
#include "tf/transform_broadcaster.h"
#include "tf2_ros/static_transform_broadcaster.h"
#include <algorithm>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <visualization_msgs/Marker.h>

using namespace alglib;
#include <complex>
using namespace std;
using namespace boost::numeric::odeint;
using namespace std::chrono;
using namespace Eigen;

namespace ctr_model {
class CtrModel {
public:
  CtrModel(ros::NodeHandle &nodeHandle);
  virtual ~CtrModel();

private:
  ros::NodeHandle nodeHandle_;
  void chatterCallback(const ros::MessageEvent<tf2_msgs::TFMessage const> &msg);
  static void callback(ctr_model::MyParamsConfig &config, uint32_t level);
  Calculation calculation;
  Calculation real_model;
  void UpdateRviz();
  void ControlLoop();
  void SetInitialValues();
  void SaveTipPositionToFile();

  real_1d_array u_0;
  real_2d_array P;
  // input vectors
  real_1d_array q;
  real_1d_array f;
  real_1d_array x0;
  // control loop variables
  real_2d_array A;
  real_2d_array Q;
  real_2d_array R;
  real_2d_array Gamma;
  real_2d_array GammaT;

  vector<double> tip_pos[3];

  static double force[3];
  static double dt;
  static double total_time;
  static double Q_param;
  static double R_param;
  static double V;
  static bool start_simulation;
  static bool simulation_done;
};

} // namespace ctr_model
