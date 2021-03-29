#include "ctr_model/CtrModel.hpp"

namespace ctr_model {

ros::Publisher vis_pub;

double CtrModel::force[3];
double CtrModel::dt;
double CtrModel::total_time;
double CtrModel::Q_param;
double CtrModel::R_param;
double CtrModel::V;
bool CtrModel::start_simulation;
bool CtrModel::simulation_done;

CtrModel::CtrModel(ros::NodeHandle &nodeHandle) : nodeHandle_(nodeHandle) {

  dynamic_reconfigure::Server<ctr_model::MyParamsConfig> server;
  dynamic_reconfigure::Server<ctr_model::MyParamsConfig>::CallbackType f;
  f = boost::bind(&callback, _1, _2);
  server.setCallback(f);
  MyParamsConfig config;

  vis_pub = nodeHandle.advertise<visualization_msgs::Marker>(
      "visualization_marker", 0);
  ros::Subscriber sub =
      nodeHandle.subscribe("/tf", 10, &CtrModel::chatterCallback, this);
  ros::spin();
  start_simulation = false;
  simulation_done = false;
}

void CtrModel::chatterCallback(
    const ros::MessageEvent<tf2_msgs::TFMessage const> &msg) {
  ros::Rate r(1000);
  ros::NodeHandle n;
  tf::TransformBroadcaster broadcaster;

  if (n.ok()) {
    if (start_simulation && !simulation_done) {
      ControlLoop();
    }
    if (!start_simulation) {
      simulation_done = false;
    }
  }
}

void CtrModel::SetInitialValues() {

  u_0.setlength(3);  // u_0 = 0
  P.setlength(3, 3); // P=I
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j)
        P[i][j] = 1;
      else
        P[i][j] = 0;
    }
    u_0[i] = 0;
  }

  // ------------------------------------------------------------------------------------------
  // Initialize input vectors
  // ------------------------------------------------------------------------------------------
  q.setlength(2);

  f.setlength(3);
  x0.setlength(3);

  // ------------------------------------------------------------------------------------------
  // Initialize control loop matrices
  // ------------------------------------------------------------------------------------------
  A.setlength(3, 3);
  Q.setlength(3, 3);
  R.setlength(3, 3);
  Gamma.setlength(3, 3);
  GammaT.setlength(3, 3);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        A[i][j] = V;
        Q[i][j] = Q_param;
        R[i][j] = R_param;
      } else {
        A[i][j] = 0;
        Q[i][j] = 0;
        R[i][j] = 0;
      }
    }
  }
}

void CtrModel::ControlLoop() {
  simulation_done = true;
  SetInitialValues();

  double t = 0;
  int rounds = 0;
  auto elapsed_time = 0;
  while (t < total_time) {

    high_resolution_clock::time_point start = high_resolution_clock::now();
    // ------------------------------------------------------------------------------------------
    // Update the input variables
    // ------------------------------------------------------------------------------------------
    q[0] = V * t;
    q[1] = 2 * M_PI * t / 10;

    f[0] = force[0] * sin(2 * M_PI * t / 10);
    f[1] = force[0] * cos(2 * M_PI * t / 10);
    f[2] = force[0] * sin(2 * M_PI * t / 10);

    x0[0] = u_0[0];
    x0[1] = u_0[1];
    x0[2] = u_0[2];
    int b = 3;

    calculation.CalculateShape(x0, q, f);

    // ------------------------------------------------------------------------------------------
    // Update Gamma from ODE results
    // ------------------------------------------------------------------------------------------
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Gamma[i][j] = calculation.GammaFinal[i * 3 + j];
        GammaT[j][i] = Gamma[i][j];
      }
    }

    real_2d_array u_prev;
    u_prev.setlength(3, 1);
    u_prev[0][0] = calculation.uFinal[0];
    u_prev[1][0] = calculation.uFinal[1];
    u_prev[2][0] = calculation.uFinal[2];

    // ------------------------------------------------------------------------------------------
    // dP = -2PA - P*Gamma*R*Gamma^T*P + Q
    // ------------------------------------------------------------------------------------------
    real_2d_array dP;
    dP.setlength(3, 3);
    real_2d_array T1;
    T1.setlength(3, 3);
    real_2d_array T2_1;
    T2_1.setlength(3, 3);
    real_2d_array T2_2;
    T2_2.setlength(3, 3);

    rmatrixgemm(3, 3, 3, 1, P, 0, 0, 0, A, 0, 0, 0, 0, T1, 0, 0);
    rmatrixgemm(3, 3, 3, 1, GammaT, 0, 0, 0, P, 0, 0, 0, 0, T2_1, 0, 0);
    rmatrixgemm(3, 3, 3, 1, R, 0, 0, 0, T2_1, 0, 0, 0, 0, T2_2, 0, 0);
    rmatrixgemm(3, 3, 3, 1, Gamma, 0, 0, 0, T2_2, 0, 0, 0, 0, T2_1, 0, 0);
    rmatrixgemm(3, 3, 3, 1, P, 0, 0, 0, T2_1, 0, 0, 0, 0, T2_2, 0, 0);

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        dP[i][j] = -2 * T1[i][j] - T2_2[i][j] + Q[i][j];
      }
    }

    // ------------------------------------------------------------------------------------------
    // du0 = -R*Gamma^T*P*u
    // ------------------------------------------------------------------------------------------
    real_2d_array du_0;
    du_0.setlength(3, 1);

    real_2d_array Temp_1;
    Temp_1.setlength(3, 1);
    real_2d_array Temp_2;
    Temp_2.setlength(3, 1);

    rmatrixgemm(3, 1, 3, 1, P, 0, 0, 0, u_prev, 0, 0, 0, 0, Temp_1, 0, 0);
    rmatrixgemm(3, 1, 3, 1, GammaT, 0, 0, 0, Temp_1, 0, 0, 0, 0, Temp_2, 0, 0);
    rmatrixgemm(3, 1, 3, -1, R, 0, 0, 0, Temp_2, 0, 0, 0, 0, du_0, 0, 0);

    // ------------------------------------------------------------------------------------------
    // P = P + dt*dP
    // ------------------------------------------------------------------------------------------

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        P[i][j] += dt * dP[i][j];
      }
    }

    // ------------------------------------------------------------------------------------------
    // u0 = u0 + dt*du0
    // ------------------------------------------------------------------------------------------

    for (int i = 0; i < 3; i++) {
      u_0[i] += dt * du_0[i][0];
    }

    // ------------------------------------------------------------------------------------------
    // Store tip position values
    // ------------------------------------------------------------------------------------------
    tip_pos[0].push_back(
        calculation.result_pos[0][calculation.result_pos->size() - 1]);
    tip_pos[1].push_back(
        calculation.result_pos[1][calculation.result_pos->size() - 1]);
    tip_pos[2].push_back(
        calculation.result_pos[2][calculation.result_pos->size() - 1]);

    high_resolution_clock::time_point end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start).count();
    elapsed_time += duration;
    rounds++;
    ROS_DEBUG("Elapsed time of calculation: %ld [us]\n", duration);
    // ------------------------------------------------------------------------------------------
    // Update the robot's shape in RVIZ
    // ------------------------------------------------------------------------------------------
    UpdateRviz();
    t += dt;
  }
  double average = elapsed_time / rounds;
  ROS_INFO("Average time of a round: %f [us]\n", average);
  SaveTipPositionToFile();
}

void CtrModel::callback(ctr_model::MyParamsConfig &config, uint32_t level) {
  // ------------------------------------------------------------------------------------------
  // Update the input parameters according to dynamic reconfigure GUI
  // ------------------------------------------------------------------------------------------

  start_simulation = config.start_simulation;

  force[0] = config.f_x;
  force[0] = config.f_y;
  force[0] = config.f_z;

  total_time = config.total_time;
  dt = config.dt;

  V = config.V;
  Q_param = config.Q;
  R_param = config.R;
}

void CtrModel::UpdateRviz() {
  // ------------------------------------------------------------------------------------------
  // Go through the robot's links and update them with the new coordinates
  // ------------------------------------------------------------------------------------------
  int i = 1;
  while (i < 248) {
    string link_name = "snake_body_";
    link_name += std::to_string(i);

    static tf2_ros::StaticTransformBroadcaster static_broadcaster;
    geometry_msgs::TransformStamped static_transformStamped;

    static_transformStamped.header.stamp = ros::Time::now();
    static_transformStamped.header.frame_id = "base_link";
    static_transformStamped.child_frame_id = link_name.c_str();

    if (i < calculation.result_pos->size()) {
      static_transformStamped.transform.translation.x =
          calculation.result_pos[0][i] * 10;
      static_transformStamped.transform.translation.y =
          calculation.result_pos[1][i] * 10;
      static_transformStamped.transform.translation.z =
          calculation.result_pos[2][i] * 10;
    } else {
      static_transformStamped.transform.translation.x = 0;
      static_transformStamped.transform.translation.y = 0;
      static_transformStamped.transform.translation.z = 0;
    }
    tf2::Quaternion quat;
    quat.setRPY(0, 0, 0);
    static_transformStamped.transform.rotation.x = quat.x();
    static_transformStamped.transform.rotation.y = quat.y();
    static_transformStamped.transform.rotation.z = quat.z();
    static_transformStamped.transform.rotation.w = quat.w();
    static_broadcaster.sendTransform(static_transformStamped);

    i++;
  }
}
void CtrModel::SaveTipPositionToFile() {
  ofstream tip_file;
  tip_file.open("tip.csv");
  for (int i = 0; i < tip_pos->size() - 1; i++) {
    tip_file << tip_pos[0][i] << "," << tip_pos[1][i] << "," << tip_pos[2][i]
             << "\n";
  }
  tip_file.close();
}

CtrModel::~CtrModel() {}

} // namespace ctr_model
