#pragma once
#include "alglib/optimization.h"
#include "alglib/stdafx.h"
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

using namespace alglib;
using namespace std;
using namespace boost::numeric::odeint;
using namespace std::chrono;

namespace ctr_model {
class Calculation {
  typedef std::vector<double> state_type;

public:
  state_type y;
  vector<double> result_pos[3];
  double GammaFinal[9];
  double uFinal[3];
  vector<double> y_all[38];

  Calculation();
  virtual ~Calculation();
  void CalculateShape(real_1d_array &ui, real_1d_array &qi, real_1d_array &fi);

private:
  real_2d_array e3_hat;
  real_2d_array K;

  runge_kutta4<state_type> stepper_t;
  void SetDefaultValues(int s_n);
  void calculate_constants(int s_n);
  void UpdateResults(const state_type y, const double t);
  void Ode(state_type &dydt, state_type &y, double t);

  double l;
  double E;
  double r1;
  double r2;
  double J;
  double I;
  double G;

  double U[3];
  double q0[2];
  double force[3];
  double EI;
  double GJ;

  vector<double> dU0[3];
  vector<double> r[3];
  vector<double> du0[9];
  vector<double> u[3];
};

} // namespace ctr_model