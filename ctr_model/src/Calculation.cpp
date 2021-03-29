#include "ctr_model/Calculation.hpp"

namespace ctr_model {
Calculation::Calculation() {
  // ------------------------------------------------------------------------------------------
  // Set model parameters and default inputs
  // ------------------------------------------------------------------------------------------
  l = 0.4;
  E = 7e+10;
  r1 = 0.0015;
  r2 = 1e-3;
  J = M_PI / 2 * (pow(r1, 4) - pow(r2, 4));
  I = 0.5 * J;
  G = 10e9;

  U[0] = 14;
  U[1] = 5;
  U[2] = 0;
  q0[0] = -0.2;
  q0[1] = 0;
  force[0] = 0;
  force[1] = 0;
  force[2] = 0;
  EI = E * I;
  GJ = G * J;
}

Calculation::~Calculation() {}

void Calculation::SetDefaultValues(int s_n) {
  // ------------------------------------------------------------------------------------------
  // Initialization
  // ------------------------------------------------------------------------------------------
  y.resize(0);
  y.shrink_to_fit();
  for (int i = 0; i < 3; i++) {
    r[i].resize(0);
    r[i].shrink_to_fit();
    du0[i].resize(0);
    du0[i].shrink_to_fit();
    u[i].resize(0);
    u[i].shrink_to_fit();
    result_pos[i].resize(0);
    result_pos[i].shrink_to_fit();
    uFinal[i] = 0;
  }
  for (int i = 0; i < 9; i++) {
    GammaFinal[i] = 0;
  }
  for (int i = 0; i < 33; i++) {
    y_all[i].resize(0);
    y_all[i].shrink_to_fit();
  }
  // ------------------------------------------------------------------------------------------
  // e3_hat
  // ------------------------------------------------------------------------------------------
  e3_hat.setlength(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      e3_hat[i][j] = 0;
    }
  }
  e3_hat[0][1] = -1;
  e3_hat[1][0] = 1;
  // ------------------------------------------------------------------------------------------
  // K
  // ------------------------------------------------------------------------------------------
  K.setlength(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      K[i][j] = 0;
    }
  }
  K[0][0] = EI;
  K[1][1] = EI;
  K[2][2] = GJ;

  return;
}

void Calculation::UpdateResults(const state_type y, const double t) {

  for (int i = 0; i < 33; i++) {
    y_all[i].push_back(y[i]);
  }
  result_pos[0].push_back(y[3]);
  result_pos[1].push_back(y[4]);
  result_pos[2].push_back(y[5]);
  uFinal[0] = y[0];
  uFinal[1] = y[1];
  uFinal[2] = y[2];

  for (int i = 0; i < 9; i++) {
    GammaFinal[i] = y[15 + i];
  }

  return;
}

void Calculation::Ode(state_type &dydt, state_type &y, double t) {

  // ------------------------------------------------------------------------------------------
  // Set input values
  // ------------------------------------------------------------------------------------------
  real_1d_array u;
  u.setlength(3);
  u[0] = y[0];
  u[1] = y[1];
  u[2] = y[2];

  double temp_u_hat[3][3] = {
      {0, -u[2], u[1]}, {u[2], 0, -u[0]}, {-u[1], u[0], 0}};
  real_2d_array u_hat;
  u_hat.setlength(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      u_hat[i][j] = temp_u_hat[i][j];
    }
  }

  real_2d_array R;
  R.setlength(3, 3);
  real_2d_array Gamma;
  Gamma.setlength(3, 3);
  real_2d_array Chi;
  Chi.setlength(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      R[i][j] = y[6 + i * 3 + j];
      Gamma[i][j] = y[15 + i * 3 + j];
      Chi[i][j] = y[24 + i * 3 + j];
    }
  }

  // ------------------------------------------------------------------------------------------
  // C = R^T*f  -> c_hat
  // ------------------------------------------------------------------------------------------
  real_2d_array C;
  C.setlength(3, 1);

  for (int i = 0; i < 3; i++) {
    C[i][0] = 0;
    for (int j = 0; j < 3; j++) {
      C[i][0] += R[j][i] * force[j];
    }
  }

  double temp[3][3] = {
      {0, -C[2][0], C[1][0]}, {C[2][0], 0, -C[0][0]}, {-C[1][0], C[0][0], 0}};
  real_2d_array c_hat;
  c_hat.setlength(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      c_hat[i][j] = temp[i][j];
    }
  }
  // ------------------------------------------------------------------------------------------
  // dChi = c_hat*Gamma - u_hat*Chi
  // ------------------------------------------------------------------------------------------
  real_2d_array dChi;
  dChi.setlength(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dChi[i][j] = 0;
      for (int k = 0; k < 3; k++)
        dChi[i][j] += c_hat[i][k] * Gamma[k][j] - u_hat[i][k] * Chi[k][j];
    }
  }

  // ------------------------------------------------------------------------------------------
  // A = K*(U-Us)  -> a_hat
  // ------------------------------------------------------------------------------------------
  real_2d_array A;
  A.setlength(3, 1);
  A[0][0] = (u[0] - U[0]) * EI;
  A[1][0] = (u[1] - U[1]) * EI;
  A[2][0] = (u[2] - U[2]) * GJ;
  double temp2[3][3] = {
      {0, -A[2][0], A[1][0]}, {A[2][0], 0, -A[0][0]}, {-A[1][0], A[0][0], 0}};
  real_2d_array a_hat;
  a_hat.setlength(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      a_hat[i][j] = temp2[i][j];
    }
  }

  // ------------------------------------------------------------------------------------------
  // dGamma = -K^(-1)(-a_hat*Gamma +u_hat*K*Gamma + e3_hat*Chi)
  // ------------------------------------------------------------------------------------------
  real_2d_array dGamma;
  dGamma.setlength(3, 3);
  // Temporary arrays
  real_2d_array T1;
  T1.setlength(3, 3);
  real_2d_array T2_2;
  T2_2.setlength(3, 3);
  real_2d_array T2;
  T2.setlength(3, 3);
  real_2d_array T3;
  T3.setlength(3, 3);

  rmatrixgemm(3, 3, 3, 1, a_hat, 0, 0, 0, Gamma, 0, 0, 0, 0, T1, 0, 0);
  rmatrixgemm(3, 3, 3, 1, K, 0, 0, 0, Gamma, 0, 0, 0, 0, T2_2, 0, 0);
  rmatrixgemm(3, 3, 3, 1, u_hat, 0, 0, 0, T2_2, 0, 0, 0, 0, T2, 0, 0);
  rmatrixgemm(3, 3, 3, 1, e3_hat, 0, 0, 0, Chi, 0, 0, 0, 0, T3, 0, 0);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dGamma[i][j] = -1 / K[i][i] * (-T1[i][j] + T2[i][j] + T3[i][j]);
    }
  }

  // ------------------------------------------------------------------------------------------
  // du = -K^(-1)*(u_hat*A+e3_hat*C)
  // ------------------------------------------------------------------------------------------
  real_2d_array du;
  du.setlength(3, 1);
  rmatrixgemm(3, 1, 3, 1, u_hat, 0, 0, 0, A, 0, 0, 0, 0, T1, 0, 0);
  rmatrixgemm(3, 1, 3, 1, e3_hat, 0, 0, 0, C, 0, 0, 0, 0, T2, 0, 0);

  for (int i = 0; i < 3; i++) {
    du[i][0] = -1 / K[i][i] * (T1[i][0] + T2[i][0]);
  }

  // ------------------------------------------------------------------------------------------
  // dr = r*e3
  // ------------------------------------------------------------------------------------------
  double dr[3];
  for (int i = 0; i < 3; i++)
    dr[i] = R[i][2];

  // ------------------------------------------------------------------------------------------
  // dR = R*u_hat
  // ------------------------------------------------------------------------------------------
  real_2d_array dR;
  dR.setlength(3, 3);

  rmatrixgemm(3, 3, 3, 1, R, 0, 0, 0, u_hat, 0, 0, 0, 0, dR, 0, 0);

  // ------------------------------------------------------------------------------------------
  // Pass the calculated values to dydt
  // ------------------------------------------------------------------------------------------
  for (int i = 0; i < 3; i++) {
    dydt[i] = du[i][0];
    dydt[3 + i] = dr[i];
    for (int j = 0; j < 3; j++) {
      dydt[6 + i * 3 + j] = dR[i][j];
      dydt[15 + i * 3 + j] = dGamma[i][j];
      dydt[24 + i * 3 + j] = dChi[i][j];
    }
  }
  return;
}

void Calculation::CalculateShape(real_1d_array &ui, real_1d_array &qi,
                                 real_1d_array &fi) {
  // ------------------------------------------------------------------------------------------
  // Update input values
  // ------------------------------------------------------------------------------------------
  SetDefaultValues(33);

  double B = qi[0] + q0[0];
  double alpha = qi[1] + q0[1];

  double r0[3] = {0, 0, 0};
  double R0[9] = {cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha),
                  0,          0,           0, 1};
  double du0[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double s_span = l + B;

  state_type y_0(33);

  for (int i = 0; i < 3; i++) {
    force[i] = fi[i];
    y_0[i] = ui[i];
    y_0[3 + i] = r0[i];
  }
  for (int i = 0; i < 9; i++) {
    y_0[6 + i] = R0[i];
    y_0[15 + i] = du0[i];
    y_0[24 + i] = 0;
  }

  // ------------------------------------------------------------------------------------------
  // Solving the ODE
  // ------------------------------------------------------------------------------------------
  integrate_adaptive(
      stepper_t,
      [this](state_type y, state_type &dydt, double t) {
        this->Ode(dydt, y, t);
      },
      y_0, 0.0, s_span, 0.0050,
      [this](state_type &y, double t) { this->UpdateResults(y, t); });

  return;
}

} // namespace ctr_model