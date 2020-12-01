
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_997765728464172060) {
   out_997765728464172060[0] = delta_x[0] + nom_x[0];
   out_997765728464172060[1] = delta_x[1] + nom_x[1];
   out_997765728464172060[2] = delta_x[2] + nom_x[2];
   out_997765728464172060[3] = delta_x[3] + nom_x[3];
   out_997765728464172060[4] = delta_x[4] + nom_x[4];
   out_997765728464172060[5] = delta_x[5] + nom_x[5];
   out_997765728464172060[6] = delta_x[6] + nom_x[6];
   out_997765728464172060[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4666118761844347469) {
   out_4666118761844347469[0] = -nom_x[0] + true_x[0];
   out_4666118761844347469[1] = -nom_x[1] + true_x[1];
   out_4666118761844347469[2] = -nom_x[2] + true_x[2];
   out_4666118761844347469[3] = -nom_x[3] + true_x[3];
   out_4666118761844347469[4] = -nom_x[4] + true_x[4];
   out_4666118761844347469[5] = -nom_x[5] + true_x[5];
   out_4666118761844347469[6] = -nom_x[6] + true_x[6];
   out_4666118761844347469[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_248191370058292062) {
   out_248191370058292062[0] = 1.0;
   out_248191370058292062[1] = 0.0;
   out_248191370058292062[2] = 0.0;
   out_248191370058292062[3] = 0.0;
   out_248191370058292062[4] = 0.0;
   out_248191370058292062[5] = 0.0;
   out_248191370058292062[6] = 0.0;
   out_248191370058292062[7] = 0.0;
   out_248191370058292062[8] = 0.0;
   out_248191370058292062[9] = 1.0;
   out_248191370058292062[10] = 0.0;
   out_248191370058292062[11] = 0.0;
   out_248191370058292062[12] = 0.0;
   out_248191370058292062[13] = 0.0;
   out_248191370058292062[14] = 0.0;
   out_248191370058292062[15] = 0.0;
   out_248191370058292062[16] = 0.0;
   out_248191370058292062[17] = 0.0;
   out_248191370058292062[18] = 1.0;
   out_248191370058292062[19] = 0.0;
   out_248191370058292062[20] = 0.0;
   out_248191370058292062[21] = 0.0;
   out_248191370058292062[22] = 0.0;
   out_248191370058292062[23] = 0.0;
   out_248191370058292062[24] = 0.0;
   out_248191370058292062[25] = 0.0;
   out_248191370058292062[26] = 0.0;
   out_248191370058292062[27] = 1.0;
   out_248191370058292062[28] = 0.0;
   out_248191370058292062[29] = 0.0;
   out_248191370058292062[30] = 0.0;
   out_248191370058292062[31] = 0.0;
   out_248191370058292062[32] = 0.0;
   out_248191370058292062[33] = 0.0;
   out_248191370058292062[34] = 0.0;
   out_248191370058292062[35] = 0.0;
   out_248191370058292062[36] = 1.0;
   out_248191370058292062[37] = 0.0;
   out_248191370058292062[38] = 0.0;
   out_248191370058292062[39] = 0.0;
   out_248191370058292062[40] = 0.0;
   out_248191370058292062[41] = 0.0;
   out_248191370058292062[42] = 0.0;
   out_248191370058292062[43] = 0.0;
   out_248191370058292062[44] = 0.0;
   out_248191370058292062[45] = 1.0;
   out_248191370058292062[46] = 0.0;
   out_248191370058292062[47] = 0.0;
   out_248191370058292062[48] = 0.0;
   out_248191370058292062[49] = 0.0;
   out_248191370058292062[50] = 0.0;
   out_248191370058292062[51] = 0.0;
   out_248191370058292062[52] = 0.0;
   out_248191370058292062[53] = 0.0;
   out_248191370058292062[54] = 1.0;
   out_248191370058292062[55] = 0.0;
   out_248191370058292062[56] = 0.0;
   out_248191370058292062[57] = 0.0;
   out_248191370058292062[58] = 0.0;
   out_248191370058292062[59] = 0.0;
   out_248191370058292062[60] = 0.0;
   out_248191370058292062[61] = 0.0;
   out_248191370058292062[62] = 0.0;
   out_248191370058292062[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_5158488645983717192) {
   out_5158488645983717192[0] = state[0];
   out_5158488645983717192[1] = state[1];
   out_5158488645983717192[2] = state[2];
   out_5158488645983717192[3] = state[3];
   out_5158488645983717192[4] = state[4];
   out_5158488645983717192[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5158488645983717192[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5158488645983717192[7] = state[7];
}
void F_fun(double *state, double dt, double *out_7534776359194035090) {
   out_7534776359194035090[0] = 1;
   out_7534776359194035090[1] = 0;
   out_7534776359194035090[2] = 0;
   out_7534776359194035090[3] = 0;
   out_7534776359194035090[4] = 0;
   out_7534776359194035090[5] = 0;
   out_7534776359194035090[6] = 0;
   out_7534776359194035090[7] = 0;
   out_7534776359194035090[8] = 0;
   out_7534776359194035090[9] = 1;
   out_7534776359194035090[10] = 0;
   out_7534776359194035090[11] = 0;
   out_7534776359194035090[12] = 0;
   out_7534776359194035090[13] = 0;
   out_7534776359194035090[14] = 0;
   out_7534776359194035090[15] = 0;
   out_7534776359194035090[16] = 0;
   out_7534776359194035090[17] = 0;
   out_7534776359194035090[18] = 1;
   out_7534776359194035090[19] = 0;
   out_7534776359194035090[20] = 0;
   out_7534776359194035090[21] = 0;
   out_7534776359194035090[22] = 0;
   out_7534776359194035090[23] = 0;
   out_7534776359194035090[24] = 0;
   out_7534776359194035090[25] = 0;
   out_7534776359194035090[26] = 0;
   out_7534776359194035090[27] = 1;
   out_7534776359194035090[28] = 0;
   out_7534776359194035090[29] = 0;
   out_7534776359194035090[30] = 0;
   out_7534776359194035090[31] = 0;
   out_7534776359194035090[32] = 0;
   out_7534776359194035090[33] = 0;
   out_7534776359194035090[34] = 0;
   out_7534776359194035090[35] = 0;
   out_7534776359194035090[36] = 1;
   out_7534776359194035090[37] = 0;
   out_7534776359194035090[38] = 0;
   out_7534776359194035090[39] = 0;
   out_7534776359194035090[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7534776359194035090[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7534776359194035090[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7534776359194035090[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7534776359194035090[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7534776359194035090[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7534776359194035090[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7534776359194035090[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7534776359194035090[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7534776359194035090[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7534776359194035090[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7534776359194035090[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7534776359194035090[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7534776359194035090[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7534776359194035090[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7534776359194035090[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7534776359194035090[56] = 0;
   out_7534776359194035090[57] = 0;
   out_7534776359194035090[58] = 0;
   out_7534776359194035090[59] = 0;
   out_7534776359194035090[60] = 0;
   out_7534776359194035090[61] = 0;
   out_7534776359194035090[62] = 0;
   out_7534776359194035090[63] = 1;
}
void h_25(double *state, double *unused, double *out_5493875580489914981) {
   out_5493875580489914981[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2269981339702852985) {
   out_2269981339702852985[0] = 0;
   out_2269981339702852985[1] = 0;
   out_2269981339702852985[2] = 0;
   out_2269981339702852985[3] = 0;
   out_2269981339702852985[4] = 0;
   out_2269981339702852985[5] = 0;
   out_2269981339702852985[6] = 1;
   out_2269981339702852985[7] = 0;
}
void h_24(double *state, double *unused, double *out_3271689822600034662) {
   out_3271689822600034662[0] = state[4];
   out_3271689822600034662[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4772064521870239053) {
   out_4772064521870239053[0] = 0;
   out_4772064521870239053[1] = 0;
   out_4772064521870239053[2] = 0;
   out_4772064521870239053[3] = 0;
   out_4772064521870239053[4] = 1;
   out_4772064521870239053[5] = 0;
   out_4772064521870239053[6] = 0;
   out_4772064521870239053[7] = 0;
   out_4772064521870239053[8] = 0;
   out_4772064521870239053[9] = 0;
   out_4772064521870239053[10] = 0;
   out_4772064521870239053[11] = 0;
   out_4772064521870239053[12] = 0;
   out_4772064521870239053[13] = 1;
   out_4772064521870239053[14] = 0;
   out_4772064521870239053[15] = 0;
}
void h_30(double *state, double *unused, double *out_8045353141191146893) {
   out_8045353141191146893[0] = state[4];
}
void H_30(double *state, double *unused, double *out_7343917571246330295) {
   out_7343917571246330295[0] = 0;
   out_7343917571246330295[1] = 0;
   out_7343917571246330295[2] = 0;
   out_7343917571246330295[3] = 0;
   out_7343917571246330295[4] = 1;
   out_7343917571246330295[5] = 0;
   out_7343917571246330295[6] = 0;
   out_7343917571246330295[7] = 0;
}
void h_26(double *state, double *unused, double *out_610209706762446080) {
   out_610209706762446080[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7974574494642062258) {
   out_7974574494642062258[0] = 0;
   out_7974574494642062258[1] = 0;
   out_7974574494642062258[2] = 0;
   out_7974574494642062258[3] = 0;
   out_7974574494642062258[4] = 0;
   out_7974574494642062258[5] = 0;
   out_7974574494642062258[6] = 0;
   out_7974574494642062258[7] = 1;
}
void h_27(double *state, double *unused, double *out_492202312829069336) {
   out_492202312829069336[0] = state[3];
}
void H_27(double *state, double *unused, double *out_8631499559082955607) {
   out_8631499559082955607[0] = 0;
   out_8631499559082955607[1] = 0;
   out_8631499559082955607[2] = 0;
   out_8631499559082955607[3] = 1;
   out_8631499559082955607[4] = 0;
   out_8631499559082955607[5] = 0;
   out_8631499559082955607[6] = 0;
   out_8631499559082955607[7] = 0;
}
void h_29(double *state, double *unused, double *out_1942833705457952503) {
   out_1942833705457952503[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7279027733334389368) {
   out_7279027733334389368[0] = 0;
   out_7279027733334389368[1] = 1;
   out_7279027733334389368[2] = 0;
   out_7279027733334389368[3] = 0;
   out_7279027733334389368[4] = 0;
   out_7279027733334389368[5] = 0;
   out_7279027733334389368[6] = 0;
   out_7279027733334389368[7] = 0;
}
void h_28(double *state, double *unused, double *out_6043363231293276562) {
   out_6043363231293276562[0] = state[5];
   out_6043363231293276562[1] = state[6];
}
void H_28(double *state, double *unused, double *out_2245767956366255813) {
   out_2245767956366255813[0] = 0;
   out_2245767956366255813[1] = 0;
   out_2245767956366255813[2] = 0;
   out_2245767956366255813[3] = 0;
   out_2245767956366255813[4] = 0;
   out_2245767956366255813[5] = 1;
   out_2245767956366255813[6] = 0;
   out_2245767956366255813[7] = 0;
   out_2245767956366255813[8] = 0;
   out_2245767956366255813[9] = 0;
   out_2245767956366255813[10] = 0;
   out_2245767956366255813[11] = 0;
   out_2245767956366255813[12] = 0;
   out_2245767956366255813[13] = 0;
   out_2245767956366255813[14] = 1;
   out_2245767956366255813[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
