/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_997765728464172060);
void inv_err_fun(double *nom_x, double *true_x, double *out_4666118761844347469);
void H_mod_fun(double *state, double *out_248191370058292062);
void f_fun(double *state, double dt, double *out_5158488645983717192);
void F_fun(double *state, double dt, double *out_7534776359194035090);
void h_25(double *state, double *unused, double *out_5493875580489914981);
void H_25(double *state, double *unused, double *out_2269981339702852985);
void h_24(double *state, double *unused, double *out_3271689822600034662);
void H_24(double *state, double *unused, double *out_4772064521870239053);
void h_30(double *state, double *unused, double *out_8045353141191146893);
void H_30(double *state, double *unused, double *out_7343917571246330295);
void h_26(double *state, double *unused, double *out_610209706762446080);
void H_26(double *state, double *unused, double *out_7974574494642062258);
void h_27(double *state, double *unused, double *out_492202312829069336);
void H_27(double *state, double *unused, double *out_8631499559082955607);
void h_29(double *state, double *unused, double *out_1942833705457952503);
void H_29(double *state, double *unused, double *out_7279027733334389368);
void h_28(double *state, double *unused, double *out_6043363231293276562);
void H_28(double *state, double *unused, double *out_2245767956366255813);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
