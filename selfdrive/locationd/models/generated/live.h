/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3695872146000351190);
void inv_err_fun(double *nom_x, double *true_x, double *out_62932209403646084);
void H_mod_fun(double *state, double *out_8541261767442969313);
void f_fun(double *state, double dt, double *out_4817811066743146942);
void F_fun(double *state, double dt, double *out_3530227812499564946);
void h_3(double *state, double *unused, double *out_4404741154282027977);
void H_3(double *state, double *unused, double *out_2274951906843450349);
void h_4(double *state, double *unused, double *out_8944305960920422952);
void H_4(double *state, double *unused, double *out_1860339505109262378);
void h_9(double *state, double *unused, double *out_7923690389575555455);
void H_9(double *state, double *unused, double *out_7508195327396450417);
void h_10(double *state, double *unused, double *out_7308023086637129461);
void H_10(double *state, double *unused, double *out_914261312880024351);
void h_12(double *state, double *unused, double *out_2574059565779885573);
void H_12(double *state, double *unused, double *out_246216006629947191);
void h_31(double *state, double *unused, double *out_7433660721431237078);
void H_31(double *state, double *unused, double *out_1969483614470895498);
void h_32(double *state, double *unused, double *out_5696575710283579842);
void H_32(double *state, double *unused, double *out_745204282745597285);
void h_13(double *state, double *unused, double *out_6905421030680352239);
void H_13(double *state, double *unused, double *out_977815335710003070);
void h_14(double *state, double *unused, double *out_7923690389575555455);
void H_14(double *state, double *unused, double *out_7508195327396450417);
void h_19(double *state, double *unused, double *out_3083198749340139022);
void H_19(double *state, double *unused, double *out_2107511512226301890);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);