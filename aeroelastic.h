#ifndef MYH_
#define MYH_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

// Useful Macros
#define cols 2
#define AEROFILE "aero.txt"
#define TRAJFILE "traj.txt"
#define SETFILE "set.txt"

// struct declaration

// case_paramaeters refer to the time integration and initial conditions
typedef struct
{

    int ntime;
    double dt;
    double u_inf;
    double alpha;
    double w_eff_0;
    double a_eff_0;
    double w_eff_pu;
    double a_eff_pu;
    double w_eff_nu;
    double a_eff_nu;
    double w_eff_pw;
    double a_eff_pw;
    double w_eff_nw;
    double a_eff_nw;
    double inc;

    gsl_vector *u_0;
    gsl_vector *u_dot_0;
    gsl_vector *u_dot_dot_0;

    bool is_steady;

} case_par;

typedef struct

{
    gsl_vector *aeff;
    gsl_vector *CL_t;
    gsl_vector *CD_t;
    gsl_vector *y1;
    gsl_vector *y2;

} theodorsen;

// state_space refers to the components of the linearized state_space system
// u, derivatives of u , reference angle of attack and velocity will be changed and saved over time
// M,C,K will remain constant following the assumption for linearization
typedef struct
{
    bool is_alloc;
    bool is_init;
    gsl_matrix *u;
    gsl_matrix *u_dot;
    gsl_matrix *u_dot_dot;
    double W_ref;
    double A_ref;
    double CL_a;
    double CD;
    double CL;
    double CD_a;
    double CL_a_pu;
    double CD_pu;
    double CL_pu;
    double CD_a_pu;
    double CL_a_nu;
    double CD_nu;
    double CL_nu;
    double CD_a_nu;
    double CL_a_pw;
    double CD_pw;
    double CL_pw;
    double CD_a_pw;
    double CL_a_nw;
    double CD_nw;
    double CL_nw;
    double CD_a_nw;
    gsl_matrix *M;
    gsl_matrix *C;
    gsl_matrix *K;
    gsl_vector *Q;

} state_space;

// newmark refers to the temporary parameters of each newmark
// timestep
typedef struct
{
    gsl_matrix *A;
    gsl_vector *B;
    gsl_vector *pred;
    gsl_vector *pred_dot;
    gsl_vector *pred_dot_dot;
} newmark;

// foil refers to the structural
// properties of the foil that are needed for
// the filling of the state space m-c-k matrices
typedef struct
{
    double mass;
    double theta_s;
    double theta_t;
    double k_ksi;
    double k_zeta;
} foil;

// Function Declaration w
double deg2rad(double alpha);
double rad2deg(double alpha);
void find_aero(state_space *pstate, case_par *pcase);
// order of execution
void initialize(double time, double dt, double a_o_a, case_par *pcase_p, state_space *pstate, newmark *pnew, theodorsen *ptheo, bool is_steady);
void foil_properties(foil *pfoil, float mass, float theta_s, float theta_t, float k_ksi, float k_zeta);
void initial_conditions(case_par *pcase_p, foil *pfoil, theodorsen *ptheo, double v_0, double v_dot_0, double v_dot_dot_0, double w_0, double w_dot_0, double w_dot_dot_0, double u_inf, double incr);
void make_foil(foil *pfoil, state_space *pstate, case_par *pcase_p);
void newmark_integration(state_space *pstate, newmark *pnew, case_par *pcase_p);

void runge_kutta(state_space *pstate, case_par *pcase_p, theodorsen *ptheo, foil *pfoil);

void terminate(case_par *pcase, state_space *pstate, newmark *pnew);

#endif