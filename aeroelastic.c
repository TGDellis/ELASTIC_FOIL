#include "aeroelastic.h"

double rad2deg(double alpha)
{
    return 180 * alpha / M_PI;
}
double deg2rad(double alpha)
{
    return alpha * M_PI / 180;
}
void initialize(double time, double dt, double a_o_a, case_par *pcase_p, state_space *pstate, newmark *pnew, theodorsen *ptheo, bool is_steady)
{
    int i;

    pcase_p->dt = dt;
    pcase_p->ntime = (int)round(time / dt);
    printf("\n ntime = %d\n", pcase_p->ntime);
    pcase_p->u_0 = gsl_vector_alloc(2);
    pcase_p->u_dot_0 = gsl_vector_alloc(2);
    pcase_p->u_dot_dot_0 = gsl_vector_alloc(2);
    pcase_p->alpha = a_o_a;
    pcase_p->is_steady = is_steady;

    pstate->u = gsl_matrix_calloc(pcase_p->ntime, cols);
    pstate->u_dot = gsl_matrix_calloc(pcase_p->ntime, cols);
    pstate->u_dot_dot = gsl_matrix_calloc(pcase_p->ntime, cols);

    pstate->M = gsl_matrix_calloc(cols, cols);
    if (pcase_p->is_steady == true)
    {
        pstate->C = gsl_matrix_calloc(cols, cols);
        printf("\nsteady calculation\n");
    }
    else
    {
        ptheo->aeff = gsl_vector_calloc(pcase_p->ntime);
        ptheo->CD_t = gsl_vector_calloc(pcase_p->ntime);
        ptheo->CL_t = gsl_vector_calloc(pcase_p->ntime);
        ptheo->y1 = gsl_vector_calloc(pcase_p->ntime);
        ptheo->y2 = gsl_vector_calloc(pcase_p->ntime);
        printf("\nunsteady calculation\n");
    }
    pstate->K = gsl_matrix_calloc(cols, cols);
    pstate->Q = gsl_vector_calloc(cols);

    pnew->A = gsl_matrix_calloc(cols, cols);
    pnew->B = gsl_vector_calloc(cols);
    pnew->pred = gsl_vector_calloc(cols);
    pnew->pred_dot = gsl_vector_calloc(cols);
    pnew->pred_dot_dot = gsl_vector_calloc(cols);
    pstate->is_alloc = true;
}

void terminate(case_par *pcase_p, state_space *pstate, newmark *pnew)
{

    gsl_matrix_free(pstate->u);
    gsl_matrix_free(pstate->u_dot);
    gsl_matrix_free(pstate->u_dot_dot);

    gsl_matrix_free(pstate->M);
    // gsl_matrix_free(pstate->C);
    gsl_matrix_free(pstate->K);
    // gsl_vector_free(pstate->Q);

    // gsl_matrix_free(pnew->A);
    // gsl_vector_free(pnew->B);
    // gsl_vector_free(pnew->pred);
    // gsl_vector_free(pnew->pred_dot);
    // gsl_vector_free(pnew->pred_dot_dot);
}
void foil_properties(foil *pfoil, float mass, float theta_s, float theta_t, float k_ksi, float k_zeta)
{
    pfoil->mass = mass;
    pfoil->k_ksi = k_ksi;
    pfoil->k_zeta = k_zeta;
    pfoil->theta_t = theta_t;
    pfoil->theta_s = theta_s;
}
void make_foil(foil *pfoil, state_space *pstate, case_par *pcase_p)
{
    // Airfoil mass
    gsl_matrix_set(pstate->M, 0, 0, pfoil->mass);
    gsl_matrix_set(pstate->M, 1, 1, pfoil->mass);
    gsl_matrix_set(pstate->M, 0, 1, 0.0);
    gsl_matrix_set(pstate->M, 1, 0, 0.0);

    // K matrix - Stiffness depending on the airfoil's elastic properties
    double k_xx = pfoil->k_ksi * pow(cos(pfoil->theta_s), 2) + pfoil->k_zeta * pow(sin(pfoil->theta_s), 2);
    double k_xz = -pfoil->k_ksi * cos(pfoil->theta_s) * sin(pfoil->theta_s) + pfoil->k_zeta * cos(pfoil->theta_s) * sin(pfoil->theta_s);
    double k_zz = pfoil->k_zeta * pow(cos(pfoil->theta_s), 2) + pfoil->k_ksi * pow(sin(pfoil->theta_s), 2);
    gsl_matrix_set(pstate->K, 0, 0, k_xx);
    gsl_matrix_set(pstate->K, 0, 1, k_xz);
    gsl_matrix_set(pstate->K, 1, 0, k_xz);
    gsl_matrix_set(pstate->K, 1, 1, k_zz);

    if (pcase_p->is_steady == true)
    {
        // C matrix - Aerodynamic damping depend on the angle of attack where the system is linearized
        find_aero(pstate, pcase_p);
        printf("\nCL_a= %g\n", pstate->CL_a);
        printf("\nCD_a= %g\n", pstate->CD_a);
        printf("\nu_inf= %g\n", pcase_p->u_inf);
        printf("\nw_eff_0= %g\n", pcase_p->w_eff_0);
        printf("\naoa= %lf\n", rad2deg(pcase_p->alpha));
        printf("\na_eff_0= %lf\n", rad2deg(pcase_p->a_eff_0));
        double dW_dv_dot = 2 * (gsl_vector_get(pcase_p->u_dot_0, 0) + pcase_p->u_inf * cos(pfoil->theta_t + pcase_p->alpha));
        double dW_dw_dot = 2 * (gsl_vector_get(pcase_p->u_dot_0, 1) - pcase_p->u_inf * sin(pfoil->theta_t + pcase_p->alpha));
        double da_dv_dot = (gsl_vector_get(pcase_p->u_dot_0, 1) - pcase_p->u_inf * sin(pfoil->theta_t + pcase_p->alpha)) / (pow((gsl_vector_get(pcase_p->u_dot_0, 0) + pcase_p->u_inf * cos(pfoil->theta_t + pcase_p->alpha)), 2) + pow((gsl_vector_get(pcase_p->u_dot_0, 1) - pcase_p->u_inf * sin(pfoil->theta_t + pcase_p->alpha)), 2));
        double da_dw_dot = -(gsl_vector_get(pcase_p->u_dot_0, 0) + pcase_p->u_inf * cos(pfoil->theta_t + pcase_p->alpha)) / (pow((gsl_vector_get(pcase_p->u_dot_0, 0) + pcase_p->u_inf * cos(pfoil->theta_t + pcase_p->alpha)), 2) + pow((gsl_vector_get(pcase_p->u_dot_0, 1) - pcase_p->u_inf * sin(pfoil->theta_t + pcase_p->alpha)), 2));
        printf("\ndw-dv = %lf", da_dw_dot);
        double Fx_v_dot = pstate->CL_a * da_dv_dot * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * sin(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CL * 0.5 * 1 * 1.225 * dW_dv_dot * sin(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CL * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * cos(pcase_p->a_eff_0 + pfoil->theta_t) * da_dv_dot -
                          pstate->CD_a * da_dv_dot * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * cos(pcase_p->a_eff_0 + pfoil->theta_t) -
                          pstate->CD * 0.5 * 1 * 1.225 * dW_dv_dot * cos(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CD * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * sin(pcase_p->a_eff_0 + pfoil->theta_t) * da_dv_dot;
        double Fx_w_dot = pstate->CL_a * da_dw_dot * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * sin(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CL * 0.5 * 1 * 1.225 * dW_dw_dot * sin(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CL * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * cos(pcase_p->a_eff_0 + pfoil->theta_t) * da_dw_dot -
                          pstate->CD_a * da_dw_dot * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * cos(pcase_p->a_eff_0 + pfoil->theta_t) -
                          pstate->CD * 0.5 * 1 * 1.225 * dW_dw_dot * cos(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CD * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * sin(pcase_p->a_eff_0 + pfoil->theta_t) * da_dw_dot;
        double Fz_v_dot = pstate->CL_a * da_dv_dot * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * cos(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CL * 0.5 * 1 * 1.225 * dW_dv_dot * cos(pcase_p->a_eff_0 + pfoil->theta_t) -
                          pstate->CL * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * sin(pcase_p->a_eff_0 + pfoil->theta_t) * da_dv_dot +
                          pstate->CD_a * da_dv_dot * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * sin(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CD * 0.5 * 1 * 1.225 * dW_dv_dot * sin(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CD * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * cos(pcase_p->a_eff_0 + pfoil->theta_t) * da_dv_dot;
        double Fz_w_dot = pstate->CL_a * da_dw_dot * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * cos(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CL * 0.5 * 1 * 1.225 * dW_dw_dot * cos(pcase_p->a_eff_0 + pfoil->theta_t) -
                          pstate->CL * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * sin(pcase_p->a_eff_0 + pfoil->theta_t) * da_dw_dot +
                          pstate->CD_a * da_dw_dot * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * sin(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CD * 0.5 * 1 * 1.225 * dW_dw_dot * sin(pcase_p->a_eff_0 + pfoil->theta_t) +
                          pstate->CD * 0.5 * 1 * 1.225 * pow(pcase_p->w_eff_0, 2) * cos(pcase_p->a_eff_0 + pfoil->theta_t) * da_dw_dot;
        // double L_pu = pstate->CL_pu * 1 * 1.225 * 0.5 * pow(pcase_p->w_eff_pu, 2);
        // double D_pu = pstate->CD_pu * 1 * 1.225 * 0.5 * pow(pcase_p->w_eff_pu, 2);
        // double Fx_pu = L_pu * sin(pcase_p->a_eff_pu + pfoil->theta_t) - D_pu * cos(pcase_p->a_eff_pu + pfoil->theta_t);
        // double Fz_pu = L_pu * cos(pcase_p->a_eff_pu + pfoil->theta_t) + D_pu * sin(pcase_p->a_eff_pu + pfoil->theta_t);

        // double L_nu = pstate->CL_nu * 1 * 1.225 * 0.5 * pow(pcase_p->w_eff_nu, 2);
        // double D_nu = pstate->CD_nu * 1 * 1.225 * 0.5 * pow(pcase_p->w_eff_nu, 2);
        // double Fx_nu = L_nu * sin(pcase_p->a_eff_nu + pfoil->theta_t) - D_nu * cos(pcase_p->a_eff_nu + pfoil->theta_t);
        // double Fz_nu = L_nu * cos(pcase_p->a_eff_nu + pfoil->theta_t) + D_nu * sin(pcase_p->a_eff_nu + pfoil->theta_t);

        // double L_nw = pstate->CL_nw * 1 * 1.225 * 0.5 * pow(pcase_p->w_eff_nw, 2);
        // double D_nw = pstate->CD_nw * 1 * 1.225 * 0.5 * pow(pcase_p->w_eff_nw, 2);
        // double Fx_nw = L_nw * sin(pcase_p->a_eff_nw + pfoil->theta_t) - D_nw * cos(pcase_p->a_eff_nw + pfoil->theta_t);
        // double Fz_nw = L_nw * cos(pcase_p->a_eff_nw + pfoil->theta_t) + D_nw * sin(pcase_p->a_eff_nw + pfoil->theta_t);

        // double L_pw = pstate->CL_pw * 1 * 1.225 * 0.5 * pow(pcase_p->w_eff_pw, 2);
        // double D_pw = pstate->CD_pw * 1 * 1.225 * 0.5 * pow(pcase_p->w_eff_pw, 2);
        // double Fx_pw = L_pw * sin(pcase_p->a_eff_pw + pfoil->theta_t) - D_pw * cos(pcase_p->a_eff_pw + pfoil->theta_t);
        // double Fz_pw = L_pw * cos(pcase_p->a_eff_pw + pfoil->theta_t) + D_pw * sin(pcase_p->a_eff_pw + pfoil->theta_t);

        // double Fx_v_dot = (Fx_nu-Fx_pu)/(2*pcase_p->inc);
        // double Fz_v_dot = (Fz_nu-Fz_pu)/(2*pcase_p->inc);
        // double Fx_w_dot = (Fx_nw-Fx_pw)/(2*pcase_p->inc);
        // double Fz_w_dot = (Fz_nw-Fz_pw)/(2*pcase_p->inc);
        // printf("\nFxv=%g\n",Fx_v_dot);
        gsl_matrix_set(pstate->C, 0, 0, -Fx_v_dot);
        gsl_matrix_set(pstate->C, 0, 1, -Fx_w_dot);
        gsl_matrix_set(pstate->C, 1, 0, -Fz_v_dot);
        gsl_matrix_set(pstate->C, 1, 1, -Fz_w_dot);
        // (-0.006731229802335476, 0.1343171571519688)
        // (0.40009695312515886, -0.1442743264659525)

        // Q source vector - depends on linearization point
        double L_0 = pstate->CL * 1 * 1.225 * 0.5 * pow(pcase_p->w_eff_0, 2);
        double D_0 = pstate->CD * 1 * 1.225 * 0.5 * pow(pcase_p->w_eff_0, 2);
        double Fx_0 = L_0 * sin(pcase_p->alpha + pfoil->theta_t) - D_0 * cos(pcase_p->alpha + pfoil->theta_t);
        double Fz_0 = L_0 * cos(pcase_p->alpha + pfoil->theta_t) + D_0 * sin(pcase_p->alpha + pfoil->theta_t);
        gsl_vector *m_u_dot_dot = gsl_vector_calloc(2);
        gsl_vector *c_u_dot = gsl_vector_calloc(2);
        gsl_vector *k_u = gsl_vector_calloc(2);
        gsl_blas_dgemv(CblasNoTrans, 1.0, pstate->M, pcase_p->u_dot_dot_0, 0.0, m_u_dot_dot);
        gsl_blas_dgemv(CblasNoTrans, 1.0, pstate->C, pcase_p->u_dot_0, 0.0, c_u_dot);
        gsl_blas_dgemv(CblasNoTrans, 1.0, pstate->K, pcase_p->u_0, 0.0, k_u);
        gsl_vector_sub(m_u_dot_dot, c_u_dot);
        gsl_vector_sub(m_u_dot_dot, k_u);
        gsl_vector_set(pstate->Q, 0, 1 * Fx_0);
        gsl_vector_set(pstate->Q, 1, 1 * Fz_0);
        gsl_vector_sub(pstate->Q, m_u_dot_dot);
        printf("\nQ_0 =%lf \n", gsl_vector_get(pstate->Q, 0));
        gsl_vector_free(m_u_dot_dot);
        gsl_vector_free(c_u_dot);
        gsl_vector_free(k_u);

        pstate->is_init = true;
    }
}

void find_aero(state_space *pstate, case_par *pcase)

{

    double angle, CL_interp, CD_interp;
    int len, idx;
    int i;
    gsl_vector *aoa_diff;
    if (pcase->is_steady == true)
    {
        FILE *f_aero = fopen(AEROFILE, "r");
        fscanf(f_aero, "%d", &len);

        aoa_diff = gsl_vector_calloc(len);
        double CL[len], CD[len], CL_a[len], CD_a[len];
        printf("parsing aerofile\n");
        char buffer[100];
        fgets(buffer, 100, f_aero);
        fgets(buffer, 100, f_aero);
        for (i = 0; i < len; i++)
        {

            fscanf(f_aero, "%lf %lf %lf %lf %lf", &angle, &CL[i], &CD[i], &CL_a[i], &CD_a[i]);
            gsl_vector_set(aoa_diff, i, fabs(rad2deg(pcase->a_eff_0) - angle));
        }
        idx = gsl_vector_min_index(aoa_diff);
        fclose(f_aero);

        printf("\nidx is %d\n", idx);
        pstate->CL = CL[idx];
        pstate->CD = CD[idx];
        pstate->CL_a = CL_a[idx];
        pstate->CD_a = CD_a[idx];

        FILE *f_aero_1 = fopen(AEROFILE, "r");
        printf("parsing aerofile\n");
        fgets(buffer, 100, f_aero);
        fgets(buffer, 100, f_aero);
        for (i = 0; i < len; i++)
        {

            fscanf(f_aero_1, "%lf %lf %lf %lf %lf", &angle, &CL[i], &CD[i], &CL_a[i], &CD_a[i]);
            gsl_vector_set(aoa_diff, i, fabs(rad2deg(pcase->a_eff_nu) - angle));
        }
        idx = gsl_vector_min_index(aoa_diff);
        fclose(f_aero_1);
        printf("\nidx is %d\n", idx);
        pstate->CL_nu = CL[idx];
        pstate->CD_nu = CD[idx];
        pstate->CL_a_nu = CL_a[idx];
        pstate->CD_a_nu = CD_a[idx];

        FILE *f_aero_2 = fopen(AEROFILE, "r");
        printf("parsing aerofile\n");
        fgets(buffer, 100, f_aero);
        fgets(buffer, 100, f_aero);
        for (i = 0; i < len; i++)
        {

            fscanf(f_aero_2, "%lf %lf %lf %lf %lf", &angle, &CL[i], &CD[i], &CL_a[i], &CD_a[i]);
            gsl_vector_set(aoa_diff, i, fabs(rad2deg(pcase->a_eff_pu) - angle));
        }
        idx = gsl_vector_min_index(aoa_diff);
        fclose(f_aero_2);
        printf("\nidx is %d\n", idx);
        pstate->CL_pu = CL[idx];
        pstate->CD_pu = CD[idx];
        pstate->CL_a_pu = CL_a[idx];
        pstate->CD_a_pu = CD_a[idx];

        FILE *f_aero_3 = fopen(AEROFILE, "r");
        printf("parsing aerofile\n");
        fgets(buffer, 100, f_aero);
        fgets(buffer, 100, f_aero);
        for (i = 0; i < len; i++)
        {

            fscanf(f_aero_3, "%lf %lf %lf %lf %lf", &angle, &CL[i], &CD[i], &CL_a[i], &CD_a[i]);
            gsl_vector_set(aoa_diff, i, fabs(rad2deg(pcase->a_eff_nw) - angle));
        }
        idx = gsl_vector_min_index(aoa_diff);
        fclose(f_aero_3);
        printf("\nidx is %d\n", idx);
        pstate->CD_nw = CD[idx];
        pstate->CL_nw = CL[idx];
        pstate->CL_nw = CL_a[idx];
        pstate->CD_a_nw = CD_a[idx];

        FILE *f_aero_4 = fopen(AEROFILE, "r");
        printf("parsing aerofile\n");
        fgets(buffer, 100, f_aero);
        fgets(buffer, 100, f_aero);
        for (i = 0; i < len; i++)
        {

            fscanf(f_aero_4, "%lf %lf %lf %lf %lf", &angle, &CL[i], &CD[i], &CL_a[i], &CD_a[i]);
            gsl_vector_set(aoa_diff, i, fabs(rad2deg(pcase->a_eff_pw) - angle));
        }
        idx = gsl_vector_min_index(aoa_diff);
        fclose(f_aero_4);
        printf("\nidx is %d\n", idx);
        pstate->CD_pw = CD[idx];
        pstate->CL_pw = CL[idx];
        pstate->CL_pw = CL_a[idx];
        pstate->CD_a_pw = CD_a[idx];
    }
    else
    {
        FILE *f_aero_5 = fopen(AEROFILE, "r");
        fscanf(f_aero_5, "%d", &len);

        aoa_diff = gsl_vector_calloc(len);
        double CL[len], CD[len], CL_a[len], CD_a[len];
        printf("parsing aerofile\n");
        char buffer[100];
        fgets(buffer, 100, f_aero_5);
        fgets(buffer, 100, f_aero_5);
        for (i = 0; i < len; i++)
        {

            fscanf(f_aero_5, "%lf %lf %lf %lf %lf", &angle, &CL[i], &CD[i], &CL_a[i], &CD_a[i]);
            gsl_vector_set(aoa_diff, i, fabs(rad2deg(pcase->a_eff_0) - angle));
        }
        idx = gsl_vector_min_index(aoa_diff);
        fclose(f_aero_5);

        printf("\nidx is %d\n", idx);
        pstate->CL = CL[idx];
        pstate->CD = CD[idx];
        pstate->CL_a = CL_a[idx];
        pstate->CD_a = CD_a[idx];
    }

    gsl_vector_free(aoa_diff);
}
void initial_conditions(case_par *pcase_p, foil *pfoil, theodorsen *ptheo, double v_0, double v_dot_0, double v_dot_dot_0, double w_0, double w_dot_0, double w_dot_dot_0, double u_inf, double incr)
{
    gsl_vector_set(pcase_p->u_0, 0, v_0);
    gsl_vector_set(pcase_p->u_0, 1, w_0);
    gsl_vector_set(pcase_p->u_dot_0, 0, v_dot_0);
    gsl_vector_set(pcase_p->u_dot_0, 1, w_dot_0);
    gsl_vector_set(pcase_p->u_dot_dot_0, 0, v_dot_dot_0);
    gsl_vector_set(pcase_p->u_dot_dot_0, 1, w_dot_dot_0);
    pcase_p->inc = incr;
    pcase_p->u_inf = u_inf;
    pcase_p->a_eff_0 = -atan((pcase_p->u_inf * sin(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 1)) / (-pcase_p->u_inf * cos(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 0))) - pfoil->theta_t;
    pcase_p->w_eff_0 = sqrt((pow((pcase_p->u_inf * sin(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 1)), 2) + pow((-pcase_p->u_inf * cos(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 0)), 2)));
    pcase_p->a_eff_pu = -atan((pcase_p->u_inf * sin(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 1)) / (-pcase_p->u_inf * cos(pcase_p->alpha + pfoil->theta_t) - (gsl_vector_get(pcase_p->u_dot_0, 0) - pcase_p->inc))) - pfoil->theta_t;
    pcase_p->w_eff_pu = sqrt((pow((pcase_p->u_inf * sin(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 1)), 2) + pow((-pcase_p->u_inf * cos(pcase_p->alpha + pfoil->theta_t) - (gsl_vector_get(pcase_p->u_dot_0, 0) - pcase_p->inc)), 2)));
    pcase_p->a_eff_nu = -atan((pcase_p->u_inf * sin(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 1)) / (-pcase_p->u_inf * cos(pcase_p->alpha + pfoil->theta_t) - (gsl_vector_get(pcase_p->u_dot_0, 0) + pcase_p->inc))) - pfoil->theta_t;
    pcase_p->w_eff_nu = sqrt((pow((pcase_p->u_inf * sin(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 1)), 2) + pow((-pcase_p->u_inf * cos(pcase_p->alpha + pfoil->theta_t) - (gsl_vector_get(pcase_p->u_dot_0, 0) + pcase_p->inc)), 2)));
    pcase_p->a_eff_pw = -atan((pcase_p->u_inf * sin(pcase_p->alpha + pfoil->theta_t) - (gsl_vector_get(pcase_p->u_dot_0, 1) - pcase_p->inc)) / (-pcase_p->u_inf * cos(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 0))) - pfoil->theta_t;
    pcase_p->w_eff_pw = sqrt((pow((pcase_p->u_inf * sin(pcase_p->alpha + pfoil->theta_t) - (gsl_vector_get(pcase_p->u_dot_0, 1) - pcase_p->inc)), 2) + pow((-pcase_p->u_inf * cos(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 0)), 2)));
    pcase_p->a_eff_nw = -atan((pcase_p->u_inf * sin(pcase_p->alpha + pfoil->theta_t) - (gsl_vector_get(pcase_p->u_dot_0, 1) + pcase_p->inc)) / (-pcase_p->u_inf * cos(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 0))) - pfoil->theta_t;
    pcase_p->w_eff_nw = sqrt((pow((pcase_p->u_inf * sin(pcase_p->alpha + pfoil->theta_t) - (gsl_vector_get(pcase_p->u_dot_0, 1) + pcase_p->inc)), 2) + pow((-pcase_p->u_inf * cos(pcase_p->alpha + pfoil->theta_t) - gsl_vector_get(pcase_p->u_dot_0, 0)), 2)));
    if (pcase_p->is_steady == false)
    {
        gsl_vector_set(ptheo->aeff, 0, pcase_p->a_eff_0);
        gsl_vector_set(ptheo->y1, 0, 0.);
        gsl_vector_set(ptheo->y2, 0, 0.);
    }
}
void newmark_integration(state_space *pstate, newmark *pnew, case_par *pcase_p)
{
    double gamma = 0.5;
    double beta = 0.25;
    int i, j, k;

    if (pstate->is_alloc == false)
    {
        printf("\nunallocated variables\n");
    }
    else if (pstate->is_init == false)
    {
        printf("\nuninitialized conditions\n");
    }
    else
    {
        gsl_vector *u_0 = gsl_vector_calloc(cols);
        gsl_vector *u_dot_0 = gsl_vector_calloc(cols);
        gsl_vector *u_dot_dot_0 = gsl_vector_calloc(cols);
        gsl_vector_set(u_0, 0, gsl_vector_get(pcase_p->u_0, 0));
        gsl_vector_set(u_dot_0, 0, gsl_vector_get(pcase_p->u_dot_0, 0));
        gsl_vector_set(u_dot_dot_0, 0, gsl_vector_get(pcase_p->u_dot_dot_0, 0));
        gsl_vector_set(u_0, 1, gsl_vector_get(pcase_p->u_0, 1));
        gsl_vector_set(u_dot_0, 1, gsl_vector_get(pcase_p->u_dot_0, 1));
        gsl_vector_set(u_dot_dot_0, 1, gsl_vector_get(pcase_p->u_dot_dot_0, 1));

        gsl_matrix_set_row(pstate->u, 0, pcase_p->u_0);
        gsl_matrix_set_row(pstate->u_dot, 0, pcase_p->u_dot_0);
        gsl_matrix_set_row(pstate->u_dot_dot, 0, pcase_p->u_dot_dot_0);

        gsl_matrix *Cgdt = gsl_matrix_calloc(cols, cols);
        gsl_matrix *Kbdt2 = gsl_matrix_calloc(cols, cols);
        for (j = 0; j < cols; j++)
        {
            for (k = 0; k < cols; k++)
            {
                gsl_matrix_set(Cgdt, j, k, gsl_matrix_get(pstate->C, j, k));
                gsl_matrix_set(Kbdt2, j, k, gsl_matrix_get(pstate->K, j, k));
            };
        }

        gsl_matrix_scale(Cgdt, gamma * pcase_p->dt);
        gsl_matrix_scale(Kbdt2, beta * pow(pcase_p->dt, 2));

        gsl_matrix_add(Cgdt, pstate->M);
        gsl_matrix_add(Kbdt2, Cgdt);
        printf("\n Kb_01 = %g\n", gsl_matrix_get(Kbdt2, 0, 1));
        for (j = 0; j < cols; j++)
        {
            for (k = 0; k < cols; k++)
            {
                gsl_matrix_set(pnew->A, j, k, gsl_matrix_get(Kbdt2, j, k));
            }
        }

        gsl_matrix_free(Cgdt);
        gsl_matrix_free(Kbdt2);
        // Perform LU decomposition
        int signum;
        gsl_matrix *A_fac = gsl_matrix_calloc(2, 2);

        for (j = 0; j < cols; j++)
        {
            for (k = 0; k < cols; k++)
            {
                gsl_matrix_set(A_fac, j, k, gsl_matrix_get(pnew->A, j, k));
            };
        }
        // Define the permutation matrix
        gsl_permutation *p = gsl_permutation_alloc(2);
        gsl_linalg_LU_decomp(A_fac, p, &signum);

        // Compute the inverse of the matrix
        gsl_matrix *A_inv = gsl_matrix_calloc(2, 2);
        gsl_linalg_LU_invert(A_fac, p, A_inv);
        gsl_permutation_free(p);
        gsl_matrix_free(A_fac);
        for (j = 0; j < cols; j++)
        {
            for (k = 0; k < cols; k++)
            {
                printf("\nC %d,%d = %g\n", j, k, gsl_matrix_get(pstate->C, j, k));
            };
        }
        FILE *ftr = fopen(TRAJFILE, "a");
        fprintf(ftr, "\nv,w\n");
        fclose(ftr);

        for (i = 1; i < pcase_p->ntime; i++)
        {
            gsl_vector *u0dt = gsl_vector_calloc(cols);
            gsl_vector *udbdt = gsl_vector_calloc(cols);

            gsl_vector_set(u0dt, 0, gsl_vector_get(u_dot_0, 0));
            gsl_vector_set(u0dt, 1, gsl_vector_get(u_dot_0, 1));
            gsl_vector_set(udbdt, 0, gsl_vector_get(u_dot_dot_0, 0));
            gsl_vector_set(udbdt, 1, gsl_vector_get(u_dot_dot_0, 1));

            gsl_vector_scale(u0dt, pcase_p->dt);
            gsl_vector_scale(udbdt, (0.5 - beta) * (pow(pcase_p->dt, 2)));
            gsl_vector_add(u0dt, u_0);
            gsl_vector_add(udbdt, u0dt);
            gsl_vector_set(pnew->pred, 0, gsl_vector_get(udbdt, 0));
            gsl_vector_set(pnew->pred, 1, gsl_vector_get(udbdt, 1));

            gsl_vector_free(u0dt);
            gsl_vector_free(udbdt);

            gsl_vector *udgt = gsl_vector_calloc(cols);
            gsl_vector_set(udgt, 0, gsl_vector_get(u_dot_dot_0, 0));
            gsl_vector_set(udgt, 1, gsl_vector_get(u_dot_dot_0, 1));
            gsl_vector_scale(udgt, (1 - gamma) * pcase_p->dt);
            gsl_vector_add(udgt, u_dot_0);

            gsl_vector_set(pnew->pred_dot, 0, gsl_vector_get(udgt, 0));
            gsl_vector_set(pnew->pred_dot, 1, gsl_vector_get(udgt, 1));
            gsl_vector_free(udgt);

            gsl_vector *Q_temp = gsl_vector_calloc(2);
            gsl_vector_set(Q_temp, 0, gsl_vector_get(pstate->Q, 0));
            gsl_vector_set(Q_temp, 1, gsl_vector_get(pstate->Q, 1));

            gsl_vector *Kp = gsl_vector_calloc(2);
            gsl_blas_dgemv(CblasNoTrans, 1.0, pstate->K, pnew->pred, 0.0, Kp);
            // printf("\nB 0 = %g\n", gsl_vector_get(Kp, 0));
            // printf("\nB 1 = %g\n", gsl_vector_get(Kp, 1));

            gsl_vector *Cpd = gsl_vector_calloc(cols);
            gsl_blas_dgemv(CblasNoTrans, 1.0, pstate->C, pnew->pred_dot, 0.0, Cpd);

            gsl_vector_sub(Q_temp, Kp);
            gsl_vector_sub(Q_temp, Cpd);
            ;

            gsl_vector_set(pnew->B, 0, gsl_vector_get(Q_temp, 0));
            gsl_vector_set(pnew->B, 1, gsl_vector_get(Q_temp, 1));
            // printf("\nB 0 = %g\n",  gsl_vector_get(pnew->B, 0));
            // printf("\nB 1 = %g\n", gsl_vector_get(pnew->B, 1));

            gsl_vector_free(Q_temp);
            gsl_vector_free(Kp);
            gsl_vector_free(Cpd);

            gsl_blas_dgemv(CblasNoTrans, 1.0, A_inv, pnew->B, 0.0, u_dot_dot_0);
            // printf("\nB 0 = %g\n",  gsl_vector_get(u_dot_dot_0, 0));
            // printf("\nB 1 = %g\n", gsl_vector_get(u_dot_dot_0, 1));

            gsl_matrix_set_row(pstate->u_dot_dot, i, u_dot_dot_0);

            gsl_vector *temp_u = gsl_vector_calloc(cols);
            gsl_vector *temp_ud = gsl_vector_calloc(cols);

            gsl_vector_set(temp_u, 0, gsl_vector_get(pnew->pred, 0));
            gsl_vector_set(temp_u, 1, gsl_vector_get(pnew->pred, 1));

            gsl_vector_set(temp_ud, 0, gsl_vector_get(u_dot_dot_0, 0));
            gsl_vector_set(temp_ud, 1, gsl_vector_get(u_dot_dot_0, 1));

            gsl_vector_scale(temp_ud, (beta) * (pow(pcase_p->dt, 2)));
            // printf("\nB 0 = %g\n",  gsl_vector_get(temp_ud, 0));
            // printf("\nB 1 = %g\n", gsl_vector_get(temp_ud, 1));
            gsl_vector_add(temp_u, temp_ud);

            gsl_vector_set(u_0, 0, gsl_vector_get(temp_u, 0));
            gsl_vector_set(u_0, 1, gsl_vector_get(temp_u, 1));
            gsl_matrix_set_row(pstate->u, i, u_0);

            gsl_vector_free(temp_u);
            gsl_vector_free(temp_ud);

            gsl_vector *temp_u_dot = gsl_vector_calloc(cols);
            gsl_vector *temp_udd = gsl_vector_calloc(cols);

            gsl_vector_set(temp_u_dot, 0, gsl_vector_get(pnew->pred_dot, 0));
            gsl_vector_set(temp_u_dot, 1, gsl_vector_get(pnew->pred_dot, 1));

            gsl_vector_set(temp_udd, 0, gsl_vector_get(u_dot_dot_0, 0));
            gsl_vector_set(temp_udd, 1, gsl_vector_get(u_dot_dot_0, 1));
            gsl_vector_scale(temp_udd, (gamma)*pcase_p->dt);
            gsl_vector_add(temp_u_dot, temp_udd);
            gsl_vector_set(u_dot_0, 0, gsl_vector_get(temp_u_dot, 0));
            gsl_vector_set(u_dot_0, 1, gsl_vector_get(temp_u_dot, 1));
            gsl_matrix_set_row(pstate->u_dot, i, u_dot_0);
            FILE *ftr = fopen(TRAJFILE, "a");
            fprintf(ftr, "\n%g,%g\n", gsl_matrix_get(pstate->u, i, 0), gsl_matrix_get(pstate->u, i, 1));
            fclose(ftr);
            gsl_vector_free(temp_u_dot);
            gsl_vector_free(temp_udd);
        }
        gsl_vector_free(u_0);
        gsl_vector_free(u_dot_0);
        gsl_vector_free(u_dot_dot_0);
        gsl_matrix_free(A_inv);
    }
}

void runge_kutta(state_space *pstate, case_par *pcase_p, theodorsen *ptheo, foil *pfoil)
{
    FILE *ftr = fopen(TRAJFILE, "a");
    fprintf(ftr, "\nv,w\n");
    fclose(ftr);
    int i, j, l;
    double rk[4] = {0.1084, 0.2602, 0.5052, 1.};

    double u_0[2];
    double u_dot_0[2];
    double u_dot_dot_0[2];
    // gsl_vector *Ku = gsl_matrix_calloc(cols,cols);
    u_0[0] = gsl_vector_get(pcase_p->u_0, 0);
    u_dot_0[0] = gsl_vector_get(pcase_p->u_dot_0, 0);
    u_dot_dot_0[0] = gsl_vector_get(pcase_p->u_dot_dot_0, 0);
    u_0[1] = gsl_vector_get(pcase_p->u_0, 1);
    u_dot_0[1] = gsl_vector_get(pcase_p->u_dot_0, 1);
    u_dot_dot_0[1] = gsl_vector_get(pcase_p->u_dot_dot_0, 1);

    double y1_0 =0.1;
    double y2_0 = -0.;

    gsl_matrix_set_row(pstate->u, 0, pcase_p->u_0);
    gsl_matrix_set_row(pstate->u_dot, 0, pcase_p->u_dot_0);
    gsl_matrix_set_row(pstate->u_dot_dot, 0, pcase_p->u_dot_dot_0);

    double u_dot_temp[2];
    double u_temp[2];
    double u_dot_dot_temp[2];
    u_dot_temp[0] = 0;
    u_dot_temp[1] = 0;
    u_temp[0] = 0;
    u_temp[1] = 0;
    u_dot_dot_temp[0] = 0;
    u_dot_dot_temp[1] = 0;

    double y1_dot_temp =0;
    double y2_dot_temp =0;
    double y1_temp =0.;
    double y2_temp =0.;
    double a_e_temp = 0;
    double a_eff_temp =0;
    find_aero(pstate, pcase_p);
    double CL_t_temp = pstate->CL;
    double CD_t_temp = pstate->CD;
    double Fx;
    double Fz;
    double L;
    double D;
    double W_eff = pcase_p->w_eff_0;

    gsl_vector_set(ptheo->CL_t, 0, pstate->CL);
    gsl_vector_set(ptheo->CD_t, 0, pstate->CD);
    double A_1 = 0.165;
    double A_2 = 0.335;
    double b_1 = 0.0455;
    double b_2 = 0.3;
    //  A 0.165, 0.335, b 0.0455, b 0.3000      .
    for (i = 1; i < pcase_p->ntime; i++)

    {
        for (l = 0; l < 4; l++)
        {
            a_eff_temp = pcase_p->a_eff_0 - (u_dot_temp[1] / W_eff);
            a_e_temp = a_eff_temp * (1 - A_1 - A_2) + y1_temp + y2_temp;
            W_eff = sqrt(pow(-pcase_p->w_eff_0*cos(a_eff_temp)-u_dot_temp[0],2)+pow(pcase_p->w_eff_0*sin(a_eff_temp)-u_dot_temp[1],2));
            CL_t_temp = 2 * M_PI * (a_e_temp - 0.0314) - u_dot_dot_temp[1] * M_PI / (2 * pow(W_eff, 2));
            printf("\n CL_t_temp = %lf",CL_t_temp);
            printf("\n W_eff = %lf",W_eff);
            L = CL_t_temp * 0.5 * 1.225 * pow(W_eff, 2);
            D = CD_t_temp * 0.5 * 1.225 * pow(W_eff, 2);

            Fx = L * sin(a_eff_temp + pfoil->theta_t) - D * cos(a_eff_temp + pfoil->theta_t);
            Fz = L * cos(a_eff_temp + pfoil->theta_t) + D * sin(a_eff_temp + pfoil->theta_t);

            u_dot_dot_temp[0] = (Fx - gsl_matrix_get(pstate->K, 0, 0) * u_temp[0] - gsl_matrix_get(pstate->K, 0, 1) * u_temp[1]) / pfoil->mass;
            u_dot_dot_temp[1] = (Fz - gsl_matrix_get(pstate->K, 1, 0) * u_temp[0] - gsl_matrix_get(pstate->K, 1, 1) * u_temp[1]) / pfoil->mass;

            y1_dot_temp = -b_1 * 2 * W_eff * y1_temp + b_1 * A_1 * 2 * W_eff * a_eff_temp;
            y2_dot_temp = -b_2 * 2 * W_eff * y2_temp + b_2 * A_2 * 2 * W_eff * a_eff_temp;

            y1_temp = y1_0 + rk[l] * pcase_p->dt * y1_dot_temp;
            y2_temp = y2_0 + rk[l] * pcase_p->dt * y2_dot_temp;

            u_temp[0] = u_0[0] + rk[l] *pcase_p->dt* u_dot_temp[0];
            u_temp[1] = u_0[1] + rk[l] * pcase_p->dt * u_dot_temp[1];

            u_dot_temp[0] = u_dot_0[0] + rk[l] * pcase_p->dt * u_dot_dot_temp[0];
            u_dot_temp[1] = u_dot_0[1] + rk[l] * pcase_p->dt * u_dot_dot_temp[1];
        }
        printf("\ntime_step\n");
        y1_0 = y1_temp;
        y2_0 = y2_temp;

        u_dot_0[0] = u_dot_temp[0];
        u_dot_0[1] = u_dot_temp[1];

        u_0[0] = u_temp[0];
        u_0[1] = u_temp[1];

        gsl_matrix_set(pstate->u, i, 0, u_0[0]);
        gsl_matrix_set(pstate->u, i, 1, u_0[1]);
        gsl_matrix_set(pstate->u_dot, i, 0, u_dot_0[0]);
        gsl_matrix_set(pstate->u_dot, i, 1, u_dot_0[1]);
        gsl_matrix_set(pstate->u_dot_dot, i, 0, u_dot_dot_temp[0]);
        gsl_matrix_set(pstate->u_dot_dot, i, 1, u_dot_dot_temp[1]);
        FILE *ftr = fopen(TRAJFILE, "a");
        fprintf(ftr, "\n%g,%g\n", gsl_matrix_get(pstate->u, i, 0), gsl_matrix_get(pstate->u, i, 1));
        fclose(ftr);
    }
}