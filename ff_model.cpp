/**************************************************************************
* Copyright (C) 2011,2013-2014 Dr. Bogdan Tanygin<b.m.tanygin@gmail.com>
* All rights reserved.
* 
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*************************************************************************/

//system inclusions
///////////////////
#include <windows.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

//#include <Eigen/Dense>
//using namespace std;
//using namespace Eigen;

#include "ff_model.h"
#include "ff_model_parameters.h"
#include "ff_model_graphics.h"
#include "ff_model_io.h"

// working variables 
////////////////////
boost::mt19937 rng;
boost::normal_distribution<>* nd; // TODO need delete
boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >* var_nor;
boost::uniform_01<>* ud;
boost::variate_generator<boost::mt19937&, boost::uniform_01<> >* var_uni;

ff_vect_t r[pN + 1];  //particles positions
double Rp_to_c[pN + 1];
ff_vect_t m[pN + 1];  //particles magnetic moment direction
ff_vect_t mt[pN + 1];
ff_vect_t m_prev[pN + 1];
ff_vect_t m_prev_before_r[pN + 1]; // magnetic moment direction before random motion
int m_sat[pN + 1];

ff_vect_t F[pN + 1]; // mean force
ff_vect_t F1[pN + 1]; // force at the t = t
ff_vect_t F2[pN + 1]; // force at the t = t + dt
//ff_vect_t P[pN + 1];
ff_vect_t tau[pN + 1]; // mean torque
ff_vect_t tau1[pN + 1]; // torque at the t = t
ff_vect_t tau2[pN + 1]; // // torque at the t = t + dt
//ff_vect_t tau_r[pN + 1]; // random torque

ff_vect_t v[pN + 1];
//ff_vect_t v_r[pN + 1]; // random heat component
ff_vect_t w[pN + 1]; // angular velocity vector
//ff_vect_t w_r[pN + 1]; // random heat component

ff_vect_t drt[pN + 1];
ff_vect_t drt_r[pN + 1]; // instantiated random translation
ff_vect_t dvt[pN + 1];
//ff_vect_t dvt_r[pN + 1]; // instantiated random velocity
//double dv_r[pN + 1]; // extra random velocity magnitude
ff_vect_t dphi[pN + 1];
ff_vect_t dphi_r[pN + 1]; // instantiated random rotation
//ff_vect_t dm[pN + 1];
ff_vect_t dw[pN + 1];
//double dw_r[pN + 1]; // extra random angular velocity magnitude

//ff_vect_t dir110[13];

long i_min = 1;
double V0_tot = 0; // total volume of the dispersed phase
double V0_tot_oleic = 0; // total volume of the dispersed phase inside the oleic droplet
double V0_largest_EV = 0; // mathematical expected value of largest particles total volume // see is_large_mode variable
double V0_tot_EV = 0; // mathematical expected value of particles total volume
double I_glob = 0;

int exist_p[pN + 1]; // particle existence; number of primary aggregate inside
int is_neel[pN + 1]; // Neel relaxation
int is_temp_sat[pN + 1]; // temperature saturation flag
int is_inside_oleic[pN + 1];
//int aggregated_p[pN + 1][pN + 1]; // map of particles aggregation, in case of dW > G_barrier
double Rp0[pN + 1];
double Rp[pN + 1];
double Vp0[pN + 1];
double Vpfull[pN + 1]; // including steric layer
double m0p[pN + 1];
double M0p[pN + 1];
double I0p[pN + 1]; // particle moment of inertia
double r0modp[pN + 1];
//double Vselfp[pN + 1];
double C2[pN + 1];
double gamma_rot[pN + 1];

int time_go = 0;

double B_hyst[21];
double Mz_hyst[21];
double B_hyst_n[21];
double Mz_hyst_n[21];

double dt = 0;
double t = 0; // time
double t_temp = 0; // temperature, [C]
double dT = 0;
double dT_prev = 0;
double T_basic = 0;
double T_mean = 0;
double T_mean_loc = 0;
double T_mean_loc_prev = 0;
double T_mean_loc_prev_revert = 0;
long k_mean = 0;

double dT_p_x[pN + 1];
double dT_p_y[pN + 1];
double dT_p_z[pN + 1];
double dT_p_rot_x[pN + 1];
double dT_p_rot_y[pN + 1];
double dT_p_rot_z[pN + 1];
double dT_prev_p_x[pN + 1];
double dT_prev_p_y[pN + 1];
double dT_prev_p_z[pN + 1];
double dT_prev_p_rot_x[pN + 1];
double dT_prev_p_rot_y[pN + 1];
double dT_prev_p_rot_z[pN + 1];
double T_basic_p_x[pN + 1];
double T_basic_p_y[pN + 1];
double T_basic_p_z[pN + 1];
double T_basic_p_rot_x[pN + 1];
double T_basic_p_rot_y[pN + 1];
double T_basic_p_rot_z[pN + 1];
double T_mean_p_x[pN + 1];
double T_mean_p_y[pN + 1];
double T_mean_p_z[pN + 1];
double T_mean_p_rot_x[pN + 1];
double T_mean_p_rot_y[pN + 1];
double T_mean_p_rot_z[pN + 1];
double T_mean_loc_p_x[pN + 1];
double T_mean_loc_p_y[pN + 1];
double T_mean_loc_p_z[pN + 1];
double T_mean_loc_p_rot_x[pN + 1];
double T_mean_loc_p_rot_y[pN + 1];
double T_mean_loc_p_rot_z[pN + 1];
//double T_mean_loc_prev_p[pN + 1];
//double T_mean_loc_prev_revert_p[pN + 1];
long k_mean_p[pN + 1];

double Ek = 0;
double Ek_rot = 0;
double Ek_tr = 0;
double Ekp_x[pN + 1];
double Ekp_y[pN + 1];
double Ekp_z[pN + 1];
double Ekp_rot_x[pN + 1];
double Ekp_rot_y[pN + 1];
double Ekp_rot_z[pN + 1];
double G_dd[pN + 1]; // dipole-dipole energy
//double dW[pN + 1]; // work
double g_Bz_prev;
long step = 0;

double kB = 0;
int hyst_mode = 1;
double mz_tot;
long glob_start_step = 1;
long glob_start_step_susc = 1;
double glob_start_t;
double mz_glob = 0; // global mean average start from the hyst. point switch
ff_vect_t m_tot_glob;

//double g_theta[pN + 1], g_phi[pN + 1];

int g_hyst_start_line;
int g_hyst_up_line;
int g_hyst_bottom_line;

double BmanX = 0;
double BmanY = 0;
double BmanZ = 0;

ff_vect_t m_tot;

ff_vect_t dir110[13];

double k_force_adapt_p_0[pN + 1];

double k_force_adapt;
double k_force_adapt_p_x[pN + 1];
double k_force_adapt_p_y[pN + 1];
double k_force_adapt_p_z[pN + 1];
double k_force_adapt_p_rot_x[pN + 1];
double k_force_adapt_p_rot_y[pN + 1];
double k_force_adapt_p_rot_z[pN + 1];
double k_force_adapt_mean = 0;
double k_force_adapt_mean_print = 0;

ff_vect_t r_brown_valid_0;

long pN_oleic_drop = 0;
long pN_oleic_drop_I = 0;
long pN_oleic_drop_II = 0;
long pN_oleic_drop_III = 0;
double d[14 + 1];

double dt_red = 0; // reducing time indicator for the random translation
double k_delta_force_rel[pN + 1]; // relation of random force to attractive forces sum
double k_delta_force_rel_tot = 0;
double k_delta_force_rel_p = 0;
double k_delta_torque_rel[pN + 1];
double k_delta_torque_rel_tot = 0;
double k_delta_torque_rel_p = 0;

double phi_vol_fract_oleic = 0;
double phi_vol_fract_oleic_0 = 0;
double R_oleic;
long pN0 = pN;
double R0_min;

void ff_model_upgrade_ext_field(void)
{
    g_Bz_prev = Bext(0,0,0).z;
    switch (hyst_mode)
    {
    case 1:
        kB = 0;
        g_hyst_start_line = 1;
        g_hyst_up_line = 0;
        g_hyst_bottom_line = 0;
        break;

    case 2:
        kB = 0.25;
        g_hyst_start_line = 1;
        g_hyst_up_line = 0;
        g_hyst_bottom_line = 0;
        break;

    case 3:
        kB = 0.5;
        g_hyst_start_line = 1;
        g_hyst_up_line = 0;
        g_hyst_bottom_line = 0;
        break;

    case 4:
        kB = 0.75;
        g_hyst_start_line = 1;
        g_hyst_up_line = 0;
        g_hyst_bottom_line = 0;
        break;

    case 5:
        kB = 1;
        g_hyst_start_line = 1;
        g_hyst_up_line = 1;
        g_hyst_bottom_line = 0;
        break;

    case 6:
        kB = 0.75;
        g_hyst_start_line = 0;
        g_hyst_up_line = 1;
        g_hyst_bottom_line = 0;
        break;

    case 7:
        kB = 0.5;
        g_hyst_start_line = 0;
        g_hyst_up_line = 1;
        g_hyst_bottom_line = 0;
        break;

    case 8:
        kB = 0.25;
        g_hyst_start_line = 0;
        g_hyst_up_line = 1;
        g_hyst_bottom_line = 0;
        break;

    case 9:
        kB = 0;
        g_hyst_start_line = 0;
        g_hyst_up_line = 1;
        g_hyst_bottom_line = 0;
        break;

    case 10:
        kB = -0.25;
        g_hyst_start_line = 0;
        g_hyst_up_line = 1;
        g_hyst_bottom_line = 0;
        break;

    case 11:
        kB = -0.5;
        g_hyst_start_line = 0;
        g_hyst_up_line = 1;
        g_hyst_bottom_line = 0;
        break;

    case 12:
        kB = -0.75;
        g_hyst_start_line = 0;
        g_hyst_up_line = 1;
        g_hyst_bottom_line = 0;
        break;

    case 13:
        kB = -1;
        g_hyst_start_line = 0;
        g_hyst_up_line = 1;
        g_hyst_bottom_line = 1;
        break;

    case 14:
        kB = -0.75;
        g_hyst_start_line = 0;
        g_hyst_up_line = 0;
        g_hyst_bottom_line = 1;
        break;

    case 15:
        kB = -0.5;
        g_hyst_start_line = 0;
        g_hyst_up_line = 0;
        g_hyst_bottom_line = 1;
        break;

    case 16:
        kB = -0.25;
        g_hyst_start_line = 0;
        g_hyst_up_line = 0;
        g_hyst_bottom_line = 1;
        break;

    case 17:
        kB = 0;
        g_hyst_start_line = 0;
        g_hyst_up_line = 0;
        g_hyst_bottom_line = 1;
        break;

    case 18:
        kB = 0.25;
        g_hyst_start_line = 0;
        g_hyst_up_line = 0;
        g_hyst_bottom_line = 1;
        break;

    case 19:
        kB = 0.5;
        g_hyst_start_line = 0;
        g_hyst_up_line = 0;
        g_hyst_bottom_line = 1;
        break;

    case 20:
        kB = 0.75;
        g_hyst_start_line = 0;
        g_hyst_up_line = 0;
        g_hyst_bottom_line = 1;
        break;
    }
}

void ff_model_auto_hyst(void)
{
    int tmp;

    if ((glob_start_step == 1) && (t - glob_start_t >= 0))
        glob_start_step = step;


    if ((t - glob_start_t) >= 0.25 /*0.25 of the 1/4.0 of period*/ * 0.25 * (1 / nu_ext))
    {

        hyst_mode++;
        tmp = 0;
        if (hyst_mode == 21) {hyst_mode = 5; tmp = 1;}
        ff_model_upgrade_ext_field();
        if (tmp) {g_hyst_up_line = 1;g_hyst_start_line = 0;}

        //       ff_io_save_hyst();

        mz_glob = 0;
        glob_start_t = t;
        glob_start_step = step;

    }
}

int ff_model_check_smooth_dr(long p)
{
    double dr, rmod, dphimag;
    long ps;
    int res = 1;
    /*ff_vect_t dr_tmp;
    ff_vect_t dphi_tmp;

    dr_tmp.x = drt[p].x + drt_r[p].x;
    dr_tmp.y = drt[p].y + drt_r[p].y;
    dr_tmp.z = drt[p].z + drt_r[p].z;

    dphi_tmp.x = dphi[p].x + dphi_r[p].x;
    dphi_tmp.y = dphi[p].y + dphi_r[p].y;
    dphi_tmp.z = dphi[p].z + dphi_r[p].z;*/
    
    dr = sqrt(MUL(drt[p], drt[p]));
    //dr = sqrt(MUL(dr_tmp, dr_tmp));
    //rmod = d[4];
    rmod = 2 * Rp[p];
    //rmod = delta;
    //rmod = 2 * R0_min;
    //rmod = 0.5 * delta_r;

    dphimag = sqrt(MUL(dphi[p], dphi[p]));
    //dphimag = sqrt(MUL(dphi_tmp, dphi_tmp));

    if (rmod > 0)
        if ((dr / rmod > smooth_r) || ((dphimag / pi > smooth_r) && (!(is_neel[p]))))
        {
            //if ((dr / rmod > smooth_r)) printf("\n DEBUG SMOOTH dr = %e", dr);
            //if ((dt < 1E-15) && (dr / rmod > smooth_r)) printf("\n DEBUG SMOOTH dr = %e, (dr / rmod > smooth_r), k_force_adapt_p[p] = %e, Fx[p] = %e, Fy[p] = %e, Fz[p] = %e, Px[p] = %e, Py[p] = %e, Pz[p] = %e, vx[p] = %e, vy[p] = %e, vz[p] = %e", dr, k_force_adapt_p[p], F[p].x, F[p].y, F[p].z, P[p].x, P[p].y, P[p].z, v[p].x, v[p].y, v[p].z);
            //if ((dt < 1E-15) && (dphimag / pi > smooth_r)) printf("\n DEBUG SMOOTH dr = %e, (dphimag / pi > smooth_r), taux[p] = %e, tauy[p] = %e, tauz[p] = %e, taux_r[p] = %e, tauy_r[p] = %e, tauz_r[p] = %e, wx[p] = %e, wy[p] = %e, wz[p] = %e", dr, tau[p].x, tau[p].y, tau[p].z, tau_r[p].x, tau_r[p].y, tau_r[p].z, w[p].x, w[p].y, w[p].z);
            
            for(ps = 1; ps <= p; ps ++)
            {
                r[ps].x -= drt[ps].x;
                r[ps].y -= drt[ps].y;
                r[ps].z -= drt[ps].z;

                /*m[ps].x -= dm[ps].x;
                m[ps].y -= dm[ps].y;
                m[ps].z -= dm[ps].z;*/

                if ((!(is_neel[ps])) && (ps < p))
                {
                    m[ps] = m_prev[ps];
                }
            }

            for(ps = 1; ps <= pN; ps ++)
            {
                r[ps].x -= drt_r[ps].x;
                r[ps].y -= drt_r[ps].y;
                r[ps].z -= drt_r[ps].z;

                if (!(is_neel[ps]))
                {
                    m[ps] = m_prev[ps];
                }
            }

            dt /= 2.0;
            slow_steps += 10;
            res = 0;
            step--;

            T_mean -= T_basic;
            k_mean --;
            T_mean_loc -= T_basic;
            k_bm_inst --;
            if (k_bm_inst == k_bm_inst_max - 1)
            {
                /*if (dT <= 0) k_force_adapt *= k_force_adapt_0;
                else k_force_adapt /= k_force_adapt_0;*/

                dT= dT_prev;
                T_mean_loc_prev = T_mean_loc_prev_revert;
            }

            for(ps = 1; ps <= pN; ps ++)
            {
                T_mean_p_x[ps] -= T_basic_p_x[ps];
                T_mean_p_y[ps] -= T_basic_p_y[ps];
                T_mean_p_z[ps] -= T_basic_p_z[ps];
                T_mean_p_rot_x[ps] -= T_basic_p_rot_x[ps];
                T_mean_p_rot_y[ps] -= T_basic_p_rot_y[ps];
                T_mean_p_rot_z[ps] -= T_basic_p_rot_z[ps];
                k_mean_p[ps] --;
                T_mean_loc_p_x[ps] -= T_basic_p_x[ps];
                T_mean_loc_p_y[ps] -= T_basic_p_y[ps];
                T_mean_loc_p_z[ps] -= T_basic_p_z[ps];
                T_mean_loc_p_rot_x[ps] -= T_basic_p_rot_x[ps];
                T_mean_loc_p_rot_y[ps] -= T_basic_p_rot_y[ps];
                T_mean_loc_p_rot_z[ps] -= T_basic_p_rot_z[ps];
                //k_bm_inst --;
                if (k_bm_inst == k_bm_inst_max - 1)
                {
                    if (dT_p_x[ps] <= 0) k_force_adapt_p_x[ps] *= k_force_adapt_0;
                    else k_force_adapt_p_x[ps] /= k_force_adapt_0;
                    if (dT_p_y[ps] <= 0) k_force_adapt_p_y[ps] *= k_force_adapt_0;
                    else k_force_adapt_p_y[ps] /= k_force_adapt_0;
                    if (dT_p_z[ps] <= 0) k_force_adapt_p_z[ps] *= k_force_adapt_0;
                    else k_force_adapt_p_z[ps] /= k_force_adapt_0;
                    if (dT_p_rot_x[ps] <= 0) k_force_adapt_p_rot_x[ps] *= k_force_adapt_0;
                    else k_force_adapt_p_rot_x[ps] /= k_force_adapt_0;
                    if (dT_p_rot_y[ps] <= 0) k_force_adapt_p_rot_y[ps] *= k_force_adapt_0;
                    else k_force_adapt_p_rot_y[ps] /= k_force_adapt_0;
                    if (dT_p_rot_z[ps] <= 0) k_force_adapt_p_rot_z[ps] *= k_force_adapt_0;
                    else k_force_adapt_p_rot_z[ps] /= k_force_adapt_0;

                    dT_p_x[ps] = dT_prev_p_x[ps];
                    dT_p_y[ps] = dT_prev_p_y[ps];
                    dT_p_z[ps] = dT_prev_p_z[ps];
                    dT_p_rot_x[ps] = dT_prev_p_rot_x[ps];
                    dT_p_rot_y[ps] = dT_prev_p_rot_y[ps];
                    dT_p_rot_z[ps] = dT_prev_p_rot_z[ps];

                    //T_mean_loc_prev_p_x[ps] = T_mean_loc_prev_revert_p[ps];
                }
            }
        }

        return res;
}

/*
int ff_model_check_smooth_dv(long p)
{
double dv, vmod, dm, mmag;
long ps;
int res = 1;
ff_vect_t v_prev;

v_prev.x = v[p].x - dvt[p].x;
v_prev.y = v[p].y - dvt[p].y;
v_prev.z = v[p].z - dvt[p].z;

dv = sqrt(MUL(dvt[p], dvt[p]));
vmod = sqrt(MUL(v_prev, v_prev));

//printf("\n !!! %e", vmod * dt / R0);

if (vmod > 0)
if (dv / vmod > smooth_v)
{
for(ps = 1; ps <= p; ps ++)
{
v[ps].x -= dvt[ps].x;
v[ps].y -= dvt[ps].y;
v[ps].z -= dvt[ps].z;

m[p].x = m_prev[p].x;
m[p].y = m_prev[p].y;
m[p].z = m_prev[p].z;

F[p] = F0[p];
}

for(ps = 1; ps <= pN; ps ++) // ALL particles because dr was first in the flow
{
r[ps].x -= drt[ps].x;
r[ps].y -= drt[ps].y;
r[ps].z -= drt[ps].z;
}

dt /= 2.0;
slow_steps += 10;
res = 0;
}

return res;
}
*/

/*double ff_model_distrib_f(double theta, long p)
{
double ret_f;
ff_vect_t B1; // excluding self field
double Bmag;
double Ct; 
double dtheta = pi / 10.0; // hardcode

B1.x = B[p].x - Bself * m[p].x / m0;
B1.y = B[p].y - Bself * m[p].y / m0;
B1.z = B[p].z - Bself * m[p].z / m0;

Bmag = sqrt(MUL(B1,B1));

Ct = m0 * Bmag / (kb * T);

ret_f = (Ct / (exp(Ct) - exp(-Ct))) * exp(Ct * cos(theta)) * sin(theta) * dtheta;

return ret_f;
}*/

/*double ff_model_distrib_f_r(double theta, long p)
{
double ret_f;
double G0;
double Ct; 
double dtheta = pi / 10.0; // hardcode

G0 = sqrt(MUL(F[p],F[p])) * r0mod * sqrt(dt);
Ct = G0 / (kb * T);

ret_f = (Ct / (exp(Ct) - exp(-Ct))) * exp(Ct * cos(theta)) * sin(theta) * dtheta;

return ret_f;
}*/

/*
void __deprecated__ff_model_set_rand_dir(long p)
{
long i;
double dtheta; // step in [R]andom long
double R_points[10 + 1];
long R;
double theta, phi;
double Fmag, tmag, bmag;
ff_vect_t Ft, t, b;

/*
dtheta = pi / 10.0;
for (i = 1; i <= 10; i++)
{
R_points[i] = 32768.0 * ff_model_distrib_f_r((i - 0.5) * dtheta, p);
if (i > 1) R_points[i] += R_points[i - 1];
//printf("\n i = %d, Rpoints = %f", i, R_points[i]);
}

R = rand();

theta = 0.5 * dtheta; // default for zero
for (i = 1; i <= 9; i++)
if ((R > R_points[i])&&(R <= R_points[i + 1]))
{
theta = (i - 0.5) * dtheta;
break;
}
*//*
phi = 2 * pi * rand() / 32768.0;
theta = pi * rand() / 32768.0;

Fmag = sqrt(MUL(F[p],F[p]));

if (Fmag == 0) Fmag = 0.00000000001;

Ft = F[p];
if ((F[p].x == 0)&&(F[p].y == 0)) {Ft.x = 1E-5; Ft.y = 1E-5;}

t.x = - Ft.y;
t.y =   Ft.x;
t.z =   0;

tmag = sqrt(MUL(t,t));

b.x = t.y * Ft.z - t.z * Ft.y;
b.y = t.z * Ft.x - t.x * Ft.z;
b.z = t.x * Ft.y - t.y * Ft.x;

bmag = sqrt(MUL(b,b));

drt[p].x += r0modp[p] * sqrt(dt) * (Ft.x * cos(theta) / Fmag + t.x * sin(theta) * cos(phi) / tmag + b.x * sin(theta) * sin(phi) / bmag);
drt[p].y += r0modp[p] * sqrt(dt) * (Ft.y * cos(theta) / Fmag + t.y * sin(theta) * cos(phi) / tmag + b.y * sin(theta) * sin(phi) / bmag);
drt[p].z += r0modp[p] * sqrt(dt) * (Ft.z * cos(theta) / Fmag + t.z * sin(theta) * cos(phi) / tmag + b.z * sin(theta) * sin(phi) / bmag);
}*/

/*void ff_model_set_m(long p)
{
long i;
double dtheta; // step in [R]andom long
double R_points[10 + 1];
long R;
double theta, phi;
double Bmag, tmag, bmag;
ff_vect_t Bt, t, b;

dtheta = pi / 10.0;
for (i = 1; i <= 10; i++)
{
R_points[i] = 32768.0 * ff_model_distrib_f((i - 0.5) * dtheta, p);
if (i > 1) R_points[i] += R_points[i - 1];
//printf("\n i = %d, Rpoints = %f", i, R_points[i]);
}

R = rand();

theta = 0.5 * dtheta; // default for zero
for (i = 1; i <= 9; i++)
if ((R > R_points[i])&&(R <= R_points[i + 1]))
{
theta = (i - 0.5) * dtheta;
break;
}
phi = 2 * pi * rand() / 32768.0;

Bmag = sqrt(MUL(B[p],B[p]));

if (Bmag == 0) Bmag = 0.00000000001;

Bt = B[p];
if ((B[p].x == 0)&&(B[p].y == 0)) {Bt.x = 1E-10; Bt.y = 1E-10;}

t.x = - Bt.y;
t.y =   Bt.x;
t.z =   0;

tmag = sqrt(MUL(t,t));

b.x = t.y * Bt.z - t.z * Bt.y;
b.y = t.z * Bt.x - t.x * Bt.z;
b.z = t.x * Bt.y - t.y * Bt.x;

bmag = sqrt(MUL(b,b));

m[p].x = m0 * (Bt.x * cos(theta) / Bmag + t.x * sin(theta) * cos(phi) / tmag + b.x * sin(theta) * sin(phi) / bmag);
m[p].y = m0 * (Bt.y * cos(theta) / Bmag + t.y * sin(theta) * cos(phi) / tmag + b.y * sin(theta) * sin(phi) / bmag);
m[p].z = m0 * (Bt.z * cos(theta) / Bmag + t.z * sin(theta) * cos(phi) / tmag + b.z * sin(theta) * sin(phi) / bmag);
}*/

ff_vect_t ff_model_nonloc_torque(long p)
{
    //long p;
    register long ps;
    double mxs, mys, mzs;
    double dx, dy, dz;
    double dx_q_plus, dy_q_plus, dz_q_plus;
    double dx_q_minus, dy_q_minus, dz_q_minus, dR_q_plus, dR_q_minus;
    double eps, max_eps;

    double MUL2mod__dR5mod1, dR2__dR5mod1;
    double dR, dR2, dR5, dR5mod1, MUL2mod;

    double tBx, tBy, tBz;
    double dtBx, dtBy, dtBz;    // field of single ps-th particle
    double dtBqx, dtBqy, dtBqz;
    double tBmag, tm, tms;
    double tk;

    double esx, esy, esz;

    double mag_qs;

    long p_prev;

    double sign;
    int is_last;

    ff_vect_t tBext;
    ff_vect_t ttau;

    ttau.x = ttau.y = ttau.z = 0;

    int di, dj, dk;

    //	do
    //	{
    //for(p = 1; p <= pN; p ++)
    if (exist_p[p])
    {
        tBext = Bext(r[p].x, r[p].y, r[p].z);
        tBx = tBext.x;
        tBy = tBext.y;
        tBz = tBext.z;

        for(ps = 1; ps <= pN; ps ++)
            if (exist_p[ps])
                if (p != ps)
                {
                    //#ifdef SECONDARY
                    //for(di = -1; di <= 1; di++)
                    //for(dj = -1; dj <= 1; dj++)
                    //for(dk = -1; dk <= 1; dk++)
                    {
                        //#endif*/
                        dx = - r[ps].x + r[p].x;// + di * Lx;
                        dy = - r[ps].y + r[p].y;// + dj * Ly;
                        dz = - r[ps].z + r[p].z;// + dk * Lz;

                        if (is_periodic)
						{
							if (dx > Lx / 2.0) dx -= Lx;
							if (dx < - Lx / 2.0) dx += Lx;
							
							if (dy > Ly / 2.0) dy -= Ly;
							if (dy < - Ly / 2.0) dy += Ly;
							
							if (dz > Lz / 2.0) dz -= Lz;
							if (dz < - Lz / 2.0) dz += Lz;
						}

                        dR = sqrt(dx * dx + dy * dy + dz * dz);

                        dR2 = dR * dR;

                        if (1)
                        {

                            mxs = m[ps].x;
                            mys = m[ps].y;
                            mzs = m[ps].z;

                            MUL2mod = 3 * (mxs * dx + mys * dy + mzs * dz);

                            dR5 = pow(dR,5);
                            dR5mod1 = dR5 / C5;

                            MUL2mod__dR5mod1 = MUL2mod / dR5mod1;
                            dR2__dR5mod1 = dR2 / dR5mod1;

                            dtBx = dx * MUL2mod__dR5mod1 - mxs * dR2__dR5mod1;
                            dtBy = dy * MUL2mod__dR5mod1 - mys * dR2__dR5mod1;
                            dtBz = dz * MUL2mod__dR5mod1 - mzs * dR2__dR5mod1;

                            if (dtBx != dtBx) printf("\n DEBUG 12 p = %d ps = %d dtBx = %e r[p].x = %e r[ps].x = %e mxs = %e mys = %e mzs = %e dR = %e", p, ps, dtBx, r[p].x, r[ps].x, mxs, mys, mzs, dR);

                        } 
                        else
                        {

                        }

                        /*tms = sqrt(MUL(m[ps],m[ps]));
                        esx = m[ps].x / tms;
                        esy = m[ps].y / tms;
                        esz = m[ps].z / tms;

                        dx_q_plus  = dx - esx * Rp[ps]; // "plus" is mag. charge
                        dx_q_minus = dx + esx * Rp[ps];

                        dy_q_plus  = dy - esy * Rp[ps];
                        dy_q_minus = dy + esy * Rp[ps];

                        dz_q_plus  = dz - esz * Rp[ps];
                        dz_q_minus = dz + esz * Rp[ps];

                        dR_q_plus = sqrt(dx_q_plus * dx_q_plus + dy_q_plus * dy_q_plus + dz_q_plus * dz_q_plus);
                        dR_q_minus = sqrt(dx_q_minus * dx_q_minus + dy_q_minus * dy_q_minus + dz_q_minus * dz_q_minus);

                        mag_qs = tms / (4.0 * Rp[ps] / 3.0);

                        dtBqx = C5 * mag_qs * (dx_q_plus / pow(dR_q_plus, 3) - dx_q_minus / pow(dR_q_minus, 3));
                        dtBqy = C5 * mag_qs * (dy_q_plus / pow(dR_q_plus, 3) - dy_q_minus / pow(dR_q_minus, 3));
                        dtBqz = C5 * mag_qs * (dz_q_plus / pow(dR_q_plus, 3) - dz_q_minus / pow(dR_q_minus, 3));

                        tBx += dtBqx * tms / m0p[ps] + dtBx * (1 - tms / m0p[ps]);
                        tBy += dtBqy * tms / m0p[ps] + dtBy * (1 - tms / m0p[ps]);
                        tBz += dtBqz * tms / m0p[ps] + dtBz * (1 - tms / m0p[ps]);*/

                        if (dR >= Rp0[p] + Rp0[ps]) // soft-sphere model related correction
                        {
                            tBx += dtBx;
                            tBy += dtBy;
                            tBz += dtBz;
                        }

                        //#ifdef SECONDARY
                    } // di and dj
                    //#endif*/
                } // if (p != ps)

                tBmag = sqrt(tBx * tBx + tBy * tBy + tBz * tBz);
                if (is_uniform_field_test) tBx = tBy = tBz = 0;

                /*tm = sqrt(MUL(m[p],m[p]));

                eps = acos((m[p].x * tBx + m[p].y * tBy + m[p].z * tBz) / (tBmag * tm));

                if (p == 1) max_eps = eps;
                else if (eps > max_eps) max_eps = eps;*/

#ifndef SECONDARY
                /*m[p].x = m0p[p] * tBx / tBmag;
                m[p].y = m0p[p] * tBy / tBmag;
                m[p].z = m0p[p] * tBz / tBmag;*/

                if (tBmag)
				{
					ttau.x =   m[p].y * tBz - m[p].z * tBy;
					ttau.y = - m[p].x * tBz + m[p].z * tBx;
					ttau.z =   m[p].x * tBy - m[p].y * tBx;
				}

                if ((is_neel[p]) && (!is_uniform_field_test) && (tBmag))
                {
                    //printf("\n p = %d !!!", p);
					m[p].x = m0p[p] * tBx / tBmag;
                    m[p].y = m0p[p] * tBy / tBmag;
                    m[p].z = m0p[p] * tBz / tBmag;
                }

				G_dd[p] = - (m[p].x * (tBx - tBext.x) + m[p].y * (tBy - tBext.y) + m[p].z * (tBz - tBext.z));
#endif

    }// if (exist_p[p])
    //	}
    //	while(max_eps > m_h_eff_tol);
    return ttau;
}

ff_vect_t ff_model_nonloc_force(long p)
{
    register long ps;
    double tFx, tFy, tFz;
    double dtFx, dtFy, dtFz;
    double dreptFx, dreptFy, dreptFz; // repulsion
    double dr, drtemp, l, dr2, dr3, dr2mod, dr5, dr5mod, MUL1, MUL2, MUL3;
    //double R1 = 0.3 * R0;
    double Cmod;
    //double Cmod1 = Ch * m0 * m0;
    //double Cmod1_repulsion = 5 * m0 * m0;
    double dx, dy, dz;
    double tk;
    double sec_pow;
    ff_vect_t ttF, dBmaggrad;

    double MUL1__dr5mod, MUL2__dr5mod, MUL3__dr5mod, MUL1_MUL2__dr2mod__dr5mod;

    double mx, my, mz;
    double mxs, mys, mzs;
    long p_prev;

    int is_last;
    double sign, tms, mag_qs;

    //ff_vect_t ttau;

    double dd, tt, dr_min, F_steric_mag;

    //long di, dj;

    tFx = 0;
    tFy = 0;
    tFz = 0;

    for(ps = 1; ps <= pN; ps ++)
        if (exist_p[ps])
            /*#ifdef SECONDARY
            for(di = -2; di <= 2; di++)
            for(dj = -2; dj <= 2; dj++) // including zero!
            {
            #endif*/
            if (p != ps)
            {
                dr_min = Rp0[p] + Rp0[ps] + delta_r;

                dx = r[ps].x - r[p].x;// + di * Lx;
                dy = r[ps].y - r[p].y;// + dj * Ly;
                dz = r[ps].z - r[p].z;

                if (is_periodic)
				{
					if (dx > Lx / 2.0) dx -= Lx;
					if (dx < - Lx / 2.0) dx += Lx;
					
					if (dy > Ly / 2.0) dy -= Ly;
					if (dy < - Ly / 2.0) dy += Ly;
					
					if (dz > Lz / 2.0) dz -= Lz;
					if (dz < - Lz / 2.0) dz += Lz;
				}

                mx = m[p].x;
                my = m[p].y;
                mz = m[p].z;

                mxs = m[ps].x;
                mys = m[ps].y;
                mzs = m[ps].z;

                dr = sqrt(dr2 = dx * dx + dy * dy + dz * dz); 
                dr3 = dr2 * dr;

                //if (dr == 0) printf ("\n \n \n dr == 0 !!!");

                // Non local force and self-confid. magnetic field
                dr2mod = dr2 / 5.0; //modified
                dr5mod = (dr5 = pow(dr,5)) / C1; //modified

                MUL2 = mxs * dx + mys * dy + mzs * dz;

                //#ifndef SECONDARY
                //if (p < ps)  !!!!!
                //{

                {
                    //#endif
                    MUL1 = mx * dx + my * dy + mz * dz;
                    MUL3 = mx * mxs + my * mys + mz * mzs;

                    MUL1__dr5mod = MUL1 / dr5mod;
                    MUL2__dr5mod = MUL2 / dr5mod;
                    MUL3__dr5mod = MUL3 / dr5mod;
                    MUL1_MUL2__dr2mod__dr5mod = MUL1 * MUL2 / (dr2mod * dr5mod);

                    dtFx = -(MUL1__dr5mod * mxs + MUL2__dr5mod * mx
                        + MUL3__dr5mod * dx - MUL1_MUL2__dr2mod__dr5mod * dx);

                    dtFy = -(MUL1__dr5mod * mys + MUL2__dr5mod * my
                        + MUL3__dr5mod * dy - MUL1_MUL2__dr2mod__dr5mod * dy);

                    dtFz = -(MUL1__dr5mod * mzs + MUL2__dr5mod * mz
                        + MUL3__dr5mod * dz - MUL1_MUL2__dr2mod__dr5mod * dz);

                    //if (p == 50) printf("\n DEBUG 7 ps = %d F.x = %e", ps, dtFx);
                    //if (p == 50) printf("\n DEBUG 7 ps = %d m[p].x = %e", ps, m[p].x);
                    //if (p == 50) printf("\n DEBUG 7 ps = %d m[ps].x = %e", ps, m[ps].x);

                } 

                if (dr >= dr_min)
                {
                    tFx += dtFx;
                    tFy += dtFy;
                    tFz += dtFz;
                }

                /*
                if (dr <= 5 * R0) // exp phenomenology
                {
                tFx += - (EPS * dx / dr) * (exp(-(dr - 2 * R00) / ro1) / ro1 - exp(-(dr - 2 * R00) / ro2) / ro2);
                tFy += - (EPS * dy / dr) * (exp(-(dr - 2 * R00) / ro1) / ro1 - exp(-(dr - 2 * R00) / ro2) / ro2);
                tFz += - (EPS * dz / dr) * (exp(-(dr - 2 * R00) / ro1) / ro1 - exp(-(dr - 2 * R00) / ro2) / ro2);
                }
                */

//#ifndef SECONDARY
                if ((dr < dr_min) && (!isMicroDrop)) //soft sphere condition
                {
                    Cmod = 10 * m0p[p] * m0p[ps] * (C1 / dr5);
					
                    if (dr)
					{
						tFx += - dx * Cmod;
						tFy += - dy * Cmod;
						tFz += - dz * Cmod;
					}

                    /*if ((dW[p] > G_barrier) && (!aggregated_p[p][ps]))
                    {
                    aggregated_p[p][ps] = 1;
                    //printf("\n aggregated_p");
                    }*/

                    /*drtemp = Rp0[p] + Rp0[ps];
                    dd = Rp0[p] + Rp0[ps];
                    l = 2 * (drtemp - dd) / dd;
                    tt = 2 * delta / dd;

                    tFx += - Ch_ss * (dx / drtemp) * (2 * pow(dd, 2) * kb * T * N_o * pi * log((tt + 1) / (l / 2 + 1)) / tt);
                    tFy += - Ch_ss * (dy / drtemp) * (2 * pow(dd, 2) * kb * T * N_o * pi * log((tt + 1) / (l / 2 + 1)) / tt);
                    tFz += - Ch_ss * (dz / drtemp) * (2 * pow(dd, 2) * kb * T * N_o * pi * log((tt + 1) / (l / 2 + 1)) / tt);*/
                }

                // Entropic repulsion
                if ((dr >= dr_min) && (dr <= Rp0[p] + Rp0[ps] + 2 * delta)) 
                {
                    /*dd = Rp0[p] + Rp0[ps];
                    l = 2 * (dr - dd) / dd;
                    tt = 2 * delta / dd;
                    tFx += - (dx / dr) * (2 * pow(dd, 2) * kb * T * N_o * pi * log((tt + 1) / (l / 2 + 1)) / tt);
                    tFy += - (dy / dr) * (2 * pow(dd, 2) * kb * T * N_o * pi * log((tt + 1) / (l / 2 + 1)) / tt);
                    tFz += - (dz / dr) * (2 * pow(dd, 2) * kb * T * N_o * pi * log((tt + 1) / (l / 2 + 1)) / tt);*/

                    //T=273;
					F_steric_mag = ((-pi)*kb*T*N_o*(dr-2*delta-Rp0[ps]-Rp0[p])*((Rp0[ps]+Rp0[p])*(2*(dr*dr+delta*delta)+(-pow(Rp0[ps],2)+Rp0[p]*Rp0[ps]-pow(Rp0[p],2))*(dr/(Rp0[ps]+Rp0[p])+1)+delta*dr)+(2*Rp0[p]*Rp0[ps]-pow((Rp0[p]-Rp0[ps]),2))*delta))/(6*delta*dr*dr);
                    //T=0;

                    //printf("\n F_steric_mag = %e", F_steric_mag);

                    tFx += - (dx / dr) * F_steric_mag;
                    tFy += - (dy / dr) * F_steric_mag;
                    tFz += - (dz / dr) * F_steric_mag;
                }
//#endif

//#ifndef SECONDARY
                // attraction
                //if ((dr > Rp[p] + Rp[ps] )&&(dr < 3 * (Rp[p] + Rp[ps]) / 2.0 )) // the Heaviside step function  and dr5 dependence finally is similar to the well-known exp. phenomenology
                //if (dr > Rp[p] + Rp[ps] + 2 * smooth_r * delta)
                //if (dr > (Rp[p] + Rp[ps]) * (1 + smooth_r))
                //if (dr > (Rp0[p] + Rp0[ps] + a0))
                //if ((dr > (Rp0[p] + Rp0[ps] + 2 * delta)) && (!isMicroDrop))
                if (dr >= dr_min)// && (!isMicroDrop))
                {
                    /*Cmod = Ch * m0p[p] * m0p[ps] * (C1 / dr5);

                    tFx += dx * Cmod;
                    tFy += dy * Cmod;
                    tFz += dz * Cmod;*/

                    tFx += - (dx / dr) * (- 64 * A_H * dr * pow(Rp0[p] * Rp0[ps], 3) / (6 * pow(dr * dr - pow(Rp0[p] - Rp0[ps], 2), 2) * pow(dr * dr - pow(Rp0[p] + Rp0[ps], 2), 2)));
                    tFy += - (dy / dr) * (- 64 * A_H * dr * pow(Rp0[p] * Rp0[ps], 3) / (6 * pow(dr * dr - pow(Rp0[p] - Rp0[ps], 2), 2) * pow(dr * dr - pow(Rp0[p] + Rp0[ps], 2), 2)));
                    tFz += - (dz / dr) * (- 64 * A_H * dr * pow(Rp0[p] * Rp0[ps], 3) / (6 * pow(dr * dr - pow(Rp0[p] - Rp0[ps], 2), 2) * pow(dr * dr - pow(Rp0[p] + Rp0[ps], 2), 2)));

                    /*if (aggregated_p[p][ps])
                    {
                    Cmod = 5 * m0p[p] * m0p[ps] * (C1 / dr5);

                    tFx += dx * Cmod;
                    tFy += dy * Cmod;
                    tFz += dz * Cmod;
                    }*/
                }
//#endif

                /*Fnonloc[p][ps].x = tFx;
                Fnonloc[p][ps].y = tFy;
                Fnonloc[p][ps].z = tFz;*/

                //#ifndef SECONDARY
                //} // end of if p < ps
                /*else
                {
                tFx -= Fnonloc[ps][p].x; // third Newton law
                tFy -= Fnonloc[ps][p].y;
                tFz -= Fnonloc[ps][p].z;
                }*/ // end of if p > ps
                //#endif
                //tBx += dx * MUL2mod__dr5mod1 - mxs * dr2__dr5mod1;
                //tBy += dy * MUL2mod__dr5mod1 - mys * dr2__dr5mod1;
                //tBz += dz * MUL2mod__dr5mod1 - mzs * dr2__dr5mod1;
                /*#ifdef SECONDARY
                }
                #endif*/
            } // end of the particles loop


            //tBmag = sqrt(tBx * tBx + tBy * tBy + tBz * tBz);

            // here we should add self magnetic field inside the particle
            //B[p].x = tBx + Bself * m[p].x / m0;
            //B[p].y = tBy + Bself * m[p].y / m0;
            //B[p].z = tBz + Bself * m[p].z / m0;

            /*m_prev[p].x = m[p].x;
            m_prev[p].y = m[p].y;
            m_prev[p].z = m[p].z;*/

            //m[p].x = m0 * tBx / tBmag;
            //m[p].y = m0 * tBy / tBmag;
            //m[p].z = m0 * tBz / tBmag;
            //ff_model_set_m(p);

            if (ext_field_is_homo == 0)
            {
                dBmaggrad = dBz_het_bem(r[p]);
                tFx += m[p].z * dBmaggrad.x;
                tFy += m[p].z * dBmaggrad.y;
                tFz += m[p].z * dBmaggrad.z;
            }

            // force
            ttF.x = tFx;
            ttF.y = tFy;
            ttF.z = tFz;

            //if (ttF.x != ttF.x) printf("\n DEBUG 10 p = %d F.x = %e", p, ttF.x);
            //if (ttF.y != ttF.y) printf("\n DEBUG 10 p = %d F.y = %e", p, ttF.y);
            //if (ttF.z != ttF.z) printf("\n DEBUG 10 p = %d F.z = %e", p, ttF.z);

            return ttF;
}

double ff_model_G_zeeman(long p)
{
	return - MUL(Bext(r[p].x, r[p].y, r[p].z), m[p]);
}

double ff_model_G_steric(long p, long ps)
{
    double G_steric = 0;
    double z = 0;
    double dx = 0; 
    double dy = 0; 
    double dz = 0;
    double dR = 0;
    double dRmax = 0;
	double dr_min = 0;

	dr_min = Rp0[p] + Rp0[ps] + delta_r;

    dx = r[p].x - r[ps].x;
    dy = r[p].y - r[ps].y;
    dz = r[p].z - r[ps].z;
    dR = sqrt(dx * dx + dy * dy + dz * dz);

    dRmax = Rp0[ps] + Rp0[p] + 2 * delta;

    // Entropic repulsion
    if ((dR <= dRmax) && (dR >= dr_min))
        G_steric = (pi * kb * pow(dR - (2 * delta + Rp0[ps] + Rp0[p]), 2) * ((Rp0[ps] + Rp0[p]) * 
        (dR + delta) - (pow(Rp0[ps], 3) + pow(Rp0[p], 3)) / (Rp0[p] + Rp0[ps])) * N_o * T) / ( 6 * delta * dR);
    else G_steric = 0;

    return G_steric;
}

double ff_model_G_london(long p, long ps)
{
    double G_london = 0;
    double z = 0;
    double dx = 0; 
    double dy = 0; 
    double dz = 0;
    double dR = 0;
    double dRmax = 0;
	double dr_min = 0;

	dr_min = Rp0[p] + Rp0[ps] + delta_r;

    dx = r[p].x - r[ps].x;
    dy = r[p].y - r[ps].y;
    dz = r[p].z - r[ps].z;
    dR = sqrt(dx * dx + dy * dy + dz * dz);

    //dRmax = Rp0[ps] + Rp0[p] + 2 * delta;

    // London dispersion energy
    if ((1 /*dR <= dRmax*/) && (dR >= dr_min))
        G_london = (- A_H / 6.0) * (2 * Rp0[p] * Rp0[ps] / (dR * dR - pow(Rp0[p] + Rp0[ps], 2)) + 2 * Rp0[p] * Rp0[ps] / (dR * dR - pow(Rp0[p] - Rp0[ps], 2))
		                            + log((dR * dR - pow(Rp0[p] + Rp0[ps], 2)) / (dR * dR - pow(Rp0[p] - Rp0[ps], 2))));
    else G_london = 0;

    return G_london;
}

ff_vect_t ff_model_force(long p)
{
    ff_vect_t tF;
    ff_vect_t tddF;
    double F_nonloc_mag = 0; // nonlocal force magnitude
    double R_mag = 0; // random force magnitude

    tF.x = tF.y = tF.z = 0;

    // non-local
    tddF = ff_model_nonloc_force(p);
    tF.x += tddF.x;
    tF.y += tddF.y;
    tF.z += tddF.z;

    F_nonloc_mag = sqrt(MUL(tddF,tddF));

    // Gravitation
    //tF.z += - C3;

    // Buoyancy force
    //tF.z +=   C6;

    // oleic droplet surface tension force
    if ((fabs(Rp_to_c[p] - R_oleic) < Rp0[p]) && (is_oleic) && (R_oleic > Rp0[p]) && (!isPGasPrevail))
    {
        tF.x += - sigma_sf * 2 * pi * Rp0[p] * r[p].x / Rp_to_c[p];
        tF.y += - sigma_sf * 2 * pi * Rp0[p] * r[p].y / Rp_to_c[p];
        tF.z += - sigma_sf * 2 * pi * Rp0[p] * r[p].z / Rp_to_c[p];
    }

    /*tF.x += P[p].x;
    tF.y += P[p].y;
    tF.z += P[p].z;*/

    //R_mag = sqrt(MUL(P[p],P[p]));

    /*k_delta_force_rel[p] = R_mag / F_nonloc_mag;
    k_delta_force_rel_tot += k_delta_force_rel[p];
    k_delta_force_rel_p ++;*/

    if (is_uniform_field_test)
    {
        tF.x = tF.y = tF.z = 0;
        tF.z += - 1E0 * sigma_sf * 2 * pi * 5E-9;
    }

    return tF;
}

ff_vect_t ff_model_torque(long p)
{
    ff_vect_t ttau;
    ff_vect_t dtau;
    double ttau_nonloc_mag = 0; // nonlocal torque magnitude
    double tau_r_mag = 0; // random torque magnitude

    ttau.x = ttau.y = ttau.z = 0;

    // non-local
    dtau = ff_model_nonloc_torque(p);

    ttau_nonloc_mag = sqrt(MUL(dtau,dtau));

    if (!(is_neel[p]))
    {
        ttau.x += dtau.x;
        ttau.y += dtau.y;
        ttau.z += dtau.z;

        // Gravitation
        ttau.z += 0;

        // Buoyancy
        ttau.z += 0;
    }

    //printf("\n %e", ttau.x);

    /*ttau.x += tau_r[p].x;
    ttau.y += tau_r[p].y;
    ttau.z += tau_r[p].z;

    tau_r_mag = sqrt(MUL(tau_r[p],tau_r[p]));
    if (!(is_neel[p]))
    {
        k_delta_torque_rel[p] = tau_r_mag / ttau_nonloc_mag;
        k_delta_torque_rel_tot += k_delta_torque_rel[p];
        k_delta_torque_rel_p ++;
    }*/

    return ttau;
}

void ff_model_next_step(void)
{ 
    ff_vect_t f, ttau;
    double I, I_min_per;
	double I_per[3][3][3];
    long p;
    int chk;
    long mz_tot_n;
    ff_vect_t r0;
	ff_vect_t r0_per[3][3][3];
	ff_vect_t r_per[pN + 1][3][3][3];
    //double C2, gamma_rot;
    double M0, I0;
    ff_vect_t ex, ey, ez; // basis of rotation around ez = w
    double tmmag;
    double theta_0_r, phi_0_r;
    double Mtot;
    double t_temp_1 = 0;
    double V_oleic = 0;

	double per_x_plus, per_x_minus, per_y_plus, per_y_minus, per_z_plus, per_z_minus;

    Ek = Ek_tr = Ek_rot = 0;
    mz_tot = 0;
    m_tot.x = m_tot.y = m_tot.z = 0;
    mz_tot_n = 0;
    pN_oleic_drop = 0;
    pN_oleic_drop_I = pN_oleic_drop_II = pN_oleic_drop_III = 0;
    phi_vol_fract_oleic = 0;
    V0_tot_oleic = 0;
    k_delta_force_rel_tot = 0;
    k_delta_force_rel_p = 0;
    P_pgas = 0;

    t_temp = T + ta0;
    if (t_temp > 90) t_temp_1 = 90;
    if (t_temp < 20) t_temp_1 = 20;
    ro_oleic = 902.0 - 0.62 * T;
    sigma_sf = sigma_sf_nano * (a_sigma_sf + b_sigma_sf * t_temp_1);
    eta_oleic = a3_eta_oleic * pow(t_temp_1, 3) + a2_eta_oleic * pow(t_temp_1, 2) + a1_eta_oleic * pow(t_temp_1, 1) + a0_eta_oleic;

    if (isMicroDrop) ff_model_update_mdrop_parameters();

    for (p = 1; p <= pN; p++) if (exist_p[p])
    {
        Ek_tr  += M0p[p] * MUL(v[p],v[p]) / 2.0;
        Ek_rot += I0p[p] * MUL(w[p],w[p]) / 2.0;

        Ekp_x[p] = M0p[p] * v[p].x * v[p].x / 2.0;
        Ekp_y[p] = M0p[p] * v[p].y * v[p].y / 2.0;
        Ekp_z[p] = M0p[p] * v[p].z * v[p].z / 2.0;
        Ekp_rot_x[p] = I0p[p] * w[p].x * w[p].x / 2.0;
        Ekp_rot_y[p] = I0p[p] * w[p].y * w[p].y / 2.0;
        Ekp_rot_z[p] = I0p[p] * w[p].z * w[p].z / 2.0;
    };

    Ek = Ek_tr + Ek_rot;

    if (time_go)
    {
        step++;
        //printf("\n ============================", step);
        //printf("\n Step %d", step);

        //ff_model_m_setting();

        ff_model_update_dT();
        for (p = 1; p <= pN; p++) if (exist_p[p]) ff_model_update_dT_p(p);

        for (p = 1; p <= pN; p++) if (exist_p[p]) ff_model_effective_random_motion_update(p);
        if (dt_red <= 0)
        {
            dt_red = dt0;
            //for (p = 1; p <= pN; p++) if (exist_p[p]) ff_model_effective_random_motion_update(p);
        }

        for (p = 1; p <= pN; p++)
        if (exist_p[p])
        {
            //ff_model_update_dT_p(p);

            //Rp_to_c[p] = sqrt(MUL(r[p], r[p]));

            /*if (k_bm_inst == 1)*/
            
            //ff_model_effective_random_force_update(p);
            
            f = ff_model_force(p);
            ttau = ff_model_torque(p);
            F1[p].x = f.x; F1[p].y = f.y; F1[p].z = f.z;
            tau1[p].x = ttau.x; tau1[p].y = ttau.y; tau1[p].z = ttau.z;
        }

        for (p = 1; p <= pN; p++) if (exist_p[p])
        {
            r[p].x += drt_r[p].x;
            r[p].y += drt_r[p].y;
            r[p].z += drt_r[p].z;

            if (!(is_neel[p]))
            {
                dphi[p].x = dphi_r[p].x;
                dphi[p].y = dphi_r[p].y;
                dphi[p].z = dphi_r[p].z;
                
                //dm[p].x = dm[p].y = dm[p].z = 0;

                mt[p] = m[p];
                // turn around ex
                mt[p].y = mt[p].y * cos(dphi[p].x) - mt[p].z * sin(dphi[p].x);
                mt[p].z = mt[p].y * sin(dphi[p].x) + mt[p].z * cos(dphi[p].x);
                // turn around ey
                mt[p].z = mt[p].z * cos(dphi[p].y) - mt[p].x * sin(dphi[p].y);
                mt[p].x = mt[p].z * sin(dphi[p].y) + mt[p].x * cos(dphi[p].y);
                // turn around ez
                mt[p].x = mt[p].x * cos(dphi[p].z) - mt[p].y * sin(dphi[p].z);
                mt[p].y = mt[p].x * sin(dphi[p].z) + mt[p].y * cos(dphi[p].z);

                //dm[p].x = mt[p].x - m[p].x;
                //dm[p].y = mt[p].y - m[p].y;
                //dm[p].z = mt[p].z - m[p].z;

                //m[p].x += dm[p].x;
                //m[p].y += dm[p].y;
                //m[p].z += dm[p].z;

                m_prev_before_r[p] = m[p];
                m[p] = mt[p];

                tmmag = sqrt(MUL(m[p], m[p]));

                m[p].x *= m0p[p] / tmmag;
                m[p].y *= m0p[p] / tmmag;
                m[p].z *= m0p[p] / tmmag;
            }
        }
                
        for (p = 1; p <= pN; p++) if (exist_p[p])
        {
            f = ff_model_force(p);
            ttau = ff_model_torque(p);
            F2[p].x = f.x; F2[p].y = f.y; F2[p].z = f.z;
            tau2[p].x = ttau.x; tau2[p].y = ttau.y; tau2[p].z = ttau.z;

            F[p].x = (F1[p].x + F2[p].x) / 2.0;
            F[p].y = (F1[p].y + F2[p].y) / 2.0;
            F[p].z = (F1[p].z + F2[p].z) / 2.0;

            tau[p].x = (tau1[p].x + tau2[p].x) / 2.0;
            tau[p].y = (tau1[p].y + tau2[p].y) / 2.0;
            tau[p].z = (tau1[p].z + tau2[p].z) / 2.0;

            /*DEBUG*/ if (f.x != f.x) printf("\n DEBUG 1 p = %d f.x = %e", p, f.x);
            /*DEBUG*/ if (f.y != f.y) printf("\n DEBUG 1 p = %d f.y = %e", p, f.y);
            /*DEBUG*/ if (f.z != f.z) printf("\n DEBUG 1 p = %d f.z = %e", p, f.z);

            /*DEBUG*/ //if (tau[p].x != tau[p].x) printf("\n DEBUG 1 p = %d ttau.x = %e", p, tau[p].x);
            /*DEBUG*/ //if (ttau.y != ttau.y) printf("\n DEBUG 1 p = %d ttau.y = %e", p, ttau.y);
            /*DEBUG*/ //if (tau[p].z != tau[p].z) printf("\n DEBUG 1 p = %d ttau.z = %e", p, tau[p].z);

            /*DEBUG*/ if (v[p].x != v[p].x) printf("\n DEBUG 1.1 p = %d r[p].x = %e v[p].x = %e", p, r[p].x, v[p].x);
            /*DEBUG*/ if (v[p].y != v[p].y) printf("\n DEBUG 1.1 p = %d r[p].y = %e v[p].y = %e", p, r[p].y, v[p].y);
            /*DEBUG*/ if (v[p].z != v[p].z) printf("\n DEBUG 1.1 p = %d r[p].z = %e v[p].z = %e", p, r[p].z, v[p].z);

            /*DEBUG*/ if (w[p].x != w[p].x) printf("\n DEBUG 1.2 p = %d", p);
            /*DEBUG*/ if (w[p].y != w[p].y) printf("\n DEBUG 1.2 p = %d", p);
            /*DEBUG*/ if (w[p].z != w[p].z) printf("\n DEBUG 1.2 p = %d", p);
        }

            //ff_model_update_dT();
            k_bm_inst ++;

            for (p = 1; p <= pN; p++)
                if (exist_p[p])
                {
                    if ((Rp_to_c[p] > R_oleic) || (!is_oleic))
                    {
                        C2[p] = 6 * pi * eta_car * Rp[p];
                        gamma_rot[p] = 8 * pi * eta_car * pow(Rp[p], 3);
                    }
                    else
                    {
                        C2[p] = 6 * pi * eta_oleic * Rp[p];
                        gamma_rot[p] = 8 * pi * eta_oleic * pow(Rp[p], 3);
                    }

                    M0 = M0p[p];
                    I0 = I0p[p];

                    drt[p].x = F[p].x * dt / C2[p] +		 
                        (v[p].x - F[p].x / C2[p]) * (1 - exp(- C2[p] * dt / M0)) * M0 / C2[p];
                    //+ drt_r[p].x * dt / dt0;

                    drt[p].y = F[p].y * dt / C2[p] +		 
                        (v[p].y - F[p].y / C2[p]) * (1 - exp(- C2[p] * dt / M0)) * M0 / C2[p];
                    //+ drt_r[p].y * dt / dt0;

                    drt[p].z = F[p].z * dt / C2[p] +		 
                        (v[p].z - F[p].z / C2[p]) * (1 - exp(- C2[p] * dt / M0)) * M0 / C2[p];
                    //+ drt_r[p].z * dt / dt0;

                    //if (brownian_shifts) ff_model_set_rand_dir(p);

                    if (!(is_neel[p]))
                    {
                        dphi[p].x = tau[p].x * dt / gamma_rot[p] +		 
                            (w[p].x - tau[p].x / gamma_rot[p]) * (1 - exp(- gamma_rot[p] * dt / I0)) * I0 / gamma_rot[p];
                        //+ dphi_r[p].x * dt / dt0;

                        dphi[p].y = tau[p].y * dt / gamma_rot[p] +		 
                            (w[p].y - tau[p].y / gamma_rot[p]) * (1 - exp(- gamma_rot[p] * dt / I0)) * I0 / gamma_rot[p];
                        //+ dphi_r[p].y * dt / dt0;

                        dphi[p].z = tau[p].z * dt / gamma_rot[p] +		 
                            (w[p].z - tau[p].z / gamma_rot[p]) * (1 - exp(- gamma_rot[p] * dt / I0)) * I0 / gamma_rot[p];
                        //+ dphi_r[p].z * dt / dt0;
                    }
                    else dphi[p].x = dphi[p].y = dphi[p].z = 0;

                    /*DEBUG*/ if (dphi[p].x != dphi[p].x) printf("\n DEBUG 1 p = %d dphi[p].x = %e", p, dphi[p].x);
                    /*DEBUG*/ if (dphi[p].y != dphi[p].y) printf("\n DEBUG 1 p = %d dphi[p].y = %e", p, dphi[p].y);
                    /*DEBUG*/ if (dphi[p].z != dphi[p].z) printf("\n DEBUG 1 p = %d dphi[p].z = %e", p, dphi[p].z);

                    //ff_model_check_overlapp(p); // hard sphere condition

                    //printf("\n %e", dphi[p].y);

                    r[p].x += drt[p].x;
                    r[p].y += drt[p].y;
                    r[p].z += drt[p].z;

                    //if (m0p[p] == 0) printf("\n !!!");

                    //printf("\n %e", sqrt(MUL(m[p],m[p])) / m0p[p]);

                    if (r[p].x != r[p].x) printf("\n DEBUG 2 p = %d r[p].x = %e v[p].x = %e", p, r[p].x, v[p].x);
                    if (r[p].y != r[p].y) printf("\n DEBUG 2 p = %d r[p].y = %e v[p].y = %e", p, r[p].y, v[p].y);
                    if (r[p].z != r[p].z) printf("\n DEBUG 2 p = %d r[p].z = %e v[p].z = %e", p, r[p].z, v[p].z);

                    chk = ff_model_check_smooth_dr(p);
                    if ( chk == 0)
                    {
                        //k_bm_inst = 1;
                        goto t_end;
                    }

                    /*r[p].x += drt_r[p].x;
                    r[p].y += drt_r[p].y;
                    r[p].z += drt_r[p].z;

                    dphi[p].x += dphi_r[p].x;
                    dphi[p].y += dphi_r[p].y;
                    dphi[p].z += dphi_r[p].z;*/

                    //m_prev[p] = m[p];
                    if (!(is_neel[p]))
                    {
                        //dm[p].x = dm[p].y = dm[p].z = 0;

                        mt[p] = m[p];
                        // turn around ex
                        mt[p].y = mt[p].y * cos(dphi[p].x) - mt[p].z * sin(dphi[p].x);
                        mt[p].z = mt[p].y * sin(dphi[p].x) + mt[p].z * cos(dphi[p].x);
                        // turn around ey
                        mt[p].z = mt[p].z * cos(dphi[p].y) - mt[p].x * sin(dphi[p].y);
                        mt[p].x = mt[p].z * sin(dphi[p].y) + mt[p].x * cos(dphi[p].y);
                        // turn around ez
                        mt[p].x = mt[p].x * cos(dphi[p].z) - mt[p].y * sin(dphi[p].z);
                        mt[p].y = mt[p].x * sin(dphi[p].z) + mt[p].y * cos(dphi[p].z);

                        //dm[p].x = mt[p].x - m[p].x;
                        //dm[p].y = mt[p].y - m[p].y;
                        //dm[p].z = mt[p].z - m[p].z;

                        //m[p].x += dm[p].x;
                        //m[p].y += dm[p].y;
                        //m[p].z += dm[p].z;

                        m_prev[p] = m[p];
                        m[p] = mt[p];

                        tmmag = sqrt(MUL(m[p], m[p]));

                        m[p].x *= m0p[p] / tmmag;
                        m[p].y *= m0p[p] / tmmag;
                        m[p].z *= m0p[p] / tmmag;
                    }

                    ff_model_check_walls(p);
                } // end of loop for dr

                //k_bm_inst++;
                //if (k_bm_inst == k_bm_inst_max) k_bm_inst = 1;
                if (k_bm_inst == k_bm_inst_max)
                {
                    k_bm_inst = 1;
                    T_mean_loc = 0;
                    for (p = 1; p <= pN; p++) T_mean_loc_p_x[p] = T_mean_loc_p_y[p] = T_mean_loc_p_z[p] = T_mean_loc_p_rot_x[p] = T_mean_loc_p_rot_y[p] = T_mean_loc_p_rot_z[p] = 0;
                    k_force_adapt_mean /= (6 * pN);
                    k_force_adapt_mean_print = k_force_adapt_mean;
                    k_force_adapt_mean = 0;
                    ///printf("\n !!!", dT);
                }

                r0.x = r0.y = r0.z = 0;
				for (int t1 = 0; t1 <= 2; t1++)
					for (int t2 = 0; t2 <= 2; t2++)
						for (int t3 = 0; t3 <= 2; t3++)
						{
							r0_per[t1][t2][t3].x = r0_per[t1][t2][t3].y = r0_per[t1][t2][t3].z = 0;
							I_per[t1][t2][t3] = 0;
						}

				Mtot = 0;
				
                for (p = 1; p <= pN; p++)
                    if (exist_p[p])
                    {
                        /*if (Rp_to_c[p] > R_oleic)
                        {
                        C2 = 6 * pi * eta * Rp[p];
                        gamma_rot = 8 * pi * eta * pow(Rp[p], 3);
                        }
                        else
                        {
                        C2 = 6 * pi * eta_oleic * Rp[p];
                        gamma_rot = 8 * pi * eta_oleic * pow(Rp[p], 3);
                        }*/

                        M0 = M0p[p];
                        I0 = I0p[p];

                        // C2 is a friction
                        dvt[p].x = (F[p].x / C2[p]) + (v[p].x - F[p].x / C2[p]) * exp(- C2[p] * dt / M0) - v[p].x;
                        dvt[p].y = (F[p].y / C2[p]) + (v[p].y - F[p].y / C2[p]) * exp(- C2[p] * dt / M0) - v[p].y;
                        dvt[p].z = (F[p].z / C2[p]) + (v[p].z - F[p].z / C2[p]) * exp(- C2[p] * dt / M0) - v[p].z;

                        v[p].x += dvt[p].x;
                        v[p].y += dvt[p].y;
                        v[p].z += dvt[p].z;

                        w[p].x = (tau[p].x / gamma_rot[p]) + (w[p].x - tau[p].x / gamma_rot[p]) * exp(- gamma_rot[p] * dt / I0);
                        w[p].y = (tau[p].y / gamma_rot[p]) + (w[p].y - tau[p].y / gamma_rot[p]) * exp(- gamma_rot[p] * dt / I0);
                        w[p].z = (tau[p].z / gamma_rot[p]) + (w[p].z - tau[p].z / gamma_rot[p]) * exp(- gamma_rot[p] * dt / I0);

                        if (v[p].x != v[p].x) printf("\n DEBUG 3 p = %d v[p].x = %e", p, v[p].x);
                        if (v[p].y != v[p].y) printf("\n DEBUG 3 p = %d v[p].y = %e", p, v[p].y);
                        if (v[p].z != v[p].z) printf("\n DEBUG 3 p = %d v[p].z = %e", p, v[p].z);

                        //chk = ff_model_check_smooth_dv(p);
                        //if ( chk == 0) goto t_end;

                        r0.x += M0p[p] * r[p].x;
                        r0.y += M0p[p] * r[p].y;
                        r0.z += M0p[p] * r[p].z;
                        Mtot += M0p[p];

						if (is_periodic)
						{
							if (r[p].x > 0) per_x_plus = r[p].x - Lx;
							else per_x_plus = r[p].x;
							if (r[p].x <= 0) per_x_minus = r[p].x + Lx;
							else per_x_minus = r[p].x;

							if (r[p].y > 0) per_y_plus = r[p].y - Ly;
							else per_y_plus = r[p].y;
							if (r[p].y <= 0) per_y_minus = r[p].y + Ly;
							else per_y_minus = r[p].y;

							if (r[p].z > 0) per_z_plus = r[p].z - Lz;
							else per_z_plus = r[p].z;
							if (r[p].z <= 0) per_z_minus = r[p].z + Lz;
							else per_z_minus = r[p].z;

							r_per[p][0][0][0].x = per_x_minus;
							r_per[p][0][0][0].y = per_y_minus;
							r_per[p][0][0][0].z = per_z_minus;
							r0_per[0][0][0].x += M0p[p] * per_x_minus;
							r0_per[0][0][0].y += M0p[p] * per_y_minus;
							r0_per[0][0][0].z += M0p[p] * per_z_minus;

							r_per[p][0][0][1].x = per_x_minus;
							r_per[p][0][0][1].y = per_y_minus;
							r_per[p][0][0][1].z = r[p].z;
							r0_per[0][0][1].x += M0p[p] * per_x_minus;
							r0_per[0][0][1].y += M0p[p] * per_y_minus;
							r0_per[0][0][1].z += M0p[p] * r[p].z;

							r_per[p][0][0][2].x = per_x_minus;
							r_per[p][0][0][2].y = per_y_minus;
							r_per[p][0][0][2].z = per_z_plus;
							r0_per[0][0][2].x += M0p[p] * per_x_minus;
							r0_per[0][0][2].y += M0p[p] * per_y_minus;
							r0_per[0][0][2].z += M0p[p] * per_z_plus;

							r_per[p][0][1][0].x = per_x_minus;
							r_per[p][0][1][0].y = r[p].y;
							r_per[p][0][1][0].z = per_z_minus;
							r0_per[0][1][0].x += M0p[p] * per_x_minus;
							r0_per[0][1][0].y += M0p[p] * r[p].y;
							r0_per[0][1][0].z += M0p[p] * per_z_minus;

							r_per[p][0][1][1].x = per_x_minus;
							r_per[p][0][1][1].y = r[p].y;
							r_per[p][0][1][1].z = r[p].z;							
							r0_per[0][1][1].x += M0p[p] * per_x_minus;
							r0_per[0][1][1].y += M0p[p] * r[p].y;
							r0_per[0][1][1].z += M0p[p] * r[p].z;

							r_per[p][0][1][2].x = per_x_minus;
							r_per[p][0][1][2].y = r[p].y;
							r_per[p][0][1][2].z = per_z_plus;							
							r0_per[0][1][2].x += M0p[p] * per_x_minus;
							r0_per[0][1][2].y += M0p[p] * r[p].y;
							r0_per[0][1][2].z += M0p[p] * per_z_plus;

							r_per[p][0][2][0].x = per_x_minus;
							r_per[p][0][2][0].y = per_y_plus;
							r_per[p][0][2][0].z = per_z_minus;
							r0_per[0][2][0].x += M0p[p] * per_x_minus;
							r0_per[0][2][0].y += M0p[p] * per_y_plus;
							r0_per[0][2][0].z += M0p[p] * per_z_minus;

							r_per[p][0][2][1].x = per_x_minus;
							r_per[p][0][2][1].y = per_y_plus;
							r_per[p][0][2][1].z = r[p].z;
							r0_per[0][2][1].x += M0p[p] * per_x_minus;
							r0_per[0][2][1].y += M0p[p] * per_y_plus;
							r0_per[0][2][1].z += M0p[p] * r[p].z;

							r_per[p][0][2][2].x = per_x_minus;
							r_per[p][0][2][2].y = per_y_plus;
							r_per[p][0][2][2].z = per_z_plus;
							r0_per[0][2][2].x += M0p[p] * per_x_minus;
							r0_per[0][2][2].y += M0p[p] * per_y_plus;
							r0_per[0][2][2].z += M0p[p] * per_z_plus;

							r_per[p][1][0][0].x = r[p].x;
							r_per[p][1][0][0].y = per_y_minus;
							r_per[p][1][0][0].z = per_z_minus;
							r0_per[1][0][0].x += M0p[p] * r[p].x;
							r0_per[1][0][0].y += M0p[p] * per_y_minus;
							r0_per[1][0][0].z += M0p[p] * per_z_minus;

							r_per[p][1][0][1].x = r[p].x;
							r_per[p][1][0][1].y = per_y_minus;
							r_per[p][1][0][1].z = r[p].z;
							r0_per[1][0][1].x += M0p[p] * r[p].x;
							r0_per[1][0][1].y += M0p[p] * per_y_minus;
							r0_per[1][0][1].z += M0p[p] * r[p].z;

							r_per[p][1][0][2].x = r[p].x;
							r_per[p][1][0][2].y = per_y_minus;
							r_per[p][1][0][2].z = per_z_plus;
							r0_per[1][0][2].x += M0p[p] * r[p].x;
							r0_per[1][0][2].y += M0p[p] * per_y_minus;
							r0_per[1][0][2].z += M0p[p] * per_z_plus;

							r_per[p][1][1][0].x = r[p].x;
							r_per[p][1][1][0].y = r[p].y;
							r_per[p][1][1][0].z = per_z_minus;
							r0_per[1][1][0].x += M0p[p] * r[p].x;
							r0_per[1][1][0].y += M0p[p] * r[p].y;
							r0_per[1][1][0].z += M0p[p] * per_z_minus;

							r_per[p][1][1][1].x = r[p].x;
							r_per[p][1][1][1].y = r[p].y;
							r_per[p][1][1][1].z = r[p].z;
							r0_per[1][1][1].x += M0p[p] * r[p].x;
							r0_per[1][1][1].y += M0p[p] * r[p].y;
							r0_per[1][1][1].z += M0p[p] * r[p].z;

							r_per[p][1][1][2].x = r[p].x;
							r_per[p][1][1][2].y = r[p].y;
							r_per[p][1][1][2].z = per_z_plus;
							r0_per[1][1][2].x += M0p[p] * r[p].x;
							r0_per[1][1][2].y += M0p[p] * r[p].y;
							r0_per[1][1][2].z += M0p[p] * per_z_plus;

							r_per[p][1][2][0].x = r[p].x;
							r_per[p][1][2][0].y = per_y_plus;
							r_per[p][1][2][0].z = per_z_minus;
							r0_per[1][2][0].x += M0p[p] * r[p].x;
							r0_per[1][2][0].y += M0p[p] * per_y_plus;
							r0_per[1][2][0].z += M0p[p] * per_z_minus;

							r_per[p][1][2][1].x = r[p].x;
							r_per[p][1][2][1].y = per_y_plus;
							r_per[p][1][2][1].z = r[p].z;
							r0_per[1][2][1].x += M0p[p] * r[p].x;
							r0_per[1][2][1].y += M0p[p] * per_y_plus;
							r0_per[1][2][1].z += M0p[p] * r[p].z;

							r_per[p][1][2][2].x = r[p].x;
							r_per[p][1][2][2].y = per_y_plus;
							r_per[p][1][2][2].z = per_z_plus;
							r0_per[1][2][2].x += M0p[p] * r[p].x;
							r0_per[1][2][2].y += M0p[p] * per_y_plus;
							r0_per[1][2][2].z += M0p[p] * per_z_plus;

							r_per[p][2][0][0].x = per_x_plus;
							r_per[p][2][0][0].y = per_y_minus;
							r_per[p][2][0][0].z = per_z_minus;
							r0_per[2][0][0].x += M0p[p] * per_x_plus;
							r0_per[2][0][0].y += M0p[p] * per_y_minus;
							r0_per[2][0][0].z += M0p[p] * per_z_minus;

							r_per[p][2][0][1].x = per_x_plus;
							r_per[p][2][0][1].y = per_y_minus;
							r_per[p][2][0][1].z = r[p].z;
							r0_per[2][0][1].x += M0p[p] * per_x_plus;
							r0_per[2][0][1].y += M0p[p] * per_y_minus;
							r0_per[2][0][1].z += M0p[p] * r[p].z;

							r_per[p][2][0][2].x = per_x_plus;
							r_per[p][2][0][2].y = per_y_minus;
							r_per[p][2][0][2].z = per_z_plus;
							r0_per[2][0][2].x += M0p[p] * per_x_plus;
							r0_per[2][0][2].y += M0p[p] * per_y_minus;
							r0_per[2][0][2].z += M0p[p] * per_z_plus;

							r_per[p][2][1][0].x = per_x_plus;
							r_per[p][2][1][0].y = r[p].y;
							r_per[p][2][1][0].z = per_z_minus;
							r0_per[2][1][0].x += M0p[p] * per_x_plus;
							r0_per[2][1][0].y += M0p[p] * r[p].y;
							r0_per[2][1][0].z += M0p[p] * per_z_minus;

							r_per[p][2][1][1].x = per_x_plus;
							r_per[p][2][1][1].y = r[p].y;
							r_per[p][2][1][1].z = r[p].z;
							r0_per[2][1][1].x += M0p[p] * per_x_plus;
							r0_per[2][1][1].y += M0p[p] * r[p].y;
							r0_per[2][1][1].z += M0p[p] * r[p].z;

							r_per[p][2][1][2].x = per_x_plus;
							r_per[p][2][1][2].y = r[p].y;
							r_per[p][2][1][2].z = per_z_plus;
							r0_per[2][1][2].x += M0p[p] * per_x_plus;
							r0_per[2][1][2].y += M0p[p] * r[p].y;
							r0_per[2][1][2].z += M0p[p] * per_z_plus;

							r_per[p][2][2][0].x = per_x_plus;
							r_per[p][2][2][0].y = per_y_plus;
							r_per[p][2][2][0].z = per_z_minus;
							r0_per[2][2][0].x += M0p[p] * per_x_plus;
							r0_per[2][2][0].y += M0p[p] * per_y_plus;
							r0_per[2][2][0].z += M0p[p] * per_z_minus;

							r_per[p][2][2][1].x = per_x_plus;
							r_per[p][2][2][1].y = per_y_plus;
							r_per[p][2][2][1].z = r[p].z;
							r0_per[2][2][1].x += M0p[p] * per_x_plus;
							r0_per[2][2][1].y += M0p[p] * per_y_plus;
							r0_per[2][2][1].z += M0p[p] * r[p].z;

							r_per[p][2][2][2].x = per_x_plus;
							r_per[p][2][2][2].y = per_y_plus;
							r_per[p][2][2][2].z = per_z_plus;
							r0_per[2][2][2].x += M0p[p] * per_x_plus;
							r0_per[2][2][2].y += M0p[p] * per_y_plus;
							r0_per[2][2][2].z += M0p[p] * per_z_plus;
						}

                        mz_tot += m[p].z;
                        m_tot.x += m[p].x;
                        m_tot.y += m[p].y;
                        m_tot.z += m[p].z;
                        mz_tot_n++;

                        //if (p == 1) ff_model_brownian_validation(p);

                        /*dW[p] = MUL(F[p], drt[p]) + MUL(tau[p], dphi[p]) 
                        - C2[p] * ((v[p].x - dvt[p].x / 2) * drt[p].x + (v[p].y - dvt[p].y / 2) * drt[p].y + (v[p].z - dvt[p].z / 2) * drt[p].z)
                        - gamma_rot[p] * ((w[p].x - dw[p].x / 2) * dphi[p].x + (w[p].y - dw[p].y / 2) * dphi[p].y + (w[p].z - dw[p].z / 2) * dphi[p].z);*/
                    } // end of loop for dv

                    //printf("\n DEBUG 5 r.x = %e v.x = %e", r[50].x, v[50].x);

                    // TODO: need the number of existing (exist_p) particles
                    r0.x /= Mtot;
                    r0.y /= Mtot;
                    r0.z /= Mtot;

					if (is_periodic) for (int t1 = 0; t1 <= 2; t1++) for (int t2 = 0; t2 <= 2; t2++) for (int t3 = 0; t3 <= 2; t3++)
					{
						r0_per[t1][t2][t3].x /= Mtot;
						r0_per[t1][t2][t3].y /= Mtot;
						r0_per[t1][t2][t3].z /= Mtot;
					}

                    I = 0;
					
                    for (p = 1; p <= pN; p++)
					{
                        I += M0p[p] * (pow(r[p].x - r0.x, 2)
                        + pow(r[p].y - r0.y, 2)
                        + pow(r[p].z - r0.z, 2));

						if (is_periodic) for (int t1 = 0; t1 <= 2; t1++) for (int t2 = 0; t2 <= 2; t2++) for (int t3 = 0; t3 <= 2; t3++)
						{
							I_per[t1][t2][t3] += M0p[p] * (
														   pow(r_per[p][t1][t2][t3].x - r0_per[t1][t2][t3].x, 2)
														 + pow(r_per[p][t1][t2][t3].y - r0_per[t1][t2][t3].y, 2)
														 + pow(r_per[p][t1][t2][t3].z - r0_per[t1][t2][t3].z, 2));
						}
					}

					if (is_periodic)
					{
						I_min_per = I_per[0][0][0];
						for (int t1 = 0; t1 <= 2; t1++) for (int t2 = 0; t2 <= 2; t2++) for (int t3 = 0; t3 <= 2; t3++)
							if (I_per[t1][t2][t3] < I_min_per) I_min_per = I_per[t1][t2][t3];
					}
					//printf("\n DEBUG I000 = %e I = %e", I_per[1][1][1], I);


                    mz_tot *= pN / mz_tot_n; // in loop is interrupted then needs to increase
                    mz_glob += mz_tot;

                    m_tot_glob.x += m_tot.x;
                    m_tot_glob.y += m_tot.y;
                    m_tot_glob.z += m_tot.z;

                    t += dt;
                    dt_red -= dt;

                    /*
                    for (p = 1; p <= pN; p++) if (exist_p[p])
                    {
                    Ek_tr  += M0p[p] * MUL(v[p],v[p]) / 2.0;
                    Ek_rot += I0p[p] * MUL(w[p],w[p]) / 2.0;
                    };

                    Ek = Ek_tr + Ek_rot;
                    ff_model_update_dT();*/

                    for (p = 1; p <= pN; p++)
                    {
                        Rp_to_c[p] = sqrt(MUL(r[p], r[p]));
                        ff_model_update_conc_in_oleic(p);

                        //if (dt_red <= 0)
                        //{
                        //ff_model_inst_random_trans_update(p);

                        /*dv_r[p] = sqrt(3 * kb * dT / M0p[p]);
                        dw_r[p] = sqrt(3 * kb * dT / I0p[p]);

                        theta_0_r = (*var_uni)() * pi;   // rotation vector random direction
                        phi_0_r = (*var_uni)() * 2 * pi;
                        v[p].x += dv_r[p] * sin(theta_0_r) * cos(phi_0_r);
                        v[p].y += dv_r[p] * sin(theta_0_r) * sin(phi_0_r);
                        v[p].z += dv_r[p] * cos(theta_0_r);

                        theta_0_r = (*var_uni)() * pi;  
                        phi_0_r = (*var_uni)() * 2 * pi;
                        w[p].x += dw_r[p] * sin(theta_0_r) * cos(phi_0_r);
                        w[p].y += dw_r[p] * sin(theta_0_r) * sin(phi_0_r);
                        w[p].z += dw_r[p] * cos(theta_0_r);*/
                        //}
                    }

                    if (phi_vol_fract_oleic_0 == 0) phi_vol_fract_oleic_0 = phi_vol_fract_oleic;
                    if ((is_oleic) && (phi_vol_fract_oleic_0 > 0) && (phi_vol_fract_oleic == phi_vol_fract_oleic))
                        R_oleic = R_oleic_0 * (phi_vol_fract_oleic / phi_vol_fract_oleic_0);
                    //printf("\n TRACE-1 %e %e", phi_vol_fract_oleic_0, phi_vol_fract_oleic);
                    V_oleic = (4 / 3.0) * pi * pow(R_oleic, 3);
                    //phi_vol_fract_oleic /= V_oleic;
                    phi_vol_fract_oleic /= (Lx * Ly * Lz);
                    phi_vol_fract_oleic *= 100;
                    
                    P_pgas = kb * T * pN_oleic_drop / (V_oleic - V0_tot_oleic); // Van 't Hoff equation [Landau-V, eq. (88,3)]
                    P_sf_oleic = 4 * sigma_sf / R_oleic;
                    
                    if (isPGasMode)
                    {
                        if (P_pgas / P_sf_oleic > 1.0) isPGasPrevail = 1.0;
                        else isPGasPrevail = 0.0;
                    }
                    //printf("\n %e", P_pgas / P_sf_oleic);

                    //if (dt_red <= 0) dt_red = dt0;

                    if (auto_reversal) ff_model_auto_hyst();
                    if (step % 100000 == 0) ff_io_autosave();
                    else if ((step < 100000) && (step % 10000 == 0)) ff_io_autosave();
                    else if ((step < 10000) && (step % 1000 == 0)) ff_io_autosave();

                    if (setting_plot)
                    {
                        if (is_periodic) I = I_min_per;
						I_glob = I;
						if (step % 1000 == 0) ff_io_save_setting(m_tot,I);
                        else if ((step < 1000) && (step % 100 == 0)) ff_io_save_setting(m_tot,I);
                        else if ((step < 100) && (step % 10 == 0)) ff_io_save_setting(m_tot,I);
                    }

                    // end for case of interrupted kinetic loops
t_end:

                    if (slow_steps > 0) slow_steps--;
                    if (slow_steps%10 == 1) dt *= 2;

                    //for (p = 1; p <= pN; p++) ff_model_inst_random_trans_update(p);

    } // time_go

    ff_mgr_show_next_step();
}

int ff_model_check_walls(long p)
{
    int res = 0;
    double dr;

    if (!is_periodic)
    {
        //walls
        // Oz
        if (r[p].z < -Lz / 2.0)
        {
            r[p].z = -Lz / 2.0;
            v[p].z *= -1;
            res = 1;
        }

        if (r[p].z >  Lz / 2.0)
        {
            r[p].z =  Lz / 2.0;
            v[p].z *= -1;
            res = 1;
        }

        // Ox
        if (r[p].x < -Lx / 2.0)
        {
            r[p].x = -Lx / 2.0;
            v[p].x *= -1;
            res = 1;
        }

        if (r[p].x >  Lx / 2.0)
        {
            r[p].x =  Lx / 2.0;
            v[p].x *= -1;
            res = 1;
        }

        // Oy
        if (r[p].y < -Ly / 2.0)
        {
            r[p].y = -Ly / 2.0;
            v[p].y *= -1;
            res = 1;
        }

        if (r[p].y >  Ly / 2.0)
        {
            r[p].y =  Ly / 2.0;
            v[p].y *= -1;
            res = 1;
        }
    }
    else
    {
        //walls
        // Oz
        if (r[p].z < -Lz / 2.0)
        {
            r[p].z += Lz;
            res = 1;
        }

        if (r[p].z >  Lz / 2.0)
        {
            r[p].z -=  Lz;
            res = 1;
        }

        // Ox
        if (r[p].x < -Lx / 2.0)
        {
            r[p].x += Lx;
            res = 1;
        }

        if (r[p].x >  Lx / 2.0)
        {
            r[p].x -= Lx;
            res = 1;
        }

        // Oy
        if (r[p].y < -Ly / 2.0)
        {
            r[p].y += Ly;
            res = 1;
        }

        if (r[p].y >  Ly / 2.0)
        {
            r[p].y -= Ly;
            res = 1;
        }
    }

    return res;
}

/*void ff_model_check_overlapp(long p)
{
long ps;
double dR, dx, dy, dz;
double dr;

//for (p = 1; p <= pN; p++)
for(ps = 1; ps <= pN; ps ++)
if (p != ps)
{
dx = r[ps].x - r[p].x - drt[p].x;
dy = r[ps].y - r[p].y - drt[p].y;
dz = r[ps].z - r[p].z - drt[p].z;

dR = sqrt(dx * dx + dy * dy + dz * dz);

if (dR <= Rp[p] + Rp[ps])
{
//r[ps].x -= drt[ps].x;
//r[ps].y -= drt[ps].y;
//r[ps].z -= drt[ps].z;

//v[ps].x -= dvt[ps].x;
//v[ps].y -= dvt[ps].y;
//v[ps].z -= dvt[ps].z;
dr = Rp[p] + Rp[ps] - dR;
drt[p].x += - (dx / dR) * dr - delta * 0.1;
drt[p].y += - (dy / dR) * dr - delta * 0.1;
drt[p].z += - (dz / dR) * dr - delta * 0.1;
}
}
}*/

// External heterogeneous field. Source is shifted on X axis. It is normilized on ext. magnetic moment directed along X axis.
ff_vect_t B_het(double x1, double y1, double z1, int what_shift)
{
    double x0, y0, z0;
    double x, y, z;
    double r, mr;
    ff_vect_t B;

    x0 = - Lx * (what_shift == 1);
    y0 = - Ly * (what_shift == 2);
    z0 = - Lz * (what_shift == 3);

    x = x1 - x0;
    y = y1 - y0;
    z = z1 - z0;

    if (what_shift == 1) mr = x;
    if (what_shift == 2) mr = y;
    if (what_shift == 3) mr = z;

    r = sqrt(x * x + y * y + z * z);

    B.x = C5 * (3 * x * mr / pow(r, 5) - (what_shift == 1) / pow(r, 3));
    B.y = C5 * (3 * y * mr / pow(r, 5) - (what_shift == 2) / pow(r, 3));
    B.z = C5 * (3 * z * mr / pow(r, 5) - (what_shift == 3) / pow(r, 3));

    return B;
}

// External inhomogeneous field. The sample is located in the center of the big electromagnet consisting of 2 concentrators.
// Symmetry: [Infinite-fold symmetry axis inf_z], [Reflection plane (xy)].
ff_vect_t B_het_bem(ff_vect_t r)
{
    double r0;
    ff_vect_t B;

    r0 = sqrt(MUL(r,r));

    B.x = B.y = 0;
    B.z = B0 * (1 - gradPerc * r0 / gradL);

    return B;
}

// External inhomogeneous field. The sample is located in the center of the big electromagnet consisting of 2 concentrators.
// Symmetry: [Infinite-fold symmetry axis inf_z], [Reflection plane (xy)].
ff_vect_t dBz_het_bem(ff_vect_t r)
{
    double r0;
    ff_vect_t dBz;

    r0 = sqrt(MUL(r,r));

    if (r0 != 0)
    {
        dBz.x = - B0 * gradPerc * r.x / (gradL * r0); // it means dBz.x == dBz / dx
        dBz.y = - B0 * gradPerc * r.y / (gradL * r0);
        dBz.z = - B0 * gradPerc * r.z / (gradL * r0);
    }
    else
    {
        dBz.x = dBz.y = dBz.z = 0;
    }

    return dBz;
}

ff_vect_t Bext(double x, double y, double z)
{
    ff_vect_t tBext,r;
    double normX, normY, normZ;

    r.x = x; r.y = y; r.z = z;

    if (!manual_field_control)
    {
        tBext.x = tBext.y = 1E-50;
        tBext.z = B0 * sin(kB * pi / 2.0);
    }
    else
    {
        if (ext_field_is_homo)
        {
            tBext.x = BmanX /*Oe*/ * (1.0 / (4.0 * pi * 1E-3)) * mu0; // Tesla
            tBext.y = BmanY /*Oe*/ * (1.0 / (4.0 * pi * 1E-3)) * mu0; // Tesla
            tBext.z = BmanZ /*Oe*/ * (1.0 / (4.0 * pi * 1E-3)) * mu0; // Tesla
        }
        else
        {
            /*normX = (BmanX * (1.0 / (4.0 * pi * 1E-3)) * mu0 / B_het(0, 0, 0, 1).x);
            normY = (BmanY * (1.0 / (4.0 * pi * 1E-3)) * mu0 / B_het(0, 0, 0, 2).y);
            normZ = (BmanZ * (1.0 / (4.0 * pi * 1E-3)) * mu0 / B_het(0, 0, 0, 3).z);

            tBext.x = normX * B_het(x, y, z, 1).x + normY * B_het(x, y, z, 2).x + normZ * B_het(x, y, z, 3).x;
            tBext.y = normX * B_het(x, y, z, 1).y + normY * B_het(x, y, z, 2).y + normZ * B_het(x, y, z, 3).y;
            tBext.z = normX * B_het(x, y, z, 1).z + normY * B_het(x, y, z, 2).z + normZ * B_het(x, y, z, 3).z;*/
            tBext = B_het_bem(r);
        }
    }

    return tBext;
}

void ff_model_init(void)
{
    long p, ps, tp, p_prev;
    double theta, phi;
    ff_vect_t dr;
    double dR;
    double sigma;
    long i,j;
    double t_temp_1 = 0;
    double dr_tmp, dr_min;
    int cont_flag;
    long i_attempt;

    double start_scale = 0.99;

    dt = dt0;

    R_oleic = R_oleic_0;

    if (is_uniform_field_test) Lz *= 5;
    
    t_temp = T + ta0;
    if (t_temp > 90) t_temp_1 = 90;
    if (t_temp < 20) t_temp_1 = 20;
    sigma_sf = sigma_sf_nano * (a_sigma_sf + b_sigma_sf * t_temp_1);
    eta_oleic = a3_eta_oleic * pow(t_temp_1, 3) + a2_eta_oleic * pow(t_temp_1, 2) + a1_eta_oleic * pow(t_temp_1, 1) + a0_eta_oleic;

    N_o = N_oa * k_o / 0.5;

    //printf("\n tau0 = %e", Ms / (2 * alpha_damp * gamma_e * K1));
    
    // Brownian motion -  parameters
    ///////////////////////////////////////////////////

    k_force_adapt = 1;
    //k_force_adapt = k_force_adapt_0 / sqrt(1E-9); //sqrt(1E-9) is selected regular dt

    // Dimensionless variance (sigma^2) of the random displacement along the single axis e_x
    sigma = 1;

    srand( (unsigned)time( NULL ) );
    rng.seed(static_cast<unsigned int>(std::time(0)));

    //boost::normal_distribution<> nd(0.0, 1.0);
    nd = new boost::normal_distribution<> (0.0, sigma);
    ud = new boost::uniform_01<> ();

    //boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, *nd);
    var_nor = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > (rng, *nd);
    var_uni = new boost::variate_generator<boost::mt19937&, boost::uniform_01<> > (rng, *ud);

    //ff_model_dir_init();

    t = 0; // time

    //printf("\n m0 = %e, N = %d", m0, pN);
    //printf("\n %e", M0);
    glob_start_t = start_t;

    m_tot_glob.x = m_tot_glob.y = m_tot_glob.z = 0;

    for (long h = 1; h <= 20; h++)
    {
        B_hyst[h] = B_hyst_n[h] = Mz_hyst[h] = Mz_hyst_n[h] = 0;
    }

    ff_model_size_dispersion_init();
    p = 1;
    while (p <= pN)
    {
again:

        if (start_ideal)
        {
            r[p].x = start_scale * (Rp[p] - 0.49 * Lx + (2 * 0.49 * Lx - Rp[p]) * (*var_uni)());
            r[p].y = start_scale * (Rp[p] - 0.49 * Ly + (2 * 0.49 * Ly - Rp[p]) * (*var_uni)());
            r[p].z = start_scale * (Rp[p] - 0.49 * Lz + (2 * 0.49 * Lz - Rp[p]) * (*var_uni)());
            Rp_to_c[p] = sqrt(MUL(r[p], r[p]));

            /*r[p].x = - R_oleic + 2 * R_oleic * (*var_uni)();
            r[p].y = - R_oleic + 2 * R_oleic * (*var_uni)();
            r[p].z = - R_oleic + 2 * R_oleic * (*var_uni)();*/

            w[p].x = w[p].y = w[p].z = 0;

            //if (p == 5)
            //r[p].x = r[p].y = r[p].z = 0;

            /*if (!is_uniform_field_test)
            {
                Rp_to_c[p] = sqrt(MUL(r[p], r[p]));
                if (Rp_to_c[p] > R_oleic - Rp[p]) goto again;
            }
			*/
            for (tp = 1; tp < p; tp++)
            {
                dr.x = r[p].x - r[tp].x;
                dr.y = r[p].y - r[tp].y;
                dr.z = r[p].z - r[tp].z;

                dR = sqrt(MUL(dr,dr));

                if (dR <= Rp0[p] + Rp0[tp] + delta_r_init) goto again;
            }
        } // start_ideal
        else
        {
            p_prev = 1;
            if (p == 1) {r[p].x = r[p].y = r[p].z = 0;}
            else
            {
cont3:                
                for (i_attempt = 1; i_attempt <= 100; i_attempt++)
                {
                    dr_min = Rp0[p_prev] + Rp0[p] + delta_r_init;
                    
                    theta = pi * (*var_uni)();
                    phi = 2 * pi * (*var_uni)();
                    
                    r[p].x = r[p_prev].x + dr_min * sin(theta) * cos(phi);
                    r[p].y = r[p_prev].y + dr_min * sin(theta) * sin(phi);
                    r[p].z = r[p_prev].z + dr_min * cos(theta);

                    cont_flag = 0;
                    for (tp = 1; tp < p; tp++)
                    {
                        dr.x = r[p].x - r[tp].x;
                        dr.y = r[p].y - r[tp].y;
                        dr.z = r[p].z - r[tp].z;

                        dR = sqrt(MUL(dr,dr));

                        if (dR < Rp0[p] + Rp0[tp] + delta_r_init) {cont_flag = 1; goto cont1;}
                    }
cont1:             if (!cont_flag) goto cont2;
                }
cont2:
                if (cont_flag) {p_prev ++; goto cont3;}
            }
        }

        v[p].x = v[p].y = v[p].z = 0;
        //v_r[p].x = v_r[p].y = v_r[p].z = 0;
        w[p].x = w[p].y = w[p].z = 0;
        //w_r[p].x = w_r[p].y = w_r[p].z = 0;
        drt_r[p].x = drt_r[p].y = drt_r[p].z = 0;
        dphi_r[p].x = dphi_r[p].y = dphi_r[p].z = 0;
        //dW[p] = 0;
        //for (tp = 1; tp <= pN; tp++) aggregated_p[p][tp] = 0;

        /*theta = pi * rand() / 32768.0;
        phi = 2 * pi * rand() / 32768.0;

        m[p].x = m0p[p] * sin(theta) * cos(phi);
        m[p].y = m0p[p] * sin(theta) * sin(phi);
        m[p].z = m0p[p] * cos(theta);*/

        exist_p[p] = 1;

        //m_freeze[p] = 0;
        m_sat[p] = 0;

        k_force_adapt_p_x[p] = k_force_adapt_p_y[p] = k_force_adapt_p_z[p] = k_force_adapt_p_rot_x[p] = k_force_adapt_p_rot_y[p] = k_force_adapt_p_rot_z[p] = 1.0;
        dT_p_x[p] = dT_p_y[p] = dT_p_z[p] = dT_p_rot_x[p] = dT_p_rot_y[p] = dT_p_rot_z[p] = T;
        dT_prev_p_x[p] = dT_prev_p_y[p] = dT_prev_p_z[p] = dT_prev_p_rot_x[p] = dT_prev_p_rot_y[p] = dT_prev_p_rot_z[p] = T;
        T_basic_p_x[p] = T_basic_p_y[p] = T_basic_p_z[p] = T_basic_p_rot_x[p] = T_basic_p_rot_y[p] = T_basic_p_rot_z[p] = 0;
        T_mean_p_x[p] = T_mean_p_y[p] = T_mean_p_z[p] = T_mean_p_rot_x[p] = T_mean_p_rot_y[p] = T_mean_p_rot_z[p] = 0;
        T_mean_loc_p_x[p] = T_mean_loc_p_y[p] = T_mean_loc_p_z[p] = T_mean_loc_p_rot_x[p] = T_mean_loc_p_rot_y[p] = T_mean_loc_p_rot_z[p] = 0;
        //T_mean_loc_prev_p[p] = 0;
        //T_mean_loc_prev_revert_p[p] = 0;
        k_mean_p[p] = 0;
        k_delta_force_rel[p] = 0;

        is_inside_oleic[p] = 1;
        is_temp_sat[p] = 0;
        k_force_adapt_p_0[p] = k_force_adapt_0;

		G_dd[p] = 0;

        p++;
    }

    if (load_at_start) ff_io_load(0);

    //if (start_sediment) ff_model_init_sediment();

    r_brown_valid_0 = r[1];

    //R_oleic *= 8;
}

/*void ff_model_init_sediment(void)
{
long i, j, k;
ff_vect_t r0, r1, rb; //rb - basic; we build sphere around it
long n;
double a = (2 + 0.1) * R0;
long p, pb;
int res;
long pt;
int jdir;

rb.x = rb.y = 0;
//rb.z = -0.99 * Lz / 2.0;
rb.z = 0;

jdir = 0;

p = 1;
pb = p;

r[p].x = rb.x;
r[p].y = rb.y;
r[p].z = rb.z;

p++;

while (p <= pN)
{
jdir++;

r1.x = rb.x + a * dir110[jdir].x / sqrt(2.0);
r1.y = rb.y + a * dir110[jdir].y / sqrt(2.0);
r1.z = rb.z + a * dir110[jdir].z / sqrt(2.0);

if ((fabs(r1.x) > Lx / 2.0)||(fabs(r1.y) > Ly / 2.0)||(fabs(r1.z) > Lz / 2.0))
{res = 0; goto m1;}

res = 1;
for(pt = 1; pt < p; pt++)
{
if ((fabs(r1.x - r[pt].x) > 0.1 * R0)||(fabs(r1.y - r[pt].y) > 0.1 * R0)||(fabs(r1.z - r[pt].z) > 0.1 * R0)) res = 1;
else {res = 0; goto m1;}
}
m1:

if (res == 1)
{
r[p].x = r1.x;
r[p].y = r1.y;
r[p].z = r1.z;

p++;
}
if ( jdir == 12)
{
jdir = 0;

pb++;

rb.x = r[pb].x;
rb.y = r[pb].y;
rb.z = r[pb].z;
}
} // p for
}

void ff_model_dir_init(void)
{

long i;

i = 1;
dir110[i].x = 1; dir110[i].y = 1; dir110[i].z = 0;
i++;
dir110[i].x = 1; dir110[i].y = 0; dir110[i].z = 1;
i++;
dir110[i].x = 0; dir110[i].y = 1; dir110[i].z = 1;
i++;

dir110[i].x = -1; dir110[i].y =  1; dir110[i].z = 0;
i++;
dir110[i].x = -1; dir110[i].y =  0; dir110[i].z = 1;
i++;
dir110[i].x =  0; dir110[i].y = -1; dir110[i].z = 1;
i++;


dir110[i].x = -1; dir110[i].y = -1; dir110[i].z = 0;
i++;
dir110[i].x = -1; dir110[i].y = 0; dir110[i].z = -1;
i++;
dir110[i].x = 0; dir110[i].y = -1; dir110[i].z = -1;
i++;

dir110[i].x = 1; dir110[i].y =  -1; dir110[i].z = 0;
i++;
dir110[i].x = 1; dir110[i].y =  0; dir110[i].z = -1;
i++;
dir110[i].x =  0; dir110[i].y = 1; dir110[i].z = -1;

}*/

// Update of the random motion
void ff_model_effective_random_motion_update(long p)
{
    double Px, Py, Pz, tau_r_phi; // instantiated effective random force
    double theta_0, phi_0; // random direction of the torque 
    double dx, dy, dz, dphi; // instantiated displacements for time dt * k_bm_inst_max
    //double dt0;
    double D, D_rot, gamma, gamma_rot;
    double gamma_oleic = 0, gamma_car = 0; // oleic vs. carrier liquid
    double M0, I0;

    double speed_ballance = 1;
    double dR = 0;

    ff_vect_t e1, e2, e3; // dynamic coordinate system
    double P1, P2, P3;
    double tmp = 0;

    double fract_oleic_positive = 0, fract_oleic_negative = 0, fract_car_positive = 0, fract_car_negative = 0;
    double N_mol_col = 0;

    double sigma = 0, sigma_rot = 0;

    //dt0 = dt * k_bm_inst_max;
    //dt0 = dt;

    gamma_oleic = 6 * pi * eta_oleic * Rp[p];
    gamma_car = 6 * pi * eta_car * Rp[p];
    if ((Rp_to_c[p] > R_oleic) || (!is_oleic))
    {
        gamma = gamma_car;
        gamma_rot = 8 * pi * eta_car * pow(Rp[p], 3);
    }
    else
    {
        gamma = gamma_oleic;
        gamma_rot = 8 * pi * eta_oleic * pow(Rp[p], 3);
    }

    D = kb * T / gamma;
    D_rot = kb * T / gamma_rot;
    M0 = M0p[p];
    I0 = I0p[p];

    // instantiation of position change
    // [Langevin equation + Stokes' law]
    sigma = D * (2 * dt + (M0 / gamma) * (- 3 + 4 * exp(- gamma * dt / M0) - exp(- 2 * gamma * dt / M0)));
    dx = (*var_nor)() * sqrt(sigma);
    dy = (*var_nor)() * sqrt(sigma);
    dz = (*var_nor)() * sqrt(sigma);

    drt_r[p].x = dx;
    drt_r[p].y = dy;
    drt_r[p].z = dz;

    // instantiation of rotation of the magnetization direction
    // [Euler-Langevin equation + Stokes' law]
    sigma_rot = D_rot * (2 * dt + (I0 / gamma_rot) * (- 3 + 4 * exp(- gamma_rot * dt / I0) - exp(- 2 * gamma_rot * dt / I0)));
    dphi = (*var_nor)() * sqrt(3 * sigma_rot); // rotation magnitude
    theta_0 = (*var_uni)() * pi;   // rotation vector random direction
    phi_0 = (*var_uni)() * 2 * pi;

    dphi_r[p].x = dphi * sin(theta_0) * cos(phi_0);
    dphi_r[p].y = dphi * sin(theta_0) * sin(phi_0);
    dphi_r[p].z = dphi * cos(theta_0);

    /*Px = (gamma * dx) / (dt0 - M0 * (1 - exp(- gamma * dt0 / M0)) / gamma);
    Py = (gamma * dy) / (dt0 - M0 * (1 - exp(- gamma * dt0 / M0)) / gamma);
    Pz = (gamma * dz) / (dt0 - M0 * (1 - exp(- gamma * dt0 / M0)) / gamma);
    tau_r_phi = (gamma_rot * dphi ) / (dt0 - I0 * (1 - exp(- gamma_rot * dt0 / I0)) / gamma_rot);*/

    speed_ballance = 1;
    dR = Rp_to_c[p] + Rp0[p] - R_oleic;
    //if ((dR > 0) && (dR <= 2 * Rp[p]) && (is_oleic)) speed_ballance = 1 + (sqrt(eta / eta_oleic) - 1) * dR / (2 * Rp[p]);*/ // this is correct only for the damping mode. Inertia mode (small dt) should disable this
    
    /*if ((dR > 0) && (dR <= 2 * Rp0[p]) && (is_oleic))
    {
        tmp = sqrt(MUL(r[p], r[p]));
        e3.x = r[p].x / tmp;
        e3.y = r[p].y / tmp;
        e3.z = r[p].z / tmp;

        e1.x = - r[p].y;
        e1.y =   r[p].x;
        e1.z =   0;

        if ((r[p].x == 0) && (r[p].y == 0)) e1.x = 1;
        tmp = sqrt(MUL(e1, e1));
        e1.x /= tmp; e1.y /= tmp; e1.z /= tmp;

        e2.x = e3.y * e1.z - e3.z * e1.y;
        e2.y = e3.z * e1.x - e3.x * e1.z;
        e2.z = e3.x * e1.y - e3.y * e1.x;

        if (dR >= Rp0[p]) 
        {
            fract_oleic_positive = (2 * Rp0[p] - dR) / Rp0[p];
            fract_oleic_negative = 0;
            fract_car_positive = (dR - Rp0[p]) / Rp0[p];
            fract_car_negative = 1;
        }
        else
        {
            fract_oleic_positive = 1;
            fract_oleic_negative = (Rp0[p] - dR) / Rp0[p];
            fract_car_positive = 0;
            fract_car_negative = dR / Rp0[p];
        }

        N_mol_col = (*var_nor)();

        if (N_mol_col >= 0)
            P3 =   N_mol_col * sqrt(2 * kb * T * (fract_oleic_positive * gamma_oleic + fract_car_positive * gamma_car) / dt) *
            (sqrt(gamma_oleic) * fract_oleic_positive + sqrt(gamma_car) * fract_car_positive);
        else
            P3 =   N_mol_col * sqrt(2 * kb * T * (fract_oleic_negative * gamma_oleic + fract_car_negative * gamma_car) / dt) *
            (sqrt(gamma_oleic) * fract_oleic_negative + sqrt(gamma_car) * fract_car_negative);

        // sum of dispersions is a dispersion of sums
        P1 = (*var_nor)() * sqrt(2 * kb * T * (gamma_car * (dR / (2 * Rp0[p])) + gamma_oleic * (1 - dR / (2 * Rp0[p]))) / dt) *
            (sqrt(gamma_car / gamma_oleic) * (dR / (2 * Rp0[p])) + 1 - dR / (2 * Rp0[p])); // speed_ballance (component of the k_force_adapt_p)
        P2 = (*var_nor)() * sqrt(2 * kb * T * (gamma_car * (dR / (2 * Rp0[p])) + gamma_oleic * (1 - dR / (2 * Rp0[p]))) / dt) *
            (sqrt(gamma_car / gamma_oleic) * (dR / (2 * Rp0[p])) + 1 - dR / (2 * Rp0[p]));

        Px = e1.x * P1 + e2.x * P2 + e3.x * P3;
        Py = e1.y * P1 + e2.y * P2 + e3.y * P3;
        Pz = e1.z * P1 + e2.z * P2 + e3.z * P3;
    }
    else*/
    //{
        /*Px = (*var_nor)() * sqrt(2 * kb * T * gamma / dt0);
        Py = (*var_nor)() * sqrt(2 * kb * T * gamma / dt0);
        Pz = (*var_nor)() * sqrt(2 * kb * T * gamma / dt0);
        tau_r_phi = (*var_nor)() * sqrt(6 * kb * T * gamma_rot / dt0);*/
    //}

    //printf("\n %e", sqrt(6 * D * 1) / (2 * Rp0[p]));

    //if ((dR > 0) && (dR >  2 * Rp0[p]) && (is_oleic)) speed_ballance = sqrt(eta_car / eta_oleic); // component of the k_force_adapt_p
    //if (dt < 1.5E-12) speed_ballance = 1;

    /*Px = (*var_nor)() * sqrt(2 * kb * T * gamma / dt); 
    Py = (*var_nor)() * sqrt(2 * kb * T * gamma / dt); 
    Pz = (*var_nor)() * sqrt(2 * kb * T * gamma / dt); */
    //tau_r_phi = (*var_nor)() * sqrt(6 * kb * T * gamma_rot / dt_red);

    //printf("\n t1 = %e", M0 / gamma);
    // printf("\n Viscous limit - saturation of velocity time %e", M0 / gamma);
    // printf("\n Viscous limit - saturation of angular velocity time %e", I0 / gamma_rot);
    //printf("\n %e", gamma * dx / (M0 * v[p].x));

    /*P[p].x = Px * k_force_adapt_p_x[p] * speed_ballance;
    P[p].y = Py * k_force_adapt_p_y[p] * speed_ballance;
    P[p].z = Pz * k_force_adapt_p_z[p] * speed_ballance;

    tau_r[p].x = tau_r_phi * sin(theta_0) * cos(phi_0) * k_force_adapt_p_rot_x[p] * speed_ballance;
    tau_r[p].y = tau_r_phi * sin(theta_0) * sin(phi_0) * k_force_adapt_p_rot_y[p] * speed_ballance;
    tau_r[p].z = tau_r_phi * cos(theta_0) * k_force_adapt_p_rot_z[p] * speed_ballance;*/
}

void ff_model_update_dT(void)
{
    T_basic = (2 / 6.0) * Ek / (kb * pN); // degree of freedom number is 6
    //dT = T - T_basic;
    T_mean += T_basic;
    k_mean ++;
    T_mean_loc += T_basic;

    if (dT < 0)
    {
        //printf("\n WARNING: dT < 0");
        //dT = 0;
    }

    //printf("\n dT = %e", dT);

    if (k_bm_inst == k_bm_inst_max - 1)
    {
        dT_prev = dT;
        dT = T - T_mean_loc / k_bm_inst;
        //if (dT > 0) k_force_adapt *= k_force_adapt_0;
        //else k_force_adapt /= k_force_adapt_0;
        T_mean_loc_prev_revert = T_mean_loc_prev;
        T_mean_loc_prev = T_mean_loc;
    }
    //k_bm_inst ++;

    //if ((T_mean > T) && (step % 10 == 0)) k_force_adapt /= k_force_adapt_0;
    //if ((T_mean < T) && (step % 10 == 0)) k_force_adapt *= k_force_adapt_0;
}

void ff_model_update_dT_p(long p)
{
    double rel_T = 1;
    double T_mean_p_tol = 0;

    T_basic_p_x[p] = 2 * Ekp_x[p] / kb; // degree of freedom number is 6
    T_basic_p_y[p] = 2 * Ekp_y[p] / kb;
    T_basic_p_z[p] = 2 * Ekp_z[p] / kb;
    T_basic_p_rot_x[p] = 2 * Ekp_rot_x[p] / kb;
    T_basic_p_rot_y[p] = 2 * Ekp_rot_y[p] / kb;
    T_basic_p_rot_z[p] = 2 * Ekp_rot_z[p] / kb;
    //dT = T - T_basic;
    T_mean_p_x[p] += T_basic_p_x[p];
    T_mean_p_y[p] += T_basic_p_y[p];
    T_mean_p_z[p] += T_basic_p_z[p];
    T_mean_p_rot_x[p] += T_basic_p_rot_x[p];
    T_mean_p_rot_y[p] += T_basic_p_rot_y[p];
    T_mean_p_rot_z[p] += T_basic_p_rot_z[p];
    k_mean_p[p] ++;
    
    T_mean_p_tol = (T_mean_p_x[p] + T_mean_p_y[p] + T_mean_p_z[p] + T_mean_p_rot_x[p] + T_mean_p_rot_y[p] + T_mean_p_rot_z[p]) / (6.0 * k_mean_p[p]);
    if (T_mean_p_tol >= T) is_temp_sat[p] = 1;

    //if (is_temp_sat[p]) k_force_adapt_p_0[p] = 1 + k_bm_inst_max * (k_force_adapt_0 - 1.0) / k_mean_p[p];

    T_mean_loc_p_x[p] += T_basic_p_x[p];
    T_mean_loc_p_y[p] += T_basic_p_y[p];
    T_mean_loc_p_z[p] += T_basic_p_z[p];
    T_mean_loc_p_rot_x[p] += T_basic_p_rot_x[p];
    T_mean_loc_p_rot_y[p] += T_basic_p_rot_y[p];
    T_mean_loc_p_rot_z[p] += T_basic_p_rot_z[p];

    //printf("\n dT = %e", dT);

    if (k_bm_inst == k_bm_inst_max - 1)
    {
        dT_prev_p_x[p] = dT_p_x[p];
        dT_prev_p_y[p] = dT_p_y[p];
        dT_prev_p_z[p] = dT_p_z[p];
        dT_prev_p_rot_x[p] = dT_p_rot_x[p];
        dT_prev_p_rot_y[p] = dT_p_rot_y[p];
        dT_prev_p_rot_z[p] = dT_p_rot_z[p];

        if (!(is_temp_sat[p]))
        {
            dT_p_x[p] = T - T_mean_loc_p_x[p] / k_bm_inst;
            dT_p_y[p] = T - T_mean_loc_p_y[p] / k_bm_inst;
            dT_p_z[p] = T - T_mean_loc_p_z[p] / k_bm_inst;
            dT_p_rot_x[p] = T - T_mean_loc_p_rot_x[p] / k_bm_inst;
            dT_p_rot_y[p] = T - T_mean_loc_p_rot_y[p] / k_bm_inst;
            dT_p_rot_z[p] = T - T_mean_loc_p_rot_z[p] / k_bm_inst;
        }
        else
        {
            dT_p_x[p] = T - T_mean_p_x[p] / k_mean_p[p];
            dT_p_y[p] = T - T_mean_p_y[p] / k_mean_p[p];
            dT_p_z[p] = T - T_mean_p_z[p] / k_mean_p[p];
            dT_p_rot_x[p] = T - T_mean_p_rot_x[p] / k_mean_p[p];
            dT_p_rot_y[p] = T - T_mean_p_rot_y[p] / k_mean_p[p];
            dT_p_rot_z[p] = T - T_mean_p_rot_z[p] / k_mean_p[p];
        }

        //if (dT_p[p] < - 5 * T) k_force_adapt_p[p] = 1; //rel_T = (T_mean_loc_p[p] / k_bm_inst) / T;
        if (dT_p_x[p] > 0) k_force_adapt_p_x[p] *= k_force_adapt_p_0[p];
        else k_force_adapt_p_x[p] /= k_force_adapt_p_0[p];

        k_force_adapt_mean += k_force_adapt_p_x[p];

        if (dT_p_y[p] > 0) k_force_adapt_p_y[p] *= k_force_adapt_p_0[p];
        else k_force_adapt_p_y[p] /= k_force_adapt_p_0[p];

        k_force_adapt_mean += k_force_adapt_p_y[p];

        if (dT_p_z[p] > 0) k_force_adapt_p_z[p] *= k_force_adapt_p_0[p];
        else k_force_adapt_p_z[p] /= k_force_adapt_p_0[p];

        k_force_adapt_mean += k_force_adapt_p_z[p];

        if (dT_p_rot_x[p] > 0) k_force_adapt_p_rot_x[p] *= k_force_adapt_p_0[p];
        else k_force_adapt_p_rot_x[p] /= k_force_adapt_p_0[p];

        k_force_adapt_mean += k_force_adapt_p_rot_x[p];

        if (dT_p_rot_y[p] > 0) k_force_adapt_p_rot_y[p] *= k_force_adapt_p_0[p];
        else k_force_adapt_p_rot_y[p] /= k_force_adapt_p_0[p];

        k_force_adapt_mean += k_force_adapt_p_rot_y[p];

        if (dT_p_rot_z[p] > 0) k_force_adapt_p_rot_z[p] *= k_force_adapt_p_0[p];
        else k_force_adapt_p_rot_z[p] /= k_force_adapt_p_0[p];

        k_force_adapt_mean += k_force_adapt_p_rot_z[p];
        
        //T_mean_loc_prev_revert_p[p] = T_mean_loc_prev_p[p];
        //T_mean_loc_prev_p[p] = T_mean_loc_p[p];
    }

    /*if (dT_p_x[p] > 0) v_r[p].x = (*var_nor)() * sqrt(kb * dT_p_x[p] / M0p[p]);
    else			   v_r[p].x = 0;

    if (dT_p_y[p] > 0) v_r[p].y = (*var_nor)() * sqrt(kb * dT_p_y[p] / M0p[p]);
    else			   v_r[p].y = 0;

    if (dT_p_z[p] > 0) v_r[p].z = (*var_nor)() * sqrt(kb * dT_p_z[p] / M0p[p]);
    else			   v_r[p].z = 0;

    if (dT_p_rot_x[p] > 0) w_r[p].x = (*var_nor)() * sqrt(kb * dT_p_rot_x[p] / I0p[p]);
    else			   w_r[p].x = 0;

    if (dT_p_rot_y[p] > 0) w_r[p].y = (*var_nor)() * sqrt(kb * dT_p_rot_y[p] / I0p[p]);
    else			   w_r[p].y = 0;

    if (dT_p_rot_z[p] > 0) w_r[p].z = (*var_nor)() * sqrt(kb * dT_p_rot_z[p] / I0p[p]);
    else			   w_r[p].z = 0;*/

    //k_bm_inst ++;

    //if ((T_mean > T) && (step % 10 == 0)) k_force_adapt /= k_force_adapt_0;
    //if ((T_mean < T) && (step % 10 == 0)) k_force_adapt *= k_force_adapt_0;
}

void ff_model_size_dispersion_init(void)
{
    //double d[14 + 1];
    double F[14 + 1];
    long i, p;
    long imax = 14;
    double Ftot;
    double random_points[14 + 1];
    double random_value;
    double large_fraction_tmp = 0;
    int is_set = 0;

    d[1] = 33.33333;
    F[1] = 0.0108843537;

    d[2] = 50;
    F[2] = 0.0795918367;

    d[3] = 66.6666;
    F[3] = 0.1850340136;

    d[4] = 83.333;
    F[4] = 0.1972789116;

    d[5] = 100;
    F[5] = 0.1925170068;

    d[6] = 116.66666;
    F[6] = 0.1217687075;

    d[7] = 133.3333;
    F[7] = 0.1006802721;

    d[8] = 150;
    F[8] = 0.074829932;

    d[9] = 166.666666;
    F[9] = 0.056462585;

    d[10] = 183.33333;
    F[10] = 0.0394557823;

    d[11] = 200;
    F[11] = 0.0278911565;

    d[12] = 216.66666;
    F[12] = 0.019047619;

    d[13] = 233.33333;
    F[13] = 0.0115646259;

    d[14] = 250;
    F[14] = 0.006122449;

	//for (i = 1; i <= imax; i++) d[i] = 200; // same size

    Ftot = 0;
    for (i = 1; i <= imax; i++)
    {
        d[i] *= 1E-10 * kr; // metric system
        Ftot += F[i];
    }

    large_fraction_tmp = 0;
    if (is_large_mode) for (i = imax; i >= 1; i--)
    {
        large_fraction_tmp += F[i] / Ftot;
        if (large_fraction_tmp > large_fraction)
        {
            i_min = i + 1;
            break;
        }
    }

    V0_tot_EV = V0_largest_EV = 0;
    // particle p = 1 parameters must be redefined below in next cycles
    for (i = 1; i <= imax; i++)
    {
        Rp0[1] = 0.5 * d[i];
        ff_model_size_dispersion_param_calc(Rp0[1], 1);
        V0_tot_EV += Vp0[1] * F[i] / Ftot;
        if (i >= i_min) V0_largest_EV += Vp0[1] * F[i] / Ftot;
    }

    if (is_large_mode) V0_tot_EV *= pN / large_fraction;
    else V0_tot_EV *= pN;

    if (is_large_mode) V0_largest_EV *= pN / large_fraction;
    else V0_largest_EV *= pN;

    for (i = 1; i <= imax; i++)
    {
        random_points[i] = F[i] / Ftot;
        if (i > 1) random_points[i] += random_points[i - 1];
    }
	
    //printf("\n random_points[14] = %e", random_points[14]);
    //printf("\n random_points[13] = %e", random_points[13]);

    V0_tot = 0;

    for (p = 1; p <= pN; p++)
    {
        is_set = 0;

        again_size_disp:
        random_value = (*var_uni)();

        //printf("\n random_value = %e", random_value);

        if ((random_value <= random_points[1]) && (i_min == 1))
        {
            if (is_large_mode) Rp0[p] = 0.5 * d[1] * k_large;
            else Rp0[p] = 0.5 * d[1];

            Rp[p] = Rp0[p] + delta;
            ff_model_size_dispersion_param_calc(Rp0[p], p);
            is_set = 1;

            V0_tot += Vp0[p];
        }
        
        for (i = 1; i <= imax - 1; i++)
            if ((random_value > random_points[i]) && (random_value <= random_points[i + 1]) && (i + 1 >= i_min))
            {
                if (is_large_mode) Rp0[p] = 0.5 * d[i + 1] * k_large;
                else Rp0[p] = 0.5 * d[i + 1];

                Rp[p] = Rp0[p] + delta;
                ff_model_size_dispersion_param_calc(Rp0[p], p);
                is_set = 1;

                V0_tot += Vp0[p];
                break;
            }
        if (is_set == 0) goto again_size_disp;

        // fixed size option
        /*Rp0[p] = 0.5 * d[6];
        Rp[p] = Rp0[p] + delta;
        ff_model_size_dispersion_param_calc(Rp0[p], p);*/
    }
}

void ff_model_size_dispersion_param_calc(double R0, long p)
{
    double d0;
    d0 = 2 * R0; // size of core, i.e. hard part of the particle: magnetic + nonmagnetic layers
    double V0 = pi * pow(d0, 3) / 6.0;
    double Vfull = pi * pow(d0 + 2 * delta, 3) / 6.0;
    double Vmag = pi * pow(d0 - 2 * a0, 3) / 6.0; // magentic core volume
    double tm = rop * V0; // total mass
    double tmmag = rop * Vmag; // mass of the magnetic part
    double theta, phi;

    if (2 * Rp0[p] > d_neel) is_neel[p] = 0;
    else is_neel[p] = 1;

    if (p == 1) R0_min = Rp0[p];
    if (Rp0[p] < R0_min) R0_min = Rp0[p];
    
    Vp0[p] = V0;
    Vpfull[p] = Vfull;
    M0p[p] = tm;
    I0p[p] = (2 / 5.0) * M0p[p] * R0 * R0;

    double m_mol_rel = 231.6;
    double m_mol = m_mol_rel * 1E-3 / Na;
    double N_mol = tmmag / m_mol;
    double s_mol = 4.1 * muB;
    double s = s_mol * N_mol;

    m0p[p] = s;

    if (t == 0)
    {
        theta = pi * (*var_uni)();
        phi = 2 * pi * (*var_uni)();

        m[p].x = m0p[p] * sin(theta) * cos(phi);
        m[p].y = m0p[p] * sin(theta) * sin(phi);
        m[p].z = m0p[p] * cos(theta);
    }
}

/*void ff_model_brownian_validation(long p)
{
double gamma = 6 * pi * eta_car * Rp[p];
double D = kb * T / gamma;
double dr_root_theory = 0;
double dr_root_sim = 0;

dr_root_sim = sqrt(pow(r[p].x - r_brown_valid_0.x, 2) + pow(r[p].y - r_brown_valid_0.y, 2) + pow(r[p].z - r_brown_valid_0.z, 2));
dr_root_theory = sqrt(6 * D * t);

//printf("\n ff_model_brownian_validation: %e %%", 100 * (dr_root_sim - dr_root_theory) / dr_root_theory);
}*/

void ff_model_update_conc_in_oleic(long p)
{
    //if (Rp_to_c[p] <= R_oleic)
        //if (Rp_to_c[p] <= Lx / 8.0)
    {
        is_inside_oleic[p] = 1;

        pN_oleic_drop++;
        if ((2 * Rp0[p] >= d[1] * 0.99) && (2 * Rp0[p] <= d[2] * 1.01)) pN_oleic_drop_I++;
        if ((2 * Rp0[p] >= d[3] * 0.99) && (2 * Rp0[p] <= d[7] * 1.01)) pN_oleic_drop_II++;
        if ((2 * Rp0[p] >= d[8] * 0.99) && (2 * Rp0[p] <= d[14] * 1.01)) pN_oleic_drop_III++;

        phi_vol_fract_oleic += Vp0[p];
        V0_tot_oleic += Vp0[p];
    }
    //else
    {
        //if (is_inside_oleic[p] == 1) k_force_adapt_p[p] = 1;
        is_inside_oleic[p] = 0;
    }
}

void ff_model_update_mdrop_parameters()
{
    long p = 0;
    double R_micro = 0;
    double n_p = 0;
    double mtot_p = 0;

    n_p = phi_v / (V0_tot_EV / pN);
    //printf("\n n_p = %e", n_p);

    R_micro = 2 * sigma_sf / (kb * T * n_p);
    //printf("\n R_micro = %e", R_micro);
    //printf("\n R_micro / R0 = %e", R_micro / (5E-9)); //2.6E4
    R0_min = R_micro;
    
    for (p = 1; p <= pN; p++)
    {
        is_neel[p] = 1; // to avoid false microdrop magnetization rotation inertia
        mtot_p = sqrt(MUL(m[p],m[p]));
        Rp[p] = Rp0[p] = R_micro;
        Vp0[p] = (4 * pi / 3.0) * pow(Rp0[p],3);
        m0p[p] = alpha * Ms * phi_v * Vp0[p];
        m[p].x *= m0p[p] / mtot_p;
        m[p].y *= m0p[p] / mtot_p;
        m[p].z *= m0p[p] / mtot_p;
        M0p[p] = (rop * phi_v + ro_oleic * (1 - phi_v)) * Vp0[p];
    }
}