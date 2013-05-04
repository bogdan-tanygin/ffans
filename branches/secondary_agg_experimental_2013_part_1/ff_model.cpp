/**************************************************************************
 * Copyright (C) 2011,2013 Dr. Bogdan Tanygin<b.m.tanygin@gmail.com>
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

#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

#include "ff_model.h"
#include "ff_model_parameters.h"
#include "ff_model_graphics.h"
#include "ff_model_io.h"

// working variables 
////////////////////
boost::mt19937 rng;
boost::normal_distribution<>* nd; // TODO need delete
boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >* var_nor;

ff_vect_t r[pN + 1];  //particles positions
ff_vect_t m[pN + 1];  //particles magnetic moment direction
int m_sat[pN + 1];

ff_vect_t F[pN + 1];

ff_vect_t v[pN + 1];

ff_vect_t drt[pN + 1];

ff_vect_t dvt[pN + 1];

//ff_vect_t dir110[13];

long number_of_mag_threads = 0;
long mag_thread_limit_plus[pN + 1]; // pN is a maximum number of threads. Index is id of thread.
long mag_thread_limit_minus[pN + 1]; // pN is a maximum number of threads. Index is id of thread.
long mag_thread_p_id[pN + 1];
long mag_thread_id_ordered_p[pN + 1][pN + 1]; // [thread_id][ordered position of particle inside thread] = p
long mag_thread_p_position[pN + 1];
//double **A_solve; 

double Rp[pN + 1];
double m0p[pN + 1];
double m0start_agg_p[pN + 1];
int exist_p[pN + 1]; // particle existence; number of primary aggregate inside
double r0modp[pN + 1];
double Vselfp[pN + 1];

int time_go = 0;

double B_hyst[21];
double Mz_hyst[21];
double B_hyst_n[21];
double Mz_hyst_n[21];

double t; // time
double Ek;
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

        ff_io_save_hyst();

        mz_glob = 0;
        glob_start_t = t;
        glob_start_step = step;

    }
}

int ff_model_check_smooth_dr(long p)
{
    double dr, rmod;
    long ps;
    int res = 1;

    dr = sqrt(MUL(drt[p], drt[p]));
    rmod = 2 * R0;

    if (rmod > 0)
        if (dr / rmod > smooth_r)
        {
            for(ps = 1; ps <= p; ps ++)
            {
                r[ps].x -= drt[ps].x;
                r[ps].y -= drt[ps].y;
                r[ps].z -= drt[ps].z;
            }
            dt /= 2.0;
            slow_steps += 10;
            res = 0;
            step--;
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

void ff_model_set_rand_dir(long p)
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
    */
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
}

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

void ff_model_m_setting(void)
{
    long p;
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

    long thread_p, thread_ps, thread_pos_ps, p_prev;

    ff_vect_t tau;

    double sign;
    int is_last;

    //long di, dj;

    ff_vect_t ttB;

    //	do
    //	{
    for(p = 1; p <= pN; p ++)
        if (exist_p[p])
        {
            ttB = Bext(r[p].x, r[p].y, r[p].z);
            tBx = ttB.x;
            tBy = ttB.y;
            tBz = ttB.z;

            for(ps = 1; ps <= pN; ps ++)
                if (exist_p[ps])
                    if (p != ps)
                    {
                        /*#ifdef SECONDARY
                        for(di = -2; di <= 2; di++)
                        for(dj = -2; dj <= 2; dj++) // including zero!
                        {
                        #endif*/
                        dx = - r[ps].x + r[p].x;// + di * Lx;
                        dy = - r[ps].y + r[p].y;// + dj * Ly;
                        dz = - r[ps].z + r[p].z;

                        dR = sqrt(dx * dx + dy * dy + dz * dz);
#ifdef SECONDARY
                        if (dR <= Rp[p] + Rp[ps])
                        {
                            tk = (Rp[p] + Rp[ps]) / dR;
                            dx *= tk;
                            dy *= tk;
                            dz *= tk;
                            dR = sqrt(dx * dx + dy * dy + dz * dz);
                        }
#endif

                        dR2 = dR * dR;

                        thread_p = mag_thread_p_id[p];
                        thread_ps = mag_thread_p_id[ps];

                        if (1/*(thread_p == thread_ps) || (thread_ps == 0)*/)
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

                            if (dtBx != dtBx) printf("\n DEBUG 12 p = %d ps = %d dtBx = %e r[p].x = %e r[ps].x = %e mxs = %e mys = %e mzs = %e", p, ps, dtBx, r[p].x, r[ps].x, mxs, mys, mzs);

                        } // --- if ((!mag_threads) || (thread_p == thread_ps) || (thread_ps == 0))
                        else
                        {
                            /* thread_p = mag_thread_p_id[p];
                            thread_ps = mag_thread_p_id[ps];

                            //if ((thread_p != thread_ps) && (thread_ps != 0))
                            //{
                            dtBx = dtBy = dtBz = 0;

                            thread_pos_ps = mag_thread_p_position[ps];
                            is_last = 0;
                            if (thread_pos_ps == mag_thread_limit_plus[thread_ps]  ) {p_prev = mag_thread_id_ordered_p[thread_ps][thread_pos_ps - 1];is_last = 1;}
                            if (thread_pos_ps == mag_thread_limit_minus[thread_ps] ) {p_prev = mag_thread_id_ordered_p[thread_ps][thread_pos_ps + 1];is_last = 1;}

                            if (is_last)
                            {
                            tau.x = r[ps].x - r[p_prev].x;
                            tau.y = r[ps].y - r[p_prev].y;
                            tau.z = r[ps].z - r[p_prev].z;

                            sign = MUL(tau,m[ps]);
                            if (sign > 0) sign = 1;
                            else sign = -1;

                            tms = sqrt(MUL(m[ps],m[ps]));

                            mag_qs = tms / (4.0 * Rp[ps] / 3.0);

                            dtBx = sign * C5 * mag_qs * dx / pow(dR, 3);
                            dtBy = sign * C5 * mag_qs * dy / pow(dR, 3);
                            dtBz = sign * C5 * mag_qs * dz / pow(dR, 3);
                            }
                            //}	*/   
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

                        tBx += dtBx;
                        tBy += dtBy;
                        tBz += dtBz;

                        /*#ifdef SECONDARY
                        } // di and dj
                        #endif*/
                    } // if (p != ps)

                    tBmag = sqrt(tBx * tBx + tBy * tBy + tBz * tBz);
                    /*tm = sqrt(MUL(m[p],m[p]));

                    eps = acos((m[p].x * tBx + m[p].y * tBy + m[p].z * tBz) / (tBmag * tm));

                    if (p == 1) max_eps = eps;
                    else if (eps > max_eps) max_eps = eps;*/

#ifndef SECONDARY
                    m[p].x = m0 * tBx / tBmag;
                    m[p].y = m0 * tBy / tBmag;
                    m[p].z = m0 * tBz / tBmag;
#endif

#ifdef SECONDARY
                    //if (   ( (tBmag / mu0) >= H_max_agg) || (m_freeze[p])   )
                    if ((tBmag / mu0) >= H_max_agg)
                    {
                        m[p].x = m0p[p] * tBx / tBmag;
                        m[p].y = m0p[p] * tBy / tBmag;
                        m[p].z = m0p[p] * tBz / tBmag;
                        m_sat[p] = 1;
                        //m_freeze[p] = 1;
                        //printf("\n (tBmag / mu0) >= H_max_agg, p = %d", p);
                    }
                    else
                    {
                        m[p].x = Vselfp[p] * hi_m_agg * tBx / mu0;
                        m[p].y = Vselfp[p] * hi_m_agg * tBy / mu0;
                        m[p].z = Vselfp[p] * hi_m_agg * tBz / mu0;

                        tm = sqrt(MUL(m[p],m[p]));
                        m_sat[p] = 0;

                        //if (100 * tm / m0p[p] > 30)
                        //	printf("\n REGULAR, p = %d, saturation_percentage = %f", p, 100 * tm / m0p[p]);

                        if (tm >= m0p[p])
                        {
                            m[p].x = m0p[p] * tBx / tBmag;
                            m[p].y = m0p[p] * tBy / tBmag;
                            m[p].z = m0p[p] * tBz / tBmag;
                            m_sat[p] = 1;
                            //m_freeze[p] = 1;
                            //printf("\n tm >= m0, p = %d", p);
                        }
                        if (tm <= m0start_agg_p[p])
                        {
                            m[p].x = m0start_agg_p[p] * tBx / tBmag;
                            m[p].y = m0start_agg_p[p] * tBy / tBmag;
                            m[p].z = m0start_agg_p[p] * tBz / tBmag;
                            m_sat[p] = 0;
                            //printf("\n tm <= m0start_agg, p = %d", p);
                        }
                    }
#endif

                    if (m[p].x != m[p].x) printf("\n DEBUG 11 p = %d m[p].x = %e", p, m[p].x);
                    if (m[p].y != m[p].y) printf("\n DEBUG 11 p = %d m[p].y = %e", p, m[p].y);
                    if (m[p].z != m[p].z) printf("\n DEBUG 11 p = %d m[p].z = %e", p, m[p].z);

        }// for(p = 1; p <= pN; p ++)
        //	}
        //	while(max_eps > m_h_eff_tol);
}

/*void ff_model_pp_phagocyt(long p1, long p2)
{
long size_new, size_old; // size [number of primary aggregates inside]

size_old = exist_p[p1];
size_new = exist_p[p1] += exist_p[p2];
exist_p[p2] = 0;

Rp[p1] = Rp[p1] * pow(size_new / size_old, 1 / 3.0);
m[p1].x = m[p1].x * size_new / size_old;
m[p1].y = m[p1].y * size_new / size_old;
m[p1].z = m[p1].z * size_new / size_old;
m0p[p1] = m0p[p1] * size_new / size_old;
m0start_agg_p[p1] *= size_new / size_old;
r0modp[p1] *= pow(size_new / size_old, - 1 / 6.0);

Vselfp[p1] *= size_new / size_old;
}*/

ff_vect_t ff_model_nonloc_force(long p)
{
    register long ps;
    double tFx, tFy, tFz;
    double dtFx, dtFy, dtFz;
    double dreptFx, dreptFy, dreptFz; // repulsion
    double dR, dR2, dR3, dR2mod, dR5, dR5mod, MUL1, MUL2, MUL3;
    //double R1 = 0.3 * R0;
    double Cmod;
    double Cmod1 = Ch * m0 * m0;
    double Cmod1_repulsion = 5 * m0 * m0;
    double dx, dy, dz;
    double tk;
    double sec_pow;
    ff_vect_t ttF;
    long thread_id;
    long thread_length_p, thread_length_ps, tp;
    long new_thread_length;

    double MUL1__dR5mod, MUL2__dR5mod, MUL3__dR5mod, MUL1_MUL2__dR2mod__dR5mod;

    double mx, my, mz;
    double mxs, mys, mzs;
    long thread_p, thread_ps, thread_pos_ps, p_prev;

    int is_last;
    double sign, tms, mag_qs;

    ff_vect_t tau;

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

                dx = r[ps].x - r[p].x;// + di * Lx;
                dy = r[ps].y - r[p].y;// + dj * Ly;
                dz = r[ps].z - r[p].z;

                mx = m[p].x;
                my = m[p].y;
                mz = m[p].z;

                mxs = m[ps].x;
                mys = m[ps].y;
                mzs = m[ps].z;

                //dipole-dipole
                dR = sqrt(dR2 = dx * dx + dy * dy + dz * dz);
                dR3 = dR2 * dR;

                //if (dR == 0) printf ("\n \n \n dR == 0 !!!");

                // Non local force and self-confid. magnetic field
                dR2mod = dR2 / 5.0; //modified
                dR5mod = (dR5 = pow(dR,5)) / C1; //modified

                MUL2 = mxs * dx + mys * dy + mzs * dz;

                //#ifndef SECONDARY
                //if (p < ps)  !!!!!
                //{

                thread_p = mag_thread_p_id[p];
                thread_ps = mag_thread_p_id[ps];

                if ((thread_p == thread_ps) || (thread_ps == 0))
                {
                    //#endif
                    MUL1 = mx * dx + my * dy + mz * dz;
                    MUL3 = mx * mxs + my * mys + mz * mzs;

                    MUL1__dR5mod = MUL1 / dR5mod;
                    MUL2__dR5mod = MUL2 / dR5mod;
                    MUL3__dR5mod = MUL3 / dR5mod;
                    MUL1_MUL2__dR2mod__dR5mod = MUL1 * MUL2 / (dR2mod * dR5mod);

                    dtFx = -(MUL1__dR5mod * mxs + MUL2__dR5mod * mx
                        + MUL3__dR5mod * dx - MUL1_MUL2__dR2mod__dR5mod * dx);

                    dtFy = -(MUL1__dR5mod * mys + MUL2__dR5mod * my
                        + MUL3__dR5mod * dy - MUL1_MUL2__dR2mod__dR5mod * dy);

                    dtFz = -(MUL1__dR5mod * mzs + MUL2__dR5mod * mz
                        + MUL3__dR5mod * dz - MUL1_MUL2__dR2mod__dR5mod * dz);

                    //if (p == 50) printf("\n DEBUG 7 ps = %d F.x = %e", ps, dtFx);
                    //if (p == 50) printf("\n DEBUG 7 ps = %d m[p].x = %e", ps, m[p].x);
                    //if (p == 50) printf("\n DEBUG 7 ps = %d m[ps].x = %e", ps, m[ps].x);

                } // --- if ((!mag_threads) || (thread_p == thread_ps) || (thread_ps == 0))
                else 
                {
                    dtFx = dtFy = dtFz = 0;

                    thread_pos_ps = mag_thread_p_position[ps];
                    is_last = 0;
                    if (thread_pos_ps == mag_thread_limit_plus[thread_ps]  ) {p_prev = mag_thread_id_ordered_p[thread_ps][thread_pos_ps - 1];is_last = 1;}
                    if (thread_pos_ps == mag_thread_limit_minus[thread_ps] ) {p_prev = mag_thread_id_ordered_p[thread_ps][thread_pos_ps + 1];is_last = 1;}

                    if (is_last)
                    {
                        tau.x = r[ps].x - r[p_prev].x;
                        tau.y = r[ps].y - r[p_prev].y;
                        tau.z = r[ps].z - r[p_prev].z;

                        sign = MUL(tau,m[ps]);
                        if (sign > 0) sign = 1;
                        else sign = -1;

                        tms = sqrt(MUL(m[ps],m[ps]));

                        mag_qs = tms / (4.0 * Rp[ps] / 3.0);

                        dx *= -1; dy *= -1; dz *= -1; // to make "- r[ps].x + r[p].x"

                        dtFx = sign * C5 * (mag_qs / dR3) * (m[p].x * (1 - 3 * dx * dx / dR2) - m[p].y * 3 * dx  * dy - m[p].z * 3 * dx  * dz);
                        dtFy = sign * C5 * (mag_qs / dR3) * (m[p].y * (1 - 3 * dy * dy / dR2) - m[p].x * 3 * dy  * dx - m[p].z * 3 * dy  * dz);
                        dtFz = sign * C5 * (mag_qs / dR3) * (m[p].z * (1 - 3 * dz * dz / dR2) - m[p].x * 3 * dz  * dx - m[p].y * 3 * dz  * dy);
                    }
                }

                tFx += dtFx;
                tFy += dtFy;
                tFz += dtFz;

                /*
                if (dR <= 5 * R0) // exp phenomenology
                {
                tFx += - (EPS * dx / dR) * (exp(-(dR - 2 * R00) / ro1) / ro1 - exp(-(dR - 2 * R00) / ro2) / ro2);
                tFy += - (EPS * dy / dR) * (exp(-(dR - 2 * R00) / ro1) / ro1 - exp(-(dR - 2 * R00) / ro2) / ro2);
                tFz += - (EPS * dz / dR) * (exp(-(dR - 2 * R00) / ro1) / ro1 - exp(-(dR - 2 * R00) / ro2) / ro2);
                }
                */
                // acid elasticity (repulsion)
#ifndef SECONDARY
                if (dR <= 2 * R0 ) // the Heaviside step function  and dR5 dependence finally is similar to the well-known exp. phenomenology
                {
                    if (Ch > 5)
                        Cmod = Cmod1 * (C1 / dR5);
                    else
                        Cmod = Cmod1_repulsion * (C1 / dR5);

                    tFx += -dx * Cmod;
                    tFy += -dy * Cmod;
                    tFz += -dz * Cmod;
                }
#endif

#ifdef SECONDARY
                dreptFx = dreptFy = dreptFz = 0;

                if (mag_threads && (mag_thread_p_id[p] == mag_thread_p_id[ps]) && (fabs(0.0 + mag_thread_p_position[p] - mag_thread_p_position[ps]) == 1.0)) // model adhesion
                {
                    if (dR > 1.3 * (Rp[p] + Rp[ps]) )
                    {
                        Cmod = - sqrt(MUL(m[p],m[p])) * sqrt(MUL(m[ps],m[ps])) * (C1 / (pow(Rp[p] + Rp[ps], 5)));
                        tk = (Rp[p] + Rp[ps]) / dR;

                        dreptFx = - tk * dx * Cmod;
                        dreptFy = - tk * dy * Cmod;
                        dreptFz = - tk * dz * Cmod;

                        //if (p == 50) printf("\n DEBUG 8 ps = %d F.x = %e", ps, dreptFx);
                    }
                }

                if (dR <= 1.1 * (Rp[p] + Rp[ps]) ) // repulsion: the Heaviside step function  and dR5 dependence finally is similar to the well-known exp. phenomenology
                {
                    if (dR <= Rp[p] + Rp[ps] )
                    {
                        Cmod = sqrt(MUL(m[p],m[p])) * sqrt(MUL(m[ps],m[ps])) * (C1 / (pow(Rp[p] + Rp[ps], 4) * dR)); // model repulsion
                        tk = (Rp[p] + Rp[ps]) / dR;

                        dreptFx = - tk * dx * Cmod - dtFx; // disable magnetic attraction only for given pair: model trick
                        dreptFy = - tk * dy * Cmod - dtFy;
                        dreptFz = - tk * dz * Cmod - dtFz;

                        //if (p == 50) printf("\n DEBUG 8 ps = %d F.x = %e", ps, dreptFx);
                    }

                    // magnetic threads connections
                    ////////////////////////////////
                    if (1 && (mag_threads) && ((mag_thread_p_id[p] != 0) && (mag_thread_p_id[ps] != 0) && (mag_thread_p_id[p] != mag_thread_p_id[ps])) )
                    {
                        thread_p = mag_thread_p_id[p];
                        thread_ps = mag_thread_p_id[ps];
                        thread_length_p =  mag_thread_limit_plus[thread_p] - mag_thread_limit_minus[thread_p] + 1;
                        thread_length_ps = mag_thread_limit_plus[thread_ps] - mag_thread_limit_minus[thread_ps] + 1;

                        if ((mag_thread_p_position[ps] == mag_thread_limit_plus[thread_ps]) && (mag_thread_p_position[p] == mag_thread_limit_plus[thread_p])
                            && (mag_thread_id_ordered_p[thread_p][0] != -1) && (mag_thread_id_ordered_p[thread_ps][0] != -1) && (mag_thread_p_id[p] != mag_thread_p_id[ps]))
                        {
                            number_of_mag_threads++;
                            new_thread_length = thread_length_p + thread_length_ps;
                            printf("\n DEBUG 1 %d",number_of_mag_threads);
                            for (long pos = 1; pos <= thread_length_ps; pos++)  // from minus to plus of ps
                            {
                                tp = mag_thread_id_ordered_p[thread_ps][pos];
                                mag_thread_p_id[tp] = number_of_mag_threads;
                                // mag_thread_p_position[tp] // - the same
                            }

                            mag_thread_id_ordered_p[thread_ps][0] = -1; // flag of removed thread

                            for (long pos = thread_length_p; pos >= 1; pos--)  // from plus to minus of p
                            {
                                tp = mag_thread_id_ordered_p[thread_p][pos];
                                mag_thread_p_id[tp] = number_of_mag_threads;
                                mag_thread_p_position[tp] = mag_thread_limit_plus[thread_ps] + thread_length_p - pos + 1;
                            }

                            mag_thread_id_ordered_p[thread_p][0] = -1; // flag of removed thread

                            mag_thread_limit_minus[number_of_mag_threads] = mag_thread_limit_minus[thread_ps];
                            mag_thread_limit_plus[number_of_mag_threads] = mag_thread_p_position[tp]; // latest tp from loop

                        }
                        /////////////////////////////////////////
                        if ((mag_thread_p_position[ps] == mag_thread_limit_plus[thread_ps]) && (mag_thread_p_position[p] == mag_thread_limit_minus[thread_p])
                            && (mag_thread_id_ordered_p[thread_p][0] != -1) && (mag_thread_id_ordered_p[thread_ps][0] != -1) && (mag_thread_p_id[p] != mag_thread_p_id[ps]))
                        {
                            number_of_mag_threads++;
                            new_thread_length = thread_length_p + thread_length_ps;
                            printf("\n DEBUG 2 %d",number_of_mag_threads);
                            for (long pos = 1; pos <= thread_length_ps; pos++)  // from minus to plus of ps
                            {
                                tp = mag_thread_id_ordered_p[thread_ps][pos];
                                mag_thread_p_id[tp] = number_of_mag_threads;
                                // mag_thread_p_position[tp] // - the same
                            }

                            mag_thread_id_ordered_p[thread_ps][0] = -1; // flag of removed thread

                            for (long pos = 1; pos <= thread_length_p; pos++)
                            {
                                tp = mag_thread_id_ordered_p[thread_p][pos];
                                mag_thread_p_id[tp] = number_of_mag_threads;
                                mag_thread_p_position[tp] = mag_thread_limit_plus[thread_ps] + pos;
                            }

                            mag_thread_id_ordered_p[thread_p][0] = -1; // flag of removed thread

                            mag_thread_limit_minus[number_of_mag_threads] = mag_thread_limit_minus[thread_ps];
                            mag_thread_limit_plus[number_of_mag_threads] = mag_thread_p_position[tp]; // latest tp from loop

                        }
                        /////////////////////////////////////////
                        if ((mag_thread_p_position[ps] == mag_thread_limit_minus[thread_ps]) && (mag_thread_p_position[p] == mag_thread_limit_minus[thread_p])
                            && (mag_thread_id_ordered_p[thread_p][0] != -1) && (mag_thread_id_ordered_p[thread_ps][0] != -1) && (mag_thread_p_id[p] != mag_thread_p_id[ps]))
                        {
                            number_of_mag_threads++;
                            new_thread_length = thread_length_p + thread_length_ps;
                            printf("\n DEBUG 3 %d",number_of_mag_threads);
                            for (long pos = thread_length_ps; pos >= 1; pos--)  // from minus to plus of ps
                            {
                                tp = mag_thread_id_ordered_p[thread_ps][pos];
                                mag_thread_p_id[tp] = number_of_mag_threads;
                                mag_thread_p_position[tp] =  - mag_thread_limit_plus[thread_ps] + thread_length_ps - pos;
                            }

                            mag_thread_id_ordered_p[thread_ps][0] = -1; // flag of removed thread

                            for (long pos = thread_length_p; pos >= 1; pos--)  // from plus to minus of p
                            {
                                tp = mag_thread_id_ordered_p[thread_p][pos];
                                mag_thread_p_id[tp] = number_of_mag_threads;
                                mag_thread_p_position[tp] = - mag_thread_limit_plus[thread_ps] + thread_length_ps - 1 + thread_length_p - pos + 1;
                            }

                            mag_thread_id_ordered_p[thread_p][0] = -1; // flag of removed thread

                            mag_thread_limit_minus[number_of_mag_threads] = - mag_thread_limit_plus[thread_ps];
                            mag_thread_limit_plus[number_of_mag_threads] = mag_thread_p_position[tp]; // latest tp from loop

                        }
                        ///////////////////////////////////////////

                    } // --- if ((!mag_threads) || ((mag_thread_p_id[p] != 0) && (mag_thread_p_id[ps] != 0) && (mag_thread_p_id[p] != mag_thread_p_id[ps])) )

                    tFx += dreptFx;
                    tFy += dreptFy;
                    tFz += dreptFz;

                    //if (p == 50) printf("\n DEBUG 9 ps = %d F.x = %e", ps, tFx);

                    if ((mag_thread_p_id[p] != 0) && (mag_thread_p_id[ps] == 0))
                    {
                        thread_id = mag_thread_p_id[p];
                        if (mag_thread_limit_plus[thread_id] == mag_thread_p_position[p])
                        {
                            mag_thread_p_id[ps] = thread_id;
                            mag_thread_p_position[ps] = ++(mag_thread_limit_plus[thread_id]);
                        }

                        if (mag_thread_limit_minus[thread_id] == mag_thread_p_position[p])
                        {
                            mag_thread_p_id[ps] = thread_id;
                            mag_thread_p_position[ps] = --(mag_thread_limit_minus[thread_id]);
                        } 
                    }

                    if ((mag_thread_p_id[p] == 0) && (mag_thread_p_id[ps] != 0))
                    {
                        thread_id = mag_thread_p_id[ps];
                        if (mag_thread_limit_plus[thread_id] == mag_thread_p_position[ps])
                        {
                            mag_thread_p_id[p] = thread_id;
                            mag_thread_p_position[p] = ++(mag_thread_limit_plus[thread_id]);
                        } 

                        if (mag_thread_limit_minus[thread_id] == mag_thread_p_position[ps])
                        {
                            mag_thread_p_id[p] = thread_id;
                            mag_thread_p_position[p] = --(mag_thread_limit_minus[thread_id]);
                        } 
                    }

                    if ((mag_thread_p_id[p] == 0) && (mag_thread_p_id[ps] == 0))
                    {
                        number_of_mag_threads++;
                        mag_thread_p_id[p] = mag_thread_p_id[ps] = number_of_mag_threads;
                        mag_thread_p_position[p] = 0;
                        mag_thread_p_position[ps] = 1;
                        mag_thread_limit_plus[number_of_mag_threads] = 1;
                        mag_thread_limit_minus[number_of_mag_threads] = 0;
                    }

                    /*if ((exist_p[p])&&(exist_p[ps]))
                    if ((m_sat[p])||(m_sat[ps]))
                    {
                    if (Rp[p] >= Rp[ps])
                    {
                    printf("\n === 1, N1=%d, N2=%d", exist_p[p], exist_p[ps]);
                    ff_model_pp_phagocyt(p, ps);
                    }
                    if (Rp[p] < Rp[ps])
                    {
                    printf("\n === 2, N1=%d, N2=%d", exist_p[p], exist_p[ps]);
                    ff_model_pp_phagocyt(ps, p);
                    }
                    }*/
                }

                // adhesion of secondary clusters
                /*if ((dR > Rp[p] + Rp[ps]) && (dR <= 2 * R0 + (Rp[p] + Rp[ps]))) // the Heaviside step function  and dR5 dependence finally is similar to the well-known exp. phenomenology
                {
                Cmod = sqrt(MUL(m[p],m[p])) * sqrt(MUL(m[ps],m[ps])) * (Ch * C1 / (pow(Rp[p] + Rp[ps], 5)));
                tk = (Rp[p] + Rp[ps]) / dR;

                tFx +=   tk * dx * Cmod;
                tFy +=   tk * dy * Cmod;
                tFz +=   tk * dz * Cmod;
                }*/
#endif

#ifndef SECONDARY
                // attraction
                if ((dR > 2 * R0 )&&(dR < 3 * R0 )) // the Heaviside step function  and dR5 dependence finally is similar to the well-known exp. phenomenology
                {
                    Cmod = Cmod1 * (C1 / dR5);

                    tFx += dx * Cmod;
                    tFy += dy * Cmod;
                    tFz += dz * Cmod;
                }
#endif

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
                //tBx += dx * MUL2mod__dR5mod1 - mxs * dR2__dR5mod1;
                //tBy += dy * MUL2mod__dR5mod1 - mys * dR2__dR5mod1;
                //tBz += dz * MUL2mod__dR5mod1 - mzs * dR2__dR5mod1;
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

            // force
            ttF.x = tFx;
            ttF.y = tFy;
            ttF.z = tFz;

            //if (ttF.x != ttF.x) printf("\n DEBUG 10 p = %d F.x = %e", p, ttF.x);
            //if (ttF.y != ttF.y) printf("\n DEBUG 10 p = %d F.y = %e", p, ttF.y);
            //if (ttF.z != ttF.z) printf("\n DEBUG 10 p = %d F.z = %e", p, ttF.z);

            return ttF;
}

ff_vect_t ff_model_brownian_force(long p)
{
    double theta, phi;
    double Fmag;
    ff_vect_t tF;

    Fmag = (*var_nor)();
    phi = 2 * pi * rand() / 32768.0;
    theta = pi * rand() / 32768.0;

    tF.x = Fmag * sin(theta) * cos(phi);
    tF.y = Fmag * sin(theta) * sin(phi);
    tF.z = Fmag * cos(theta);

    return tF;
}

ff_vect_t ff_model_force(long p)
{
    ff_vect_t tF;
    ff_vect_t tddF;

    tF.x = tF.y = tF.z = 0;

    // non-local
    tddF = ff_model_nonloc_force(p);
    tF.x += tddF.x;
    tF.y += tddF.y;
    tF.z += tddF.z;

    // Gravitation
    tF.z += - C3;

    // Buoyancy force
    tF.z +=   C6;

    // Random Gaussian process
    if (brownian_force) tddF = ff_model_brownian_force(p);
    tF.x += tddF.x;
    tF.y += tddF.y;
    tF.z += tddF.z;

    return tF;
}

int ff_model_check_walls(long p)
{
    int res = 0;
    double dr;

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

    return res;
}

/*void ff_model_check_overlapp(void)
{
long p, ps;
double dR, dx, dy, dz;

for (p = 1; p <= pN; p++)
for(ps = 1; ps <= pN; ps ++)
if (p != ps)
{
dx = r[ps].x - r[p].x;
dy = r[ps].y - r[p].y;
dz = r[ps].z - r[p].z;

dR = sqrt(dx * dx + dy * dy + dz * dz);

if (dR <= 2.01 * R0 )
{
r[ps].x -= drt[ps].x;
r[ps].y -= drt[ps].y;
r[ps].z -= drt[ps].z;
}
}
}*/

void ff_mag_threads_update_force(void)
{
    MatrixXd A_solve;
    VectorXd b_solve;
    VectorXd x_solve;

    long i, j, k, tp, tp_prev, pos;
    long thread_id, thread_length;
    double tau_x[pN + 1], tau_y[pN + 1], tau_z[pN + 1];
    double det_A, det_a_x, det_a_y, det_a_z, det_T, det_T_prev, tau_mod;
    double v_tau_1, v_tau_2, F_tau_1, F_tau_2;

    // --used-- long mag_thread_id_ordered_p[pN + 1][pN + 1]; // [thread_id][ordered position of particle inside thread] = p
    for (thread_id = 1; thread_id <= number_of_mag_threads; thread_id++)
        for (long p = 1; p <= pN; p++)
            if (mag_thread_p_id[p] == thread_id)
            {
                pos = mag_thread_p_position[p] - mag_thread_limit_minus[thread_id] + 1;
                mag_thread_id_ordered_p[thread_id][pos] = p;
            }	

            if (number_of_mag_threads)
            {
                for (thread_id = 1; thread_id <= number_of_mag_threads; thread_id++)
                {
                    thread_length = mag_thread_limit_plus[thread_id] - mag_thread_limit_minus[thread_id] + 1;
                    if ((thread_length >= 2) && (mag_thread_id_ordered_p[thread_id][0] != -1)) // -1 is removed thread after threads connection
                    {
                        k = 1;

                        A_solve.resize(4 * thread_length - 1, 4 * thread_length - 1);
                        b_solve.resize(4 * thread_length - 1);
                        x_solve.resize(4 * thread_length - 1);
                        A_solve.setZero();
                        b_solve.setZero();

                        for (pos = 1; pos <= thread_length; pos++)
                        {   
                            tp = mag_thread_id_ordered_p[thread_id][pos];

                            if (v[tp].x != v[tp].x) printf("\n DEBUG 1.3 p = %d r[p].x = %e v[p].x = %e", tp, r[tp].x, v[tp].x);
                            if (v[tp].y != v[tp].y) printf("\n DEBUG 1.3 p = %d r[p].y = %e v[p].y = %e", tp, r[tp].y, v[tp].y);
                            if (v[tp].z != v[tp].z) printf("\n DEBUG 1.3 p = %d r[p].z = %e v[p].z = %e", tp, r[tp].z, v[tp].z);

                            b_solve(3 * pos - 2/**/- 1) = F[tp].x / M0;
                            b_solve(3 * pos - 1/**/- 1) = F[tp].y / M0;
                            b_solve(3 * pos    /**/- 1) = F[tp].z / M0;

                            if (pos < thread_length) b_solve(3 * thread_length + pos /**/- 1) = 0;

                            A_solve(3 * k /**/- 1, 3 * k /**/- 1) = 1;  // normalized /M
                            A_solve(3 * k - 1 /**/- 1, 3 * k - 1 /**/- 1) = 1;
                            A_solve(3 * k - 2 /**/- 1, 3 * k - 2 /**/- 1) = 1;

                            if (k > 1)
                            {
                                tau_x[k - 1] = r[tp].x - r[tp_prev].x;
                                tau_y[k - 1] = r[tp].y - r[tp_prev].y;
                                tau_z[k - 1] = r[tp].z - r[tp_prev].z;

                                tau_mod = sqrt(tau_x[k - 1] * tau_x[k - 1] + tau_y[k - 1] * tau_y[k - 1] + tau_z[k - 1] * tau_z[k - 1]);

                                //printf("\n DEBUG ! tau_mod = %e", tau_mod);

                                tau_x[k - 1] /= tau_mod;
                                tau_y[k - 1] /= tau_mod;
                                tau_z[k - 1] /= tau_mod;

                                // align velocities - BEGIN
                                ///////////////////////////

                                v_tau_1 = (v[tp_prev].x * tau_x[k - 1] + v[tp_prev].y * tau_y[k - 1] + v[tp_prev].z * tau_z[k - 1]);
                                v_tau_2 = (v[tp].x * tau_x[k - 1] + v[tp].y * tau_y[k - 1] + v[tp].z * tau_z[k - 1]);

                                v[tp].x -= (v_tau_2 - v_tau_1) * tau_x[k - 1];
                                v[tp].y -= (v_tau_2 - v_tau_1) * tau_y[k - 1];
                                v[tp].z -= (v_tau_2 - v_tau_1) * tau_z[k - 1];

                                if (v[tp].x != v[tp].x) printf("\n DEBUG 1.4 p = %d tp_prev = %d r[p].x = %e v[p].x = %e taux = %e", tp, tp_prev, r[tp].x, v[tp].x, tau_x[k - 1]);
                                if (v[tp].y != v[tp].y) printf("\n DEBUG 1.4 p = %d tp_prev = %d r[p].y = %e v[p].y = %e tauy = %e", tp, tp_prev, r[tp].y, v[tp].y, tau_y[k - 1]);
                                if (v[tp].z != v[tp].z) printf("\n DEBUG 1.4 p = %d tp_prev = %d r[p].z = %e v[p].z = %e tauz = %e", tp, tp_prev, r[tp].z, v[tp].z, tau_z[k - 1]);

                                // align velocities - END
                                ///////////////////////////
                                if (k == thread_length)
                                {
                                    A_solve(3 * k /**/- 1,     3 * thread_length + k - 1 /**/- 1) = tau_z[k - 1];
                                    A_solve(3 * k - 1 /**/- 1, 3 * thread_length + k - 1 /**/- 1) = tau_y[k - 1];
                                    A_solve(3 * k - 2 /**/- 1, 3 * thread_length + k - 1 /**/- 1) = tau_x[k - 1];
                                    //printf("\n DEBUG ! tau_mod = %e", tau_mod);
                                }

                                k--; // !!! previous particle related k
                                A_solve(3 * k /**/- 1, 3 * thread_length + k /**/- 1)     = -tau_z[k];
                                A_solve(3 * k - 1 /**/- 1, 3 * thread_length + k /**/- 1) = -tau_y[k];
                                A_solve(3 * k - 2 /**/- 1, 3 * thread_length + k /**/- 1) = -tau_x[k];

                                if (k > 1)
                                {
                                    A_solve(3 * k /**/- 1, 3 * thread_length + k - 1 /**/- 1)     = tau_z[k - 1];
                                    A_solve(3 * k - 1 /**/- 1, 3 * thread_length + k - 1 /**/- 1) = tau_y[k - 1];
                                    A_solve(3 * k - 2 /**/- 1, 3 * thread_length + k - 1 /**/- 1) = tau_x[k - 1];
                                }
                                k++; // !!! previous particle related k
                            }

                            if (k == thread_length)
                            {
                                k = k - 1; // !!! previous particle related k
                                A_solve(3 * thread_length + k /**/- 1, 3 * k /**/- 1)     = tau_z[k];
                                A_solve(3 * thread_length + k /**/- 1, 3 * k - 1 /**/- 1) = tau_y[k];
                                A_solve(3 * thread_length + k /**/- 1, 3 * k - 2 /**/- 1) = tau_x[k];

                                A_solve(3 * thread_length + k /**/- 1, 3 * (k + 1) /**/- 1)     = -tau_z[k];
                                A_solve(3 * thread_length + k /**/- 1, 3 * (k + 1) - 1 /**/- 1) = -tau_y[k];
                                A_solve(3 * thread_length + k /**/- 1, 3 * (k + 1) - 2 /**/- 1) = -tau_x[k];
                                k = k + 1; // !!! previous particle related k
                            }

                            if (k > 2)
                            {
                                k = k - 2; // !!! previous particle related k
                                A_solve(3 * thread_length + k /**/- 1, 3 * k /**/- 1)     = tau_z[k];
                                A_solve(3 * thread_length + k /**/- 1, 3 * k - 1 /**/- 1) = tau_y[k];
                                A_solve(3 * thread_length + k /**/- 1, 3 * k - 2 /**/- 1) = tau_x[k];

                                A_solve(3 * thread_length + k /**/- 1, 3 * (k + 1) /**/- 1)     = -tau_z[k];
                                A_solve(3 * thread_length + k /**/- 1, 3 * (k + 1) - 1 /**/- 1) = -tau_y[k];
                                A_solve(3 * thread_length + k /**/- 1, 3 * (k + 1) - 2 /**/- 1) = -tau_x[k];
                                k = k + 2; // !!! previous particle related k
                            }

                            k++; // loop
                            tp_prev = tp;
                        }

                        //printf("\n DEBUG 0");
                        //printf("\n DEBUG ! step = %d thread_length = %d, threads number = %d", step, thread_length, number_of_mag_threads);
                        //det_A = Determinant(A_solve, n);
                        //printf("\n DEBUG ! det_A = %e %d", det_A, thread_length);
                        //printf("\n DEBUG ! A_solve[n][n] = %e %d", A_solve[n][n], thread_length);
                        //printf("\n DEBUG ! %e %d", A_solve[2][0], thread_length);

                        x_solve = A_solve.colPivHouseholderQr().solve(b_solve);

                        for (pos = 1; pos <= thread_length; pos++)
                        {
                            tp = mag_thread_id_ordered_p[thread_id][pos];

                            // System Solution - only final force
                            F[tp].x = M0 * x_solve(3 * pos - 2/**/- 1);
                            F[tp].y = M0 * x_solve(3 * pos - 1/**/- 1);
                            F[tp].z = M0 * x_solve(3 * pos/**/- 1);

                            if (pos > 1)
                            {
                                k = pos;
                                tau_x[k - 1] = r[tp].x - r[tp_prev].x;
                                tau_y[k - 1] = r[tp].y - r[tp_prev].y;
                                tau_z[k - 1] = r[tp].z - r[tp_prev].z;

                                tau_mod = sqrt(tau_x[k - 1] * tau_x[k - 1] + tau_y[k - 1] * tau_y[k - 1] + tau_z[k - 1] * tau_z[k - 1]);

                                //printf("\n DEBUG ! tau_mod = %e", tau_mod);

                                tau_x[k - 1] /= tau_mod;
                                tau_y[k - 1] /= tau_mod;
                                tau_z[k - 1] /= tau_mod;

                                F_tau_1 = (F[tp_prev].x * tau_x[pos - 1] + F[tp_prev].y * tau_y[pos - 1] + F[tp_prev].z * tau_z[pos - 1]);
                                F_tau_2 = (F[tp].x * tau_x[pos - 1] + F[tp].y * tau_y[pos - 1] + F[tp].z * tau_z[pos - 1]);

                                F[tp].x -= (F_tau_2 - F_tau_1) * tau_x[pos - 1];
                                F[tp].y -= (F_tau_2 - F_tau_1) * tau_y[pos - 1];
                                F[tp].z -= (F_tau_2 - F_tau_1) * tau_z[pos - 1];
                            }

                            //printf("\n DEBUG TEMP! thread_id = %d max = %d", thread_id, number_of_mag_threads);

                            tp_prev = tp;
                        }

                    } // ((thread_length >= 2) && (mag_thread_id_ordered_p[thread_id][0] != -1))
                } 
            } // if at least one mag-thread exists
}

/*void ff_mag_threads_align_v(void)
{
long tp, tp_prev, thread_length, pos, thread_id;

if (number_of_mag_threads)
for (thread_id = 1; thread_id <= number_of_mag_threads; thread_id++)
{
thread_length = mag_thread_limit_plus[thread_id] - mag_thread_limit_minus[thread_id] + 1;
if ((thread_length >= 2) && (mag_thread_id_ordered_p[thread_id][0] != -1)) // -1 is removed thread after threads connection
for (pos = 1; pos <= thread_length; pos++)
{
tp = mag_thread_id_ordered_p[thread_id][pos];

//TODO

tp_prev = tp;
}
}
}*/

void ff_model_next_step(void)
{ 

    ff_vect_t f;
    double I;
    long p;
    int chk;
    long mz_tot_n;
    ff_vect_t r0;

    Ek = 0;
    mz_tot = 0;
    m_tot.x = m_tot.y = m_tot.z = 0;
    mz_tot_n = 0;

    if (time_go)
    {
        printf("\n ============================", step);
        printf("\n Step %d", step);

        ff_model_m_setting();

        for (p = 1; p <= pN; p++)
            if (exist_p[p])
            {
                f = ff_model_force(p);
                F[p] = f;

                if (f.x != f.x) printf("\n DEBUG 1 p = %d f.x = %e", p, f.x);
                if (f.y != f.y) printf("\n DEBUG 1 p = %d f.y = %e", p, f.y);
                if (f.z != f.z) printf("\n DEBUG 1 p = %d f.z = %e", p, f.z);

                if (v[p].x != v[p].x) printf("\n DEBUG 1.1 p = %d r[p].x = %e v[p].x = %e", p, r[p].x, v[p].x);
                if (v[p].y != v[p].y) printf("\n DEBUG 1.1 p = %d r[p].y = %e v[p].y = %e", p, r[p].y, v[p].y);
                if (v[p].z != v[p].z) printf("\n DEBUG 1.1 p = %d r[p].z = %e v[p].z = %e", p, r[p].z, v[p].z);
            }

            if (mag_threads)
            {
                ff_mag_threads_update_force();
                //ff_mag_threads_align_v();
            }

            for (p = 1; p <= pN; p++)
                if (exist_p[p])
                {
                    drt[p].x = F[p].x * dt / C2 +		 
                        (v[p].x - F[p].x / C2) * (1 - exp(- C2 * dt / M0)) * M0 / C2;

                    drt[p].y = F[p].y * dt / C2 +		 
                        (v[p].y - F[p].y / C2) * (1 - exp(- C2 * dt / M0)) * M0 / C2;

                    drt[p].z = F[p].z * dt / C2 +		 
                        (v[p].z - F[p].z / C2) * (1 - exp(- C2 * dt / M0)) * M0 / C2;

                    if (brownian_shifts) ff_model_set_rand_dir(p);

                    r[p].x += drt[p].x;
                    r[p].y += drt[p].y;
                    r[p].z += drt[p].z;

                    if (v[p].x != v[p].x) printf("\n DEBUG 1.9 p = %d r[p].x = %e v[p].x = %e", p, r[p].x, v[p].x);
                    if (v[p].y != v[p].y) printf("\n DEBUG 1.9 p = %d r[p].y = %e v[p].y = %e", p, r[p].y, v[p].y);
                    if (v[p].z != v[p].z) printf("\n DEBUG 1.9 p = %d r[p].z = %e v[p].z = %e", p, r[p].z, v[p].z);

                    ff_model_check_walls(p);

                    if (r[p].x != r[p].x) printf("\n DEBUG 2 p = %d r[p].x = %e v[p].x = %e", p, r[p].x, v[p].x);
                    if (r[p].y != r[p].y) printf("\n DEBUG 2 p = %d r[p].y = %e v[p].y = %e", p, r[p].y, v[p].y);
                    if (r[p].z != r[p].z) printf("\n DEBUG 2 p = %d r[p].z = %e v[p].z = %e", p, r[p].z, v[p].z);

                    Ek += M0 * MUL(v[p],v[p]) / 2.0;

                    chk = ff_model_check_smooth_dr(p);
                    if ( chk == 0) goto t_end;
                } // end of loop for dr

                //ff_model_check_overlapp();

                r0.x = r0.y = r0.z = 0;

                for (p = 1; p <= pN; p++)
                    if (exist_p[p])
                    {
                        // C2 is a friction
                        dvt[p].x = (F[p].x / C2) + (v[p].x - F[p].x / C2) * exp(- C2 * dt / M0) - v[p].x;
                        dvt[p].y = (F[p].y / C2) + (v[p].y - F[p].y / C2) * exp(- C2 * dt / M0) - v[p].y;
                        dvt[p].z = (F[p].z / C2) + (v[p].z - F[p].z / C2) * exp(- C2 * dt / M0) - v[p].z;

                        v[p].x += dvt[p].x;
                        v[p].y += dvt[p].y;
                        v[p].z += dvt[p].z;

                        if (v[p].x != v[p].x) printf("\n DEBUG 3 p = %d v[p].x = %e", p, v[p].x);
                        if (v[p].y != v[p].y) printf("\n DEBUG 3 p = %d v[p].y = %e", p, v[p].y);
                        if (v[p].z != v[p].z) printf("\n DEBUG 3 p = %d v[p].z = %e", p, v[p].z);

                        //chk = ff_model_check_smooth_dv(p);
                        //if ( chk == 0) goto t_end;

                        r0.x += r[p].x;
                        r0.y += r[p].y;
                        r0.z += r[p].z;

                        mz_tot += m[p].z;
                        m_tot.x += m[p].x;
                        m_tot.y += m[p].y;
                        m_tot.z += m[p].z;
                        mz_tot_n++;
                    } // end of loop for dv

                    printf("\n DEBUG 5 r.x = %e v.x = %e", r[50].x, v[50].x);

                    // TODO: need the number of existing (exist_p) particles ! 
                    r0.x /= pN;
                    r0.y /= pN;
                    r0.z /= pN;

                    I = 0;
                    for (p = 1; p <= pN; p++)
                        I += pow(r[p].x - r0.x, 2)
                        + pow(r[p].y - r0.y, 2)
                        + pow(r[p].z - r0.z, 2);


                    mz_tot *= pN / mz_tot_n; // in loop is interrupted then needs to increase
                    mz_glob += mz_tot;

                    m_tot_glob.x += m_tot.x;
                    m_tot_glob.y += m_tot.y;
                    m_tot_glob.z += m_tot.z;

                    t += dt;

                    if (auto_reversal) ff_model_auto_hyst();
                    if ((auto_save)&&(step%10000 == 0)) ff_io_autosave();

                    if (setting_plot)
                    {
                        if (step < 1000)
                        {
                            ff_io_save_setting(m_tot,I);
                        }
                        else if (step%1000 == 0) ff_io_save_setting(m_tot,I);
                    }

                    // end for case of interrupted kinetic loops
t_end:

                    if (slow_steps > 0) slow_steps--;
                    if (slow_steps%10 == 1) dt *= 2;

    } // time_go

    ff_mgr_show_next_step();
}

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

ff_vect_t Bext(double x, double y, double z)
{
    ff_vect_t tBext;
    double normX, normY, normZ;

    if (!manual_field_control)
    {
        tBext.x = tBext.y = 1E-50;
        tBext.z = B0 * sin(kB * pi / 2.0);
    }
    else
    {
        if (ext_field_is_homo)
        {
            tBext.x = BmanX /*Oe*/ * 79.577 * mu0; // Tesla
            tBext.y = BmanY /*Oe*/ * 79.577 * mu0; // Tesla
            tBext.z = BmanZ /*Oe*/ * 79.577 * mu0; // Tesla
        }
        else
        {
            normX = (BmanX /*Oe*/ * 79.577 * mu0 / B_het(0, 0, 0, 1).x);
            normY = (BmanY /*Oe*/ * 79.577 * mu0 / B_het(0, 0, 0, 2).y);
            normZ = (BmanZ /*Oe*/ * 79.577 * mu0 / B_het(0, 0, 0, 3).z);

            tBext.x = normX * B_het(x, y, z, 1).x + normY * B_het(x, y, z, 2).x + normZ * B_het(x, y, z, 3).x;
            tBext.y = normX * B_het(x, y, z, 1).y + normY * B_het(x, y, z, 2).y + normZ * B_het(x, y, z, 3).y;
            tBext.z = normX * B_het(x, y, z, 1).z + normY * B_het(x, y, z, 2).z + normZ * B_het(x, y, z, 3).z;
        }
    }

    return tBext;
}

/*void ff_model_init_sediment(void)
{
long i, j, k;
ff_vect_t r0, r1, rb; //rb - basic; we build sphere around it
long n;
double a = (2 + gap) * R0;
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

r1.x = rb.x + a * dir110[jdir].x / sqrt(2);
r1.y = rb.y + a * dir110[jdir].y / sqrt(2);
r1.z = rb.z + a * dir110[jdir].z / sqrt(2);

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
}*/

/*void ff_model_dir_init(void)
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

void ff_model_init(void)
{
    long p, tp;
    double theta, phi;
    ff_vect_t dr;
    double dR;
    double sigma;
    long i,j;

#ifdef SECONDARY

    R00 = R0 = R00_agg;
    m0 = m0_agg; // Magnetic moment [J / T]

    Vself = Vself_agg; // [m^3]
    M0 = M0_agg;  // mass [kg]

    C2 = C2_agg;
    C3 = C3_agg;

    D = D_agg;
    r0mod = r0mod_agg; // needs extra * dt^(1/2.0)

    C6 = C6_agg;

#endif

    // Brownian motion - random force parameters
    ////////////////////////////////////////////
    sigma = sqrt(18 * kb * T * pi * eta * 2 * R0);

    srand( (unsigned)time( NULL ) );
    rng.seed(static_cast<unsigned int>(std::time(0)));

    //boost::normal_distribution<> nd(0.0, 1.0);
    nd = new boost::normal_distribution<> (0.0, sigma);

    //boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, *nd);
    var_nor = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > (rng, *nd);

    //ff_model_dir_init();

    t = 0; // time

    //printf("\n %e", m0);
    //printf("\n %e", M0);
    glob_start_t = start_t;

    m_tot_glob.x = m_tot_glob.y = m_tot_glob.z = 0;

    for (long h = 1; h <= 20; h++)
    {
        B_hyst[h] = B_hyst_n[h] = Mz_hyst[h] = Mz_hyst_n[h] = 0;
    }

    p = 1;
    while (p <= pN)
    {
again:

        if (start_ideal)
        {
            r[p].x = -0.4999 * Lx + 0.99 * Lx * rand() / 32768.0;
            r[p].y = -0.4999 * Ly + 0.99 * Ly * rand() / 32768.0;
            r[p].z = -0.4999 * Lz + 0.99 * Lz * rand() / 32768.0;

            //if (p == 5)
            //r[p].x = r[p].y = r[p].z = 0;

            for (tp = 1; tp < p; tp++)
            {
                dr.x = r[p].x - r[tp].x;
                dr.y = r[p].y - r[tp].y;
                dr.z = r[p].z - r[tp].z;

                dR = sqrt(MUL(dr,dr));

                if (dR <= 4 * R0) goto again;
            }
        } // start_ideal

        // 4 means that <|v|> must be ~ heat speed
        v[p].x =  4 * sqrt(8 * kb * 300 / (pi * M0)) * (rand() / 32768.0 - 0.5);
        v[p].y =  4 * sqrt(8 * kb * 300 / (pi * M0)) * (rand() / 32768.0 - 0.5);
        v[p].z =  4 * sqrt(8 * kb * 300 / (pi * M0)) * (rand() / 32768.0 - 0.5);

        theta = pi * rand() / 32768.0;
        phi = 2 * pi * rand() / 32768.0;

        m[p].x = m0 * sin(theta) * cos(phi);
        m[p].y = m0 * sin(theta) * sin(phi);
        m[p].z = m0 * cos(theta);


#ifdef SECONDARY
        mag_thread_p_id[p] = 0;
        mag_thread_p_position[p] = 0;
        mag_thread_limit_plus[p] = 0;
        mag_thread_limit_minus[p] = 0;

        for (tp = 1; tp <= pN; tp++) mag_thread_id_ordered_p[p][tp] = 0;

        m[p].x = m0start_agg * sin(theta) * cos(phi);
        m[p].y = m0start_agg * sin(theta) * sin(phi);
        m[p].z = m0start_agg * cos(theta);

        Rp[p] = R0; // TODO: replace in all logic !
        exist_p[p] = 1;
        m0p[p] = m0; // TODO: replace in all logic !

        m0start_agg_p[p] = m0start_agg;
        r0modp[p] = r0mod;

        Vselfp[p] = Vself;

        /*if (p == 5)
        {
        Rp[p] = 5 * R0;
        exist_p[p] = 125;
        m0p[p] = 125 * m0;

        m0start_agg_p[p] = 125 * m0start_agg;
        r0modp[p] = r0mod * pow(125, - 1 / 6.0);

        Vselfp[p] = 125 * Vself;
        }*/

#endif
        //m_freeze[p] = 0;
        m_sat[p] = 0;

        p++;
    }

    if (load_at_start) ff_io_load();

    //	if (start_sediment) ff_model_init_sediment();



}
