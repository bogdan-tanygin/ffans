/**************************************************************************
 * Copyright (C) 2011 Dr. Bogdan Tanygin<b.m.tanygin@gmail.com>
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

#include "mf_model.h"
#include "mf_model_parameters.h"
#include "mf_model_graphics.h"
#include "mf_io_pal.h"

// working variables 
////////////////////
mf_vect_t r[pN + 1];  //particles positions
mf_vect_t m[pN + 1];  //particles magnetic moment direction

mf_vect_t F[pN + 1];

mf_vect_t v[pN + 1];

mf_vect_t drt[pN + 1];

mf_vect_t dvt[pN + 1];

mf_vect_t Fnonloc[pN + 1][pN + 1];

mf_vect_t dir110[13];

int time_go = 1;

long double B_hyst[21];
long double Mz_hyst[21];
long double B_hyst_n[21];
long double Mz_hyst_n[21];

long double t; // time
long double Ek;
long double g_Bz_prev;
long step = 0;

long double kB = 0;
int hyst_mode = 1;
long double mz_tot;
long glob_start_step = 1;
long double glob_start_t;
long double mz_glob = 0; // global mean average start from the hyst. point switch 

long double g_theta[pN + 1], g_phi[pN + 1];

int g_hyst_start_line;
int g_hyst_up_line;
int g_hyst_bottom_line;

void mf_model_upgrade_ext_field(void)
{
    g_Bz_prev = Bext().z;
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

void mf_model_auto_hyst(void)
{
        int tmp;

        if ((glob_start_step == 1) && (t - glob_start_t >= 0))
            glob_start_step = step;


        if ((t - glob_start_t) >= 0.25 /*0.25 of the 1/4.0 of period*/ * 0.25 * (1 / nu_ext))
        {

        hyst_mode++;
        tmp = 0;
        if (hyst_mode == 21) {hyst_mode = 5; tmp = 1;}
        mf_model_upgrade_ext_field();
        if (tmp) {g_hyst_up_line = 1;g_hyst_start_line = 0;}

        mf_io_save_hyst();

        mz_glob = 0;
        glob_start_t = t;
        glob_start_step = step;

        }
}

int mf_model_check_smooth_dr(long p)
{
    long double dr, rmod;
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
int mf_model_check_smooth_dv(long p)
{
    long double dv, vmod, dm, mmag;
    long ps;
    int res = 1;
    mf_vect_t v_prev;

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

/*long double mf_model_distrib_f(long double theta, long p)
{
    long double ret_f;
    mf_vect_t B1; // excluding self field
    long double Bmag;
    long double Ct; // temp constant
    long double dtheta = pi / 10.0; // hardcode

    B1.x = B[p].x - Bself * m[p].x / m0;
    B1.y = B[p].y - Bself * m[p].y / m0;
    B1.z = B[p].z - Bself * m[p].z / m0;

    Bmag = sqrt(MUL(B1,B1));

    Ct = m0 * Bmag / (kb * T);

    ret_f = (Ct / (exp(Ct) - exp(-Ct))) * exp(Ct * cos(theta)) * sin(theta) * dtheta;

    return ret_f;
}*/

long double mf_model_distrib_f_r(long double theta, long p)
{
    long double ret_f;
    long double G0;
    long double Ct; // temp constant
    long double dtheta = pi / 10.0; // hardcode

    G0 = sqrt(MUL(F[p],F[p])) * r0mod * sqrt(dt);
    Ct = G0 / (kb * T);

    ret_f = (Ct / (exp(Ct) - exp(-Ct))) * exp(Ct * cos(theta)) * sin(theta) * dtheta;

    return ret_f;
}

void mf_model_set_rand_dir(long p)
{
    long i;
    long double dtheta; // step in [R]andom long
    long double R_points[10 + 1];
    long R;
    long double theta, phi;
    long double Fmag, tmag, bmag;
    mf_vect_t Ft, t, b;

    dtheta = pi / 10.0;
    for (i = 1; i <= 10; i++)
    {
        R_points[i] = 32768.0 * mf_model_distrib_f_r((i - 0.5) * dtheta, p);
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

    drt[p].x += r0mod * sqrt(dt) * (Ft.x * cos(theta) / Fmag + t.x * sin(theta) * cos(phi) / tmag + b.x * sin(theta) * sin(phi) / bmag);
    drt[p].y += r0mod * sqrt(dt) * (Ft.y * cos(theta) / Fmag + t.y * sin(theta) * cos(phi) / tmag + b.y * sin(theta) * sin(phi) / bmag);
    drt[p].z += r0mod * sqrt(dt) * (Ft.z * cos(theta) / Fmag + t.z * sin(theta) * cos(phi) / tmag + b.z * sin(theta) * sin(phi) / bmag);
}

/*void mf_model_set_m(long p)
{
    long i;
    long double dtheta; // step in [R]andom long
    long double R_points[10 + 1];
    long R;
    long double theta, phi;
    long double Bmag, tmag, bmag;
    mf_vect_t Bt, t, b;

    dtheta = pi / 10.0;
    for (i = 1; i <= 10; i++)
    {
        R_points[i] = 32768.0 * mf_model_distrib_f((i - 0.5) * dtheta, p);
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

void mf_model_m_setting(void)
{
    long p;
    register long ps;
    long double mxs, mys, mzs;
    long double dx, dy, dz;
    long double eps, max_eps;

    long double MUL2mod__dR5mod1, dR2__dR5mod1;
    long double dR, dR2, dR5, dR5mod1, MUL2mod;

    long double tBx, tBy, tBz;
    long double tBmag;
    mf_vect_t ttB;

    do
    {
    for(p = 1; p <= pN; p ++)
    {
        ttB = Bext();
        tBx = ttB.x;
        tBy = ttB.y;
        tBz = ttB.z;

        for(ps = 1; ps <= pN; ps ++)
        if (p != ps)
        {
           dx = r[ps].x - r[p].x;
           dy = r[ps].y - r[p].y;
           dz = r[ps].z - r[p].z;

           mxs = m[ps].x;
           mys = m[ps].y;
           mzs = m[ps].z;

           MUL2mod = 3 * (mxs * dx + mys * dy + mzs * dz);

           dR = sqrt(dR2 = dx * dx + dy * dy + dz * dz);
           dR5 = pow(dR,5);
           dR5mod1 = dR5 / C5;


           MUL2mod__dR5mod1 = MUL2mod / dR5mod1;
           dR2__dR5mod1 = dR2 / dR5mod1;

           tBx += dx * MUL2mod__dR5mod1 - mxs * dR2__dR5mod1;
           tBy += dy * MUL2mod__dR5mod1 - mys * dR2__dR5mod1;
           tBz += dz * MUL2mod__dR5mod1 - mzs * dR2__dR5mod1;

        } // if (p != ps)

        tBmag = sqrt(tBx * tBx + tBy * tBy + tBz * tBz);

        eps = acos((m[p].x * tBx + m[p].y * tBy + m[p].z * tBz) / (tBmag * m0));

        if (p == 1) max_eps = eps;
        else if (eps > max_eps) max_eps = eps;

        m[p].x = m0 * tBx / tBmag;
        m[p].y = m0 * tBy / tBmag;
        m[p].z = m0 * tBz / tBmag;

    }// for(p = 1; p <= pN; p ++)
    }
    while(max_eps > m_h_eff_tol);
}

mf_vect_t mf_model_nonloc_force(long p)
{
    register long ps;
    long double tFx, tFy, tFz;
    long double dR, dR2, dR2mod, dR5, dR5mod, MUL1, MUL2, MUL3;
    long double R1 = 0.3 * R0;
    long double Cmod;
    long double Cmod1 = Ch * m0 * m0;
    long double Cmod1_repulsion = 5 * m0 * m0;
    long double dx, dy, dz;
    mf_vect_t ttF;

    long double MUL1__dR5mod, MUL2__dR5mod, MUL3__dR5mod, MUL1_MUL2__dR2mod__dR5mod;

    long double mx, my, mz;
    long double mxs, mys, mzs;
    
    tFx = 0;
    tFy = 0;
    tFz = 0;
    
    for(ps = 1; ps <= pN; ps ++)
    if (p != ps)
      {

        dx = r[ps].x - r[p].x;
        dy = r[ps].y - r[p].y;
        dz = r[ps].z - r[p].z;

        mx = m[p].x;
        my = m[p].y;
        mz = m[p].z;

        mxs = m[ps].x;
        mys = m[ps].y;
        mzs = m[ps].z;
        
        //dipole-dipole
        dR = sqrt(dR2 = dx * dx + dy * dy + dz * dz);

        // Non local force and self-confid. magnetic field
        dR2mod = dR2 / 5.0; //modified
        dR5mod = (dR5 = pow(dR,5)) / C1; //modified
                
        MUL2 = mxs * dx + mys * dy + mzs * dz;
        
        if (p < ps)
        {
        MUL1 = mx * dx + my * dy + mz * dz;
        MUL3 = mx * mxs + my * mys + mz * mzs;

        MUL1__dR5mod = MUL1 / dR5mod;
        MUL2__dR5mod = MUL2 / dR5mod;
        MUL3__dR5mod = MUL3 / dR5mod;
        MUL1_MUL2__dR2mod__dR5mod = MUL1 * MUL2 / (dR2mod * dR5mod);

        tFx -= MUL1__dR5mod * mxs + MUL2__dR5mod * mx
                + MUL3__dR5mod * dx - MUL1_MUL2__dR2mod__dR5mod * dx;

        tFy -= MUL1__dR5mod * mys + MUL2__dR5mod * my
                + MUL3__dR5mod * dy - MUL1_MUL2__dR2mod__dR5mod * dy;

        tFz -= MUL1__dR5mod * mzs + MUL2__dR5mod * mz
                + MUL3__dR5mod * dz - MUL1_MUL2__dR2mod__dR5mod * dz;

/*
        if (dR <= 5 * R0)
        {
            tFx += - (EPS * dx / dR) * (exp(-(dR - 2 * R00) / ro1) / ro1 - exp(-(dR - 2 * R00) / ro2) / ro2);
            tFy += - (EPS * dy / dR) * (exp(-(dR - 2 * R00) / ro1) / ro1 - exp(-(dR - 2 * R00) / ro2) / ro2);
            tFz += - (EPS * dz / dR) * (exp(-(dR - 2 * R00) / ro1) / ro1 - exp(-(dR - 2 * R00) / ro2) / ro2);
        }
*/
        // acid elasticity (repulsion)
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

        // attraction
        if ((dR > 2 * R0 )&&(dR < 3 * R0 )) // the Heaviside step function  and dR5 dependence finally is similar to the well-known exp. phenomenology
        {
            Cmod = Cmod1 * (C1 / dR5);

            tFx += dx * Cmod;
            tFy += dy * Cmod;
            tFz += dz * Cmod;
        }
        Fnonloc[p][ps].x = tFx;
        Fnonloc[p][ps].y = tFy;
        Fnonloc[p][ps].z = tFz;
        } // if not p < ps
        else
        {
        tFx -= Fnonloc[ps][p].x; // third Newton law
        tFy -= Fnonloc[ps][p].y;
        tFz -= Fnonloc[ps][p].z;
        }

        //tBx += dx * MUL2mod__dR5mod1 - mxs * dR2__dR5mod1;
        //tBy += dy * MUL2mod__dR5mod1 - mys * dR2__dR5mod1;
        //tBz += dz * MUL2mod__dR5mod1 - mzs * dR2__dR5mod1;
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
    //mf_model_set_m(p);

    // force
    ttF.x = tFx;
    ttF.y = tFy;
    ttF.z = tFz;
  
    return ttF;
}

mf_vect_t mf_model_force(long p)
{
    mf_vect_t tF;
    mf_vect_t tddF;
    
    tF.x = tF.y = tF.z = 0;
    
    // non-local
    tddF = mf_model_nonloc_force(p);
    tF.x += tddF.x;
    tF.y += tddF.y;
    tF.z += tddF.z;
    
    // Gravitation
    tF.z += - C3;

    // Buoyancy force
    tF.z +=   C6;
    
    return tF;
}

int mf_model_check_walls(long p)
{
     int res = 0;
    
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

void mf_model_next_step(void)
{ 

    mf_vect_t f;
    long double I;
    long p;
    int chk;
    long mz_tot_n;
    mf_vect_t m_tot, r0;
        
    Ek = 0;
    mz_tot = 0;
    m_tot.x = m_tot.y = m_tot.z = 0;
    mz_tot_n = 0;

if (time_go)
{
    mf_model_m_setting();

    for (p = 1; p <= pN; p++)
    {
     f = mf_model_force(p);
     F[p] = f;

     drt[p].x = f.x * dt / C2 +         
                (v[p].x - f.x / C2) * (1 - exp(- C2 * dt / M0)) * M0 / C2;

     drt[p].y = f.y * dt / C2 +         
                (v[p].y - f.y / C2) * (1 - exp(- C2 * dt / M0)) * M0 / C2;

     drt[p].z = f.z * dt / C2 +         
                (v[p].z - f.z / C2) * (1 - exp(- C2 * dt / M0)) * M0 / C2;

     // Fluctuations are important but should be not here. They are in the r section. r^2 = D * dt in the random direction;
     // This is model approximation;
     // the molecular forces which produce fluct. are stronger and stronger than other forces. So this one time jump or decrease of the energy.

     if (brownian_enable) mf_model_set_rand_dir(p);

     r[p].x += drt[p].x;
     r[p].y += drt[p].y;
     r[p].z += drt[p].z;
         
     mf_model_check_walls(p);
     
     Ek += M0 * MUL(v[p],v[p]) / 2.0;

     chk = mf_model_check_smooth_dr(p);
     if ( chk == 0) goto t_end;
    } // end of loop for dr
    



    r0.x = r0.y = r0.z = 0;

    for (p = 1; p <= pN; p++)
    {
     // C2 is a friction
     dvt[p].x = (F[p].x / C2) + (v[p].x - F[p].x / C2) * exp(- C2 * dt / M0) - v[p].x;
     dvt[p].y = (F[p].y / C2) + (v[p].y - F[p].y / C2) * exp(- C2 * dt / M0) - v[p].y;
     dvt[p].z = (F[p].z / C2) + (v[p].z - F[p].z / C2) * exp(- C2 * dt / M0) - v[p].z;

     v[p].x += dvt[p].x;
     v[p].y += dvt[p].y;
     v[p].z += dvt[p].z;
     
     //chk = mf_model_check_smooth_dv(p);
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

    t += dt;

    if (auto_reversal) mf_model_auto_hyst();
    if ((auto_save)&&(step%1000 == 0)) mf_io_autosave();
    
    if (setting_plot)
    {
        if (step < 1000)
        {
            mf_io_save_setting(m_tot,I);
        }
        else if (step%1000 == 0) mf_io_save_setting(m_tot,I);
    }

    // end for case of interrupted kinetic loops
    t_end:

    if (slow_steps > 0) slow_steps--;
    if (slow_steps%10 == 1) dt *= 2;

} // time_go

    mf_mgr_show_next_step();
}

mf_vect_t Bext(void)
{
    mf_vect_t tBext;
    tBext.x = tBext.y = 1E-50;
    tBext.z = B0 * sin(kB * pi / 2.0);

    return tBext;
}

/*void mf_model_init_sediment(void)
{
    long i, j, k;
    mf_vect_t r0, r1, rb; //rb - basic; we build sphere around it
    long n;
    long double a = (2 + gap) * R0;
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

/*void mf_model_dir_init(void)
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

void mf_model_init(void)
{
    long p, tp;
    long double theta, phi;
    mf_vect_t dr;
    long double dR;

    //    mf_model_dir_init();

    t = 0; // time
    
    printf("\n %e", m0);
    printf("\n %e", M0);
    glob_start_t = start_t;

    for (long h = 1; h <= 20; h++)
    {
        B_hyst[h] = B_hyst_n[h] = Mz_hyst[h] = Mz_hyst_n[h] = 0;
    }

    srand( (unsigned)time( NULL ) );

    p = 1;
    while (p <= pN)
    {
        again:

        if (start_ideal)
        {
        r[p].x = -0.4999 * Lx + 0.99 * Lx * rand() / 32768.0;
        r[p].y = -0.4999 * Ly + 0.99 * Ly * rand() / 32768.0;
        r[p].z = -0.4999 * Lz + 0.99 * Lz * rand() / 32768.0;

        for (tp = 1; tp < p; tp++)
        {
            dr.x = r[p].x - r[tp].x;
            dr.y = r[p].y - r[tp].y;
            dr.z = r[p].z - r[tp].z;

            dR = sqrt(MUL(dr,dr));

            if (dR <= 2.2 * R0) goto again;
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

        p++;
    }

//    if (start_sediment) mf_model_init_sediment();
}
