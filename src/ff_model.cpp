/**************************************************************************
* Copyright (C) 2011,2013-2014,2017 Dr. Bogdan Tanygin<b.m.tanygin@gmail.com>
* Copyright (C) 2016,2017 Dmytro Matskevych<dimqqqq@mail.ru>
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

#include "ff_model.h"
#include "ff_model_parameters.h"
#include "ff_model_graphics.h"
#include "ff_model_io.h"
#include "ff_analysis.h"
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
ff_vect_t tau[pN + 1]; // mean torque
ff_vect_t tau1[pN + 1]; // torque at the t = t
ff_vect_t tau2[pN + 1]; // // torque at the t = t + dt

ff_vect_t v[pN + 1];
ff_vect_t w[pN + 1]; // angular velocity vector

ff_vect_t drt[pN + 1];
ff_vect_t drt_r[pN + 1]; // instantiated random translation
ff_vect_t dvt[pN + 1];
ff_vect_t dvt_r[pN + 1];
ff_vect_t dphi[pN + 1];
ff_vect_t dphi_r[pN + 1]; // instantiated random rotation
ff_vect_t dw[pN + 1];
ff_vect_t dw_r[pN + 1];

long i_min = 1;
double V0_tot = 0; // total volume of the dispersed phase
double V0_largest_EV = 0; // mathematical expected value of largest particles total volume // see is_large_mode variable
double V0_tot_EV = 0; // mathematical expected value of particles total volume
double I_glob = 0;
double phi_vol_fract = 0;

int exist_p[pN + 1]; // particle existence; number of primary aggregate inside
int is_neel[pN + 1]; // Neel relaxation
int is_temp_sat[pN + 1]; // temperature saturation flag
double Rp0[pN + 1];
double Rp[pN + 1];
double Vp0[pN + 1];
double Vpfull[pN + 1]; // including steric layer
double m0p[pN + 1];
double M0p[pN + 1];
double I0p[pN + 1]; // particle moment of inertia
double r0modp[pN + 1];
double C2[pN + 1];
double gamma_rot[pN + 1];

int time_go = 0;

double dt = 0;
double t = 0; // time
double t_temp = 0; // temperature, [C]
double T_basic = 0;
double T_mean = 0;
long k_mean = 0;

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
double g_Bz_prev;
long step = 0;

double kB = 0;
int hyst_mode = 1;
double mz_tot;
double mz_glob = 0; // global mean average start from the hyst. point switch
ff_vect_t m_tot_glob;

double BmanX = 0;
double BmanY = 0;
double BmanZ = 0;

ff_vect_t m_tot;

ff_vect_t dir110[13];

ff_vect_t r_brown_valid_0;

double d[14 + 1];

double dt_red = 0; // reducing time indicator for the random translation

long pN0 = pN;
double R0_min;

int ff_model_check_smooth_dr(long p)
{
    double dr, rmod, dphimag;
    long ps;
    int res = 1;
    
    dr = sqrt(MUL(drt[p], drt[p]));
    rmod = 2 * Rp[p];
    dphimag = sqrt(MUL(dphi[p], dphi[p]));

    if (rmod > 0)
        if ((dr / rmod > smooth_r) || ((dphimag / pi > smooth_r) && (!(is_neel[p]))))
        { 
            for(ps = 1; ps <= p; ps ++)
            {
                r[ps].x -= drt[ps].x;
                r[ps].y -= drt[ps].y;
                r[ps].z -= drt[ps].z;

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
        }

        return res;
}

ff_vect_t ff_model_nonloc_torque(long p)
{
    register long ps;
    double mxs, mys, mzs;
    double dx, dy, dz;

    double MUL2mod__dR5mod1, dR2__dR5mod1;
    double dR, dR2, dR5, dR5mod1, MUL2mod;

    double tBx, tBy, tBz;
    double dtBx, dtBy, dtBz;    // field of single ps-th particle
    double tBmag;

	ff_vect_t tBext;
    ff_vect_t ttau;

    ttau.x = ttau.y = ttau.z = 0;

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
                    {
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

                        if (dR >= Rp0[p] + Rp0[ps]) // soft-sphere model related correction
                        {
                            tBx += dtBx;
                            tBy += dtBy;
                            tBz += dtBz;
                        }

                    } // di and dj
                } // if (p != ps)

                tBmag = sqrt(tBx * tBx + tBy * tBy + tBz * tBz);

#ifndef SECONDARY

                if (tBmag)
				{
					ttau.x =   m[p].y * tBz - m[p].z * tBy;
					ttau.y = - m[p].x * tBz + m[p].z * tBx;
					ttau.z =   m[p].x * tBy - m[p].y * tBx;
				}

                if ((is_neel[p]) && (tBmag))
                {
					m[p].x = m0p[p] * tBx / tBmag;
                    m[p].y = m0p[p] * tBy / tBmag;
                    m[p].z = m0p[p] * tBz / tBmag;
                }

				G_dd[p] = - (m[p].x * (tBx - tBext.x) + m[p].y * (tBy - tBext.y) + m[p].z * (tBz - tBext.z));
#endif

    }
    return ttau;
}

ff_vect_t ff_model_nonloc_force(long p)
{
    register long ps;
    double tFx, tFy, tFz;
    double dtFx, dtFy, dtFz;
    double dr, dr2, dr3, dr2mod, dr5, dr5mod, MUL1, MUL2, MUL3;
    double Cmod;
    double dx, dy, dz;
    ff_vect_t ttF, dBmaggrad;

    double MUL1__dr5mod, MUL2__dr5mod, MUL3__dr5mod, MUL1_MUL2__dr2mod__dr5mod;

    double mx, my, mz;
    double mxs, mys, mzs;
	
    double dr_min, F_steric_mag;

    tFx = 0;
    tFy = 0;
    tFz = 0;

    for(ps = 1; ps <= pN; ps ++)
        if (exist_p[ps])
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

                dr2mod = dr2 / 5.0; //modified
                dr5mod = (dr5 = pow(dr,5)) / C1; //modified

                MUL2 = mxs * dx + mys * dy + mzs * dz;

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

                if (dr >= dr_min)
                {
                    tFx += dtFx;
                    tFy += dtFy;
                    tFz += dtFz;
                }

                if (dr < dr_min) //soft sphere condition
                {
                    Cmod = 10 * m0p[p] * m0p[ps] * (C1 / dr5);
					
                    if (dr)
					{
						tFx += - dx * Cmod;
						tFy += - dy * Cmod;
						tFz += - dz * Cmod;
					}
                }

                // Entropic repulsion
                if ((dr >= dr_min) && (dr <= Rp0[p] + Rp0[ps] + 2 * delta)) 
                {
					F_steric_mag = ((-pi)*kb*T*N_o*(dr-2*delta-Rp0[ps]-Rp0[p])*((Rp0[ps]+Rp0[p])*(2*(dr*dr+delta*delta)+(-pow(Rp0[ps],2)+Rp0[p]*Rp0[ps]-pow(Rp0[p],2))*(dr/(Rp0[ps]+Rp0[p])+1)+delta*dr)+(2*Rp0[p]*Rp0[ps]-pow((Rp0[p]-Rp0[ps]),2))*delta))/(6*delta*dr*dr);

                    tFx += - (dx / dr) * F_steric_mag;
                    tFy += - (dy / dr) * F_steric_mag;
                    tFz += - (dz / dr) * F_steric_mag;
                }

                if (dr >= dr_min)
                {

                    tFx += - (dx / dr) * (- 64 * A_H * dr * pow(Rp0[p] * Rp0[ps], 3) / (6 * pow(dr * dr - pow(Rp0[p] - Rp0[ps], 2), 2) * pow(dr * dr - pow(Rp0[p] + Rp0[ps], 2), 2)));
                    tFy += - (dy / dr) * (- 64 * A_H * dr * pow(Rp0[p] * Rp0[ps], 3) / (6 * pow(dr * dr - pow(Rp0[p] - Rp0[ps], 2), 2) * pow(dr * dr - pow(Rp0[p] + Rp0[ps], 2), 2)));
                    tFz += - (dz / dr) * (- 64 * A_H * dr * pow(Rp0[p] * Rp0[ps], 3) / (6 * pow(dr * dr - pow(Rp0[p] - Rp0[ps], 2), 2) * pow(dr * dr - pow(Rp0[p] + Rp0[ps], 2), 2)));
                }
            } // end of the particles loop

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
    }

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
    double M0, I0;
    double tmmag;
    double Mtot;
    double t_temp_1 = 0;

	double per_x_plus, per_x_minus, per_y_plus, per_y_minus, per_z_plus, per_z_minus;

    Ek = Ek_tr = Ek_rot = 0;
    mz_tot = 0;
    m_tot.x = m_tot.y = m_tot.z = 0;
    mz_tot_n = 0;

    t_temp = T + ta0;
    if (t_temp > 90) t_temp_1 = 90;
    if (t_temp < 20) t_temp_1 = 20;
    ro_oleic = 902.0 - 0.62 * T;
    eta_oleic = a3_eta_oleic * pow(t_temp_1, 3) + a2_eta_oleic * pow(t_temp_1, 2) + a1_eta_oleic * pow(t_temp_1, 1) + a0_eta_oleic;
	
    //Screen capturing
    //----------------
	if(step==1000000)
	{
		BmanZ=2000;
		ScreenCaptureStep = 2000;
	}
	if(v_oleic!=0)
	{
	eta_car = ff_visousity_mix(
								ff_molar_part(
									ff_mol(mass_oleic,mol_mass_oleic),
									ff_mol(mass_car,mol_mass_car)),
								eta_oleic,
								ff_molar_part(
									ff_mol(mass_car,mol_mass_car),
									ff_mol(mass_oleic,mol_mass_oleic)),
									eta_car0);
	}
	else
	{
	eta_car = eta_car0;//Oleic acid = 0 mll
	}
	//cout<<eta_car<<" "<<eta_oleic<<" "<<eta_car0<<endl;
	
	//if(step%ScreenCaptureStep==0){counterOfPosition=0;}
	if(step%ScreenCaptureStep == counterOfPosition*10 && step>=ScreenCaptureStep)
	{
		ChangePosition();
		ostringstream out;
		out<<"step ="<<step << " V_oleic = " << v_oleic<<" V_car = "<<v_car<<" Bmanz ="<<BmanZ <<".bmp";
		GetScreenShot(out.str());
		out.flush();
	}
    //----------------

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
        ff_model_update_dT();

        for (p = 1; p <= pN; p++) if (exist_p[p]) ff_model_effective_random_motion_update(p);
        if (dt_red <= 0)
        {
            dt_red = dt0;
        }

        // START: Velocity Verlet based Analytical Dissipative Integrator.
        // STEP 0. Forces for t=0 (beginning of dt interval), x(0), v(0). 
        // Translational and rotational motion of bead particles.
        for (p = 1; p <= pN; p++)
        if (exist_p[p])
        {
            f = ff_model_force(p);
            ttau = ff_model_torque(p);
            F1[p].x = f.x; F1[p].y = f.y; F1[p].z = f.z;
            tau1[p].x = ttau.x; tau1[p].y = ttau.y; tau1[p].z = ttau.z;
        }

        // STEP 1. Velocity v(0.5*dt) calculation: deterministic part.
           for (p = 1; p <= pN; p++)
                    if (exist_p[p])
                    {
                        C2[p] = 6 * pi * eta_car * Rp[p];
                        gamma_rot[p] = 8 * pi * eta_car * pow(Rp[p], 3);


                        M0 = M0p[p];
                        I0 = I0p[p];

                        // C2 is a friction
                        dvt[p].x = (F[p].x / C2[p]) + (v[p].x - F[p].x / C2[p]) * exp(- C2[p] * 0.5 * dt / M0) - v[p].x;
                        dvt[p].y = (F[p].y / C2[p]) + (v[p].y - F[p].y / C2[p]) * exp(- C2[p] * 0.5 * dt / M0) - v[p].y;
                        dvt[p].z = (F[p].z / C2[p]) + (v[p].z - F[p].z / C2[p]) * exp(- C2[p] * 0.5 * dt / M0) - v[p].z;

                        v[p].x += dvt[p].x;
                        v[p].y += dvt[p].y;
                        v[p].z += dvt[p].z;

                        w[p].x = (tau[p].x / gamma_rot[p]) + (w[p].x - tau[p].x / gamma_rot[p]) * exp(- gamma_rot[p] * 0.5 * dt / I0);
                        w[p].y = (tau[p].y / gamma_rot[p]) + (w[p].y - tau[p].y / gamma_rot[p]) * exp(- gamma_rot[p] * 0.5 * dt / I0);
                        w[p].z = (tau[p].z / gamma_rot[p]) + (w[p].z - tau[p].z / gamma_rot[p]) * exp(- gamma_rot[p] * 0.5 * dt / I0);

                        if (v[p].x != v[p].x) printf("\n DEBUG 3 p = %d v[p].x = %e", p, v[p].x);
                        if (v[p].y != v[p].y) printf("\n DEBUG 3 p = %d v[p].y = %e", p, v[p].y);
                        if (v[p].z != v[p].z) printf("\n DEBUG 3 p = %d v[p].z = %e", p, v[p].z);
                    } // end of loop for dv
        
        // STEP 2. Position x(dt) calculation: stochastic and deterministic parts.
        // Random walk
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

                m_prev_before_r[p] = m[p];
                m[p] = mt[p];

                tmmag = sqrt(MUL(m[p], m[p]));

                m[p].x *= m0p[p] / tmmag;
                m[p].y *= m0p[p] / tmmag;
                m[p].z *= m0p[p] / tmmag;
            }
        }  // end of loop for dr_r
        
        // Deterministic position update
            for (p = 1; p <= pN; p++)
                if (exist_p[p])
                {
                        C2[p] = 6 * pi * eta_car * Rp[p];
                        gamma_rot[p] = 8 * pi * eta_car * pow(Rp[p], 3);

                    M0 = M0p[p];
                    I0 = I0p[p];

                    drt[p].x = F1[p].x * dt / C2[p] +		 
                        (v[p].x - F1[p].x / C2[p]) * (1 - exp(- C2[p] * dt / M0)) * M0 / C2[p];

                    drt[p].y = F1[p].y * dt / C2[p] +		 
                        (v[p].y - F1[p].y / C2[p]) * (1 - exp(- C2[p] * dt / M0)) * M0 / C2[p];
                    
                    drt[p].z = F1[p].z * dt / C2[p] +		 
                        (v[p].z - F1[p].z / C2[p]) * (1 - exp(- C2[p] * dt / M0)) * M0 / C2[p];
                    
                    if (!(is_neel[p]))
                    {
                        dphi[p].x = tau1[p].x * dt / gamma_rot[p] +		 
                            (w[p].x - tau1[p].x / gamma_rot[p]) * (1 - exp(- gamma_rot[p] * dt / I0)) * I0 / gamma_rot[p];
                        
                        dphi[p].y = tau1[p].y * dt / gamma_rot[p] +		 
                            (w[p].y - tau1[p].y / gamma_rot[p]) * (1 - exp(- gamma_rot[p] * dt / I0)) * I0 / gamma_rot[p];
                        
                        dphi[p].z = tau1[p].z * dt / gamma_rot[p] +		 
                            (w[p].z - tau1[p].z / gamma_rot[p]) * (1 - exp(- gamma_rot[p] * dt / I0)) * I0 / gamma_rot[p];
                    }
                    else dphi[p].x = dphi[p].y = dphi[p].z = 0;

                    /*DEBUG*/ if (dphi[p].x != dphi[p].x) printf("\n DEBUG 1 p = %d dphi[p].x = %e", p, dphi[p].x);
                    /*DEBUG*/ if (dphi[p].y != dphi[p].y) printf("\n DEBUG 1 p = %d dphi[p].y = %e", p, dphi[p].y);
                    /*DEBUG*/ if (dphi[p].z != dphi[p].z) printf("\n DEBUG 1 p = %d dphi[p].z = %e", p, dphi[p].z);

                    r[p].x += drt[p].x;
                    r[p].y += drt[p].y;
                    r[p].z += drt[p].z;

                    if (r[p].x != r[p].x) printf("\n DEBUG 2 p = %d r[p].x = %e v[p].x = %e", p, r[p].x, v[p].x);
                    if (r[p].y != r[p].y) printf("\n DEBUG 2 p = %d r[p].y = %e v[p].y = %e", p, r[p].y, v[p].y);
                    if (r[p].z != r[p].z) printf("\n DEBUG 2 p = %d r[p].z = %e v[p].z = %e", p, r[p].z, v[p].z);

                    chk = ff_model_check_smooth_dr(p);
                    if ( chk == 0)
                    {
                        goto t_end;
                    }

                    if (!(is_neel[p]))
                    {
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

                        m_prev[p] = m[p];
                        m[p] = mt[p];

                        tmmag = sqrt(MUL(m[p], m[p]));

                        m[p].x *= m0p[p] / tmmag;
                        m[p].y *= m0p[p] / tmmag;
                        m[p].z *= m0p[p] / tmmag;
                    }

                    ff_model_check_walls(p);
                } // end of loop for dr
        
        // STEP 3. Forces f(dt) calculation.
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

            /*DEBUG*/ if (v[p].x != v[p].x) printf("\n DEBUG 1.1 p = %d r[p].x = %e v[p].x = %e", p, r[p].x, v[p].x);
            /*DEBUG*/ if (v[p].y != v[p].y) printf("\n DEBUG 1.1 p = %d r[p].y = %e v[p].y = %e", p, r[p].y, v[p].y);
            /*DEBUG*/ if (v[p].z != v[p].z) printf("\n DEBUG 1.1 p = %d r[p].z = %e v[p].z = %e", p, r[p].z, v[p].z);

            /*DEBUG*/ if (w[p].x != w[p].x) printf("\n DEBUG 1.2 p = %d", p);
            /*DEBUG*/ if (w[p].y != w[p].y) printf("\n DEBUG 1.2 p = %d", p);
            /*DEBUG*/ if (w[p].z != w[p].z) printf("\n DEBUG 1.2 p = %d", p);
        } // force calc
        
           // STEP 4. Velocity v(dt) calculation: stochastic (dt step) and deterministic (dt/2 step) parts.
           for (p = 1; p <= pN; p++)
                    if (exist_p[p])
                    {
                        C2[p] = 6 * pi * eta_car * Rp[p];
                        gamma_rot[p] = 8 * pi * eta_car * pow(Rp[p], 3);

                        M0 = M0p[p];
                        I0 = I0p[p];

                        // C2 is a friction
                        dvt[p].x = (F2[p].x / C2[p]) + (v[p].x - F2[p].x / C2[p]) * exp(- C2[p] * 0.5 * dt / M0) - v[p].x;
                        dvt[p].y = (F2[p].y / C2[p]) + (v[p].y - F2[p].y / C2[p]) * exp(- C2[p] * 0.5 * dt / M0) - v[p].y;
                        dvt[p].z = (F2[p].z / C2[p]) + (v[p].z - F2[p].z / C2[p]) * exp(- C2[p] * 0.5 * dt / M0) - v[p].z;

                        v[p].x += dvt[p].x + dvt_r[p].x;
                        v[p].y += dvt[p].y + dvt_r[p].y;
                        v[p].z += dvt[p].z + dvt_r[p].z;

                        w[p].x = (tau2[p].x / gamma_rot[p]) + (w[p].x - tau2[p].x / gamma_rot[p]) * exp(- gamma_rot[p] * 0.5 * dt / I0);
                        w[p].y = (tau2[p].y / gamma_rot[p]) + (w[p].y - tau2[p].y / gamma_rot[p]) * exp(- gamma_rot[p] * 0.5 * dt / I0);
                        w[p].z = (tau2[p].z / gamma_rot[p]) + (w[p].z - tau2[p].z / gamma_rot[p]) * exp(- gamma_rot[p] * 0.5 * dt / I0);

						w[p].x += dw_r[p].x;
						w[p].y += dw_r[p].y;
						w[p].z += dw_r[p].z;

                        if (v[p].x != v[p].x) printf("\n DEBUG 3 p = %d v[p].x = %e", p, v[p].x);
                        if (v[p].y != v[p].y) printf("\n DEBUG 3 p = %d v[p].y = %e", p, v[p].y);
                        if (v[p].z != v[p].z) printf("\n DEBUG 3 p = %d v[p].z = %e", p, v[p].z);
                    } // end of loop for dv

                r0.x = r0.y = r0.z = 0;
				for (int t1 = 0; t1 <= 2; t1++)
					for (int t2 = 0; t2 <= 2; t2++)
						for (int t3 = 0; t3 <= 2; t3++)
						{
							I_per[t1][t2][t3] = 0;
						}

				Mtot = 0;
				
                for (p = 1; p <= pN; p++)
                    if (exist_p[p])
                    {
                        M0 = M0p[p];
                        I0 = I0p[p];

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
						}

                        mz_tot += m[p].z;
                        m_tot.x += m[p].x;
                        m_tot.y += m[p].y;
                        m_tot.z += m[p].z;
                        mz_tot_n++;
                    } // end of loop for dv

                    // TODO: need the number of existing (exist_p) particles
                    r0.x /= Mtot;
                    r0.y /= Mtot;
                    r0.z /= Mtot;

                    I = 0;
					
                    for (p = 1; p <= pN; p++)
					{
                        I += M0p[p] * (pow(r[p].x - r0.x, 2)
                        + pow(r[p].y - r0.y, 2)
                        + pow(r[p].z - r0.z, 2));
					}

					if (is_periodic)
					{
						I_min_per = I_per[0][0][0];
						for (int t1 = 0; t1 <= 2; t1++) for (int t2 = 0; t2 <= 2; t2++) for (int t3 = 0; t3 <= 2; t3++)
							if (I_per[t1][t2][t3] < I_min_per) I_min_per = I_per[t1][t2][t3];
					}

					mz_tot *= pN / mz_tot_n; // in loop is interrupted then needs to increase
                    mz_glob += mz_tot;

                    m_tot_glob.x += m_tot.x;
                    m_tot_glob.y += m_tot.y;
                    m_tot_glob.z += m_tot.z;

                    t += dt;
                    dt_red -= dt;

                    for (p = 1; p <= pN; p++)
                    {
                        Rp_to_c[p] = sqrt(MUL(r[p], r[p]));
                    }

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
    } // time_go

    ff_mgr_show_next_step();
}

int ff_model_check_walls(long p)
{
    int res = 0;

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
			// [Archived]
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
    long p, tp, p_prev;
    double theta, phi;
    ff_vect_t dr;
    double dR;
    double sigma;
    double t_temp_1 = 0;
    double dr_min;
    int cont_flag;
    long i_attempt; 
    double start_scale = 0.99;

    dt = dt0;

    t_temp = T + ta0;
    if (t_temp > 90) t_temp_1 = 90;
    if (t_temp < 20) t_temp_1 = 20;
    eta_oleic = a3_eta_oleic * pow(t_temp_1, 3) + a2_eta_oleic * pow(t_temp_1, 2) + a1_eta_oleic * pow(t_temp_1, 1) + a0_eta_oleic;

    N_o = N_oa * k_o / 0.5;

    // Brownian motion -  parameters
    ///////////////////////////////////////////////////

    // Dimensionless variance (sigma^2) of the random displacement along the single axis e_x
    sigma = 1;

    srand( (unsigned)time( NULL ) );
    rng.seed(static_cast<unsigned int>(time(0)));

    nd = new boost::normal_distribution<> (0.0, sigma);
    ud = new boost::uniform_01<> ();

    var_nor = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > (rng, *nd);
    var_uni = new boost::variate_generator<boost::mt19937&, boost::uniform_01<> > (rng, *ud);

    t = 0; // time

    m_tot_glob.x = m_tot_glob.y = m_tot_glob.z = 0;

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

            w[p].x = w[p].y = w[p].z = 0;

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
        w[p].x = w[p].y = w[p].z = 0;
        drt_r[p].x = drt_r[p].y = drt_r[p].z = 0;
        dphi_r[p].x = dphi_r[p].y = dphi_r[p].z = 0;
		dvt_r[p].x = dvt_r[p].y = dvt_r[p].z = 0;
		dw_r[p].x = dw_r[p].y = dw_r[p].z = 0;
        exist_p[p] = 1;
        m_sat[p] = 0;

        is_temp_sat[p] = 0;

		G_dd[p] = 0;

        p++;
    }

    if (load_at_start) ff_io_load(0);

    r_brown_valid_0 = r[1];
}

// Update of the random motion
void ff_model_effective_random_motion_update(long p)
{
    double theta_0, phi_0; // random direction of the torque 
    double dx, dy, dz, dvx, dvy, dvz, dphi, dw; // instantiated displacements for time dt * k_bm_inst_max
    double D, D_rot, gamma, gamma_rot;
    double gamma_oleic = 0, gamma_car = 0; // oleic vs. carrier liquid
    double M0, I0;

    double sigma2 = 0, sigma2_rot = 0;
	double sigma2v = 0, sigma2w_rot = 0;

    gamma_oleic = 6 * pi * eta_oleic * Rp[p];
    gamma_car = 6 * pi * eta_car * Rp[p];

    gamma = gamma_car;
    gamma_rot = 8 * pi * eta_car * pow(Rp[p], 3);


    D = kb * T / gamma;
    D_rot = kb * T / gamma_rot;
    M0 = M0p[p];
    I0 = I0p[p];

    // instantiation of position change
    // [Langevin equation + Stokes' law]
    sigma2 = D * (2 * dt + (M0 / (2 * gamma)) * (- 3 + 4 * exp(- gamma * dt / M0) - exp(- 2 * gamma * dt / M0)));
    dx = (*var_nor)() * sqrt(sigma2);
    dy = (*var_nor)() * sqrt(sigma2);
    dz = (*var_nor)() * sqrt(sigma2);

    drt_r[p].x = dx;
    drt_r[p].y = dy;
    drt_r[p].z = dz;

	sigma2v = (kb * T / M0) * (1 - exp(-2 * gamma * dt / M0));
	dvx = (*var_nor)() * sqrt(sigma2v);
	dvy = (*var_nor)() * sqrt(sigma2v);
	dvz = (*var_nor)() * sqrt(sigma2v);

	dvt_r[p].x = dvx;
	dvt_r[p].y = dvy;
	dvt_r[p].z = dvz;

    // instantiation of rotation of the magnetization direction
    // [Euler-Langevin equation + Stokes' law]
	sigma2_rot = D_rot * (2 * dt + (I0 / (2 * gamma_rot)) * (- 3 + 4 * exp(- gamma_rot * dt / I0) - exp(- 2 * gamma_rot * dt / I0)));
    dphi = (*var_nor)() * sqrt(3 * sigma2_rot); // rotation magnitude
    theta_0 = (*var_uni)() * pi;   // rotation vector random direction
    phi_0 = (*var_uni)() * 2 * pi;

    dphi_r[p].x = dphi * sin(theta_0) * cos(phi_0);
    dphi_r[p].y = dphi * sin(theta_0) * sin(phi_0);
    dphi_r[p].z = dphi * cos(theta_0);

	sigma2w_rot = (kb * T / I0) * (1 - exp(-2 * gamma_rot * dt / I0));
	dw = (*var_nor)() * sqrt(3 * sigma2w_rot); // rotation angular velocity change magnitude
	theta_0 = (*var_uni)() * pi;   // rotation vector random direction
	phi_0 = (*var_uni)() * 2 * pi;

	dw_r[p].x = dw * sin(theta_0) * cos(phi_0);
	dw_r[p].y = dw * sin(theta_0) * sin(phi_0);
	dw_r[p].z = dw * cos(theta_0);
}

void ff_model_update_dT(void)
{
    T_basic = (2 / 6.0) * Ek / (kb * pN); // degree of freedom number is 6
    T_mean += T_basic;
    k_mean ++;
}

void ff_model_size_dispersion_init(void)
{
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
	
    V0_tot = 0;
    for (p = 1; p <= pN; p++)
    {
        is_set = 0;

        again_size_disp:
        random_value = (*var_uni)();
		
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
		phi_vol_fract = 100 * V0_tot / (Lx * Ly * Lz);
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

// [Archived]
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

// Let's keep this function archived. It can be treated as a hard-sphere condition.
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

// [Archived]
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