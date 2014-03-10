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

// SYSTEM SI

#include <math.h>

#define pN 50 // number of the particles

// Math constants
extern double pi;

// Physics constants
extern const double mu0;
extern const double muB;
extern const double R;
extern const double Na;
extern double g;
extern const double kb;
extern const double ta0;

extern double C1;

// Geometry
extern double Lx, Ly, Lz; //meters

// Basic physical model parameters
extern double dt;
extern long k_bm_inst_max;
extern long k_bm_inst;
extern double k_force_adapt_0;
extern long slow_steps;
//extern double smooth_v;
extern double smooth_r;
//extern double m_h_eff_tol;

extern double Ch;
extern double Ch_ss;
//extern double EPS;
//extern double ro1;
//extern double ro2;
extern int load_at_start;
extern int auto_reversal;
extern int auto_save;
extern int manual_field_control;
extern int ext_field_is_homo;
extern int setting_plot;
extern double start_t; // [s]
extern double nu_ext;

extern double start_ideal;
extern double start_sediment;

//extern int __deprecated__brownian_shifts;
//extern int __deprecated__brownian_force;

//extern double gap;

//extern double R00; // Radius of the nanoparticle [m]
extern double delta;
//extern double R0; // Radius including the acid sphere [m]
//extern double M0; // mass [kg]
//extern double m0; // Magnetic moment [J / T]

extern int is_oleic;
extern double R_oleic_0;
extern double eta_oleic;
extern double a3_eta_oleic;
extern double a2_eta_oleic;
extern double a1_eta_oleic;
extern double a0_eta_oleic;
extern double sigma_sf;
extern double a_sigma_sf;
extern double b_sigma_sf;
extern double sigma_sf_nano;

extern double T;
extern double kr;
extern double rop;
extern double A_H;
extern double N_oa;
extern double a0;
extern double eta_car;
extern double G_barrier;

//order of magnitude of the external field but exact function is hardcoded 
extern double B0; // Tesla

// Derived parameters
//extern double C2;
//extern double C3;
//extern double D;
//extern double gamma;
//extern double r0;
//extern double r0mod;
//extern double C4;
extern double C5;

extern double Vself;
//extern double Mself;

//extern double Hself;
//extern double Bself;

//extern double Vtot;

//extern double C6;

extern double gamma_e;
extern double alpha;