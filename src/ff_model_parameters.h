/**************************************************************************
* Copyright (C) 2011,2013-2014,2017 Dr. Bogdan Tanygin<b.m.tanygin@gmail.com>
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

#define pN 500 // number of the particles

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
extern const double gamma_e;

extern double C1;

// Geometry
extern double Lx, Ly, Lz; //meters
extern double delta_r;
extern double delta_r_init;

extern int is_periodic;

// Basic physical model parameters
extern const double dt0;
extern double d_neel;
extern long k_bm_inst_max;
extern long k_bm_inst;
extern double k_force_adapt_0;
extern long slow_steps;
extern double smooth_r;

extern double Ch;
extern double Ch_ss;
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

extern double delta;

extern int is_uniform_field_test;

extern int isMicroDrop;
extern double phi_v;
extern double alpha;

extern int is_large_mode;
extern double large_fraction;
extern double k_large;
extern int is_oleic;
extern double R_oleic_0;
extern int isPGasMode;
extern int isPGasPrevail;
extern double P_pgas;
extern double P_sf_oleic;
extern double eta_oleic;
extern double a3_eta_oleic;
extern double a2_eta_oleic;
extern double a1_eta_oleic;
extern double a0_eta_oleic;
extern double sigma_sf;
extern double a_sigma_sf;
extern double b_sigma_sf;
extern double ro_oleic;
extern double sigma_sf_nano;
extern double mol_mass_oleic;
extern double v_oleic;
extern double mass_oleic; 

extern double T;
extern double kr;
extern double rop;
extern double A_H;
extern double N_oa;
extern double k_o;
extern double N_o;
extern double a0;
extern double ro0;
extern double eta_car;
extern double eta_car0;
extern double G_barrier;
extern double mol_mass_car;
extern double v_car;
extern double mass_car;


//order of magnitude of the external field but exact function is hardcoded 
extern double B0; // Tesla
extern double gradPerc;
extern double gradL;

// Derived parameters
extern double C5;

extern double Vself;

extern double alpha_damp;

extern double Ms;
extern double K1;

extern int ScreenCaptureStep; 

extern void ParamInfo();
