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
extern double gl_scale;
extern int isAutoSetPosition;
// Math constants
extern double pi;
extern int isShowInfo;
// Physics constants
extern const double mu0;
extern  double muB;
extern  double R;
extern  double Na;
extern  double g;
extern  double kb;
extern  double ta0;
extern  double gamma_e;

extern double C1;

// Geometry
extern double Lx, Ly, Lz; //meters
extern double delta_r;
extern double delta_r_init;
extern double nano_size;

extern int is_periodic;

// Basic physical model parameters
extern double dt0;
extern double d_neel;
extern long slow_steps;
extern double smooth_r;
extern double dt_neel;
extern double d_min;

extern int load_at_start;
extern int auto_reversal;
extern int auto_save;
extern int manual_field_control;
extern int ext_field_is_homo;
extern int setting_plot;

extern double start_ideal;

extern double delta;

extern int is_large_mode;
extern double large_fraction;
extern double k_large;
extern double eta_oleic;
extern double a3_eta_oleic;
extern double a2_eta_oleic;
extern double a1_eta_oleic;
extern double a0_eta_oleic;
extern double ro_oleic;
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
