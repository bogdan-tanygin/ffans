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

// SYSTEM SI

#include <math.h>
#include "ff_model_parameters.h"

// Math constants
double pi = acos(-1.0);

// Physics constants
double mu0 = 4 * pi * 1E-7;
const double muB = 9.27400968 * 1E-24; // Bohr magneton
double R = 8.31;
double Na = 6.02214179 * 1E23;
double g = 9.81;
const double kb = 1.3806488 * 1E-23; // [m2 kg s-2 K-1]

// Space
double gl_scale = 1;
double volume_reduce = 0.45;
double scale = 2 * 7 * 0.15 * volume_reduce * gl_scale / pow(500.0 / pN, 1 / 3.0); // / pow(50.0, 1 / 3.0);
double Lx = 1E-6 * scale, Ly = 1E-6 * scale, Lz = 1E-6 * scale; //meters

//double kExtra = 0.27 * 20 * 1; // change only this coef. instead Lz
//double Lx = 10 * kExtra * 1E-6, Ly = kExtra * 1E-6, Lz = kExtra * 1E-6;

// Basic physical model parameters
double dt = 1.5E-5; // s
long k_bm_inst_max = 100; // coefficient of a brownian motion instantiation: dt_bm_inst = dt * k_bm_inst_max
long k_bm_inst = 1;
double k_force_adapt_0 = 1.05;
long slow_steps = 0;
//double smooth_v = 10; // disabled in code
double smooth_r = 0.2;
//double m_h_eff_tol = 1; // max. angle [rad] between m and B

double T = 273.15 + 50; // K
double kr = gl_scale; // particle size parameter []
//double R00 = 0.5 * 15E-9; // Radius of the nanoparticle [m]
double delta = 2.0E-9;
//double R0 = R00 + delta; // Radius including the acid sphere [m]
//double Vself = (4 * pi / 3.0) * pow(R00, 3); // [m^3]

// Parameters of oleic acid drop
int is_oleic = 0;
double R_oleic_0 = (Lx / 8.0);
double eta_oleic = 25.6 * 1E-3; // [Pa * s]
double sigma_sf = 32.5 * 1E-3; // [N / m]
//double sigma_sf_nano = sigma_sf * 1 * 5.017559E-04; // [N / m]
//double sigma_sf_nano = 5E-2 * sigma_sf;
double sigma_sf_nano = 1 * sigma_sf;

double rop = 5240; // mass density [kg / m^3]
//double M0 = Vself * rop;  // mass [kg]
double A_H = 1E-19; // Hamaker constant [J]
double N_oa = 10E19; // Surface density of oleic acid at the 50% coating [m-2]
//double G_barrier = pow(kr, 2) * 25 * kb * T; // [TEMP] barrier which prevents particles aggregation

double Ms_mass = 80 /* emu / g */ * (1E3) /* emu / kg */ * (1 / (9.274009 * (1E-21))); /* Bohr magnetons / kg */
//double m0 = Ms_mass * M0 /* Bohr magnetons */* 927.400915 * (1E-26); // Magnetic moment [J / T]

double Ch = 0.01; // [DEPRECATED] adhesion / magnetic relation
double Ch_ss = 1 * 1E5; // soft-sphere repulsion parameter (parameter of the numerical model)

int load_at_start = 0;
int auto_save = 1;
int manual_field_control = 1; // 0-1-2-3 keys control to skip, Bx+, By+, Bz+ control 
int ext_field_is_homo = 1;
int auto_reversal = 0;

int setting_plot = 1; // cluster creation plot m_tot / m0 and I.

double start_t =30 /* [micro_s] */ * (1E-6); // [s]
double T_ext = 0.2 * (1E6); // [micro_s] // external field period
double nu_ext = (1 / T_ext) * (1E6); // [Hz] // frequency of the external field (sin(w*t) dependence)

double start_ideal = 1; // start chaos (ideal superparam. gas)
double start_sediment = 0;
double ro0 = 0.5 * (0.78 + 0.85) * 1E3; // kerosene density
double eta = 0.00164; //Pa * s //kerosene

//int __deprecated__brownian_shifts = 0;
//int __deprecated__brownian_force = 1;

//default order of magnitude of the external field but exact function is hardcoded 
double B0 = 100 /*Oe*/ * 79.577 * mu0; // Tesla

// Derived parameters
double C1 = 3 * mu0 / (4 * pi);
//double C2 = 6 * pi * eta * R0;
//double C3 = M0 * g;
//double D = R * T / (6 * Na * pi * R0 * eta);
//double gamma = 6 * pi * R0 * eta;
//double r0 = sqrt(3 * 2 * D * dt);
//double r0mod = sqrt(3 * 2 * D); // needs extra * dt^(1/2.0)
//double C4 = r0;
double C5 = mu0 / (4 * pi);

//double Mself = m0 / Vself;

//double Hself = - (1 / 3.0) * Mself;
//double Bself = mu0 * (Hself + Mself);

// 0.740 is atomic packing factor of the fcc lattice
// total mean total vol. of the clusters
//double Vtot = (4 / 3.0) * pi * pow(R0, 3) * (1 / 0.740) * pN;

//double C6 = ro0 * Vself * g;

//double gamma_e = 1.7609 * (1E11); //s^-1 T^-1
//double alpha   = 0.05; //magnetization dynamic damping