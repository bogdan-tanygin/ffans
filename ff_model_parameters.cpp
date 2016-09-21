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
#include "ff_model_parameters.h"

double gl_scale = 1.0; // particle diameter scale parameter

// Math constants
double pi = acos(-1.0);

// Physics constants
const double mu0 = 4 * pi * 1E-7;
const double muB = 9.27400968 * 1E-24; // Bohr magneton
const double R = 8.31;
const double Na = 6.02214179 * 1E23;
double g = 9.81;
const double kb = 1.3806488 * 1E-23; // [m2 kg s-2 K-1]
const double ta0 = -273.15; // [°C]
const double gamma_e = 1.760859708 * 1E11; // [s^-1 T^-1] // electron gyromagnetic ratio

// Material
double a0 = 0.8397E-9; // [m] // magnetite unit cell size - a cubic spinel structure with space group Fd3m (above the Verwey temperature) // \cite{Fleet1981}

// Space
double volume_reduce = 0.15; //0.25; // 0.09; //0.085; // 0.1862; // initial density set
double scale = 2 * 7 * 0.15 * volume_reduce * gl_scale / pow(500.0 / pN, 1 / 3.0); // / pow(50.0, 1 / 3.0);
//double Lx = 1E-6 * scale, Ly = 1E-6 * scale, Lz = 25E-6 * scale; //meters
double nano_size = 950;
double Lx = nano_size * 1E-9 / 3.0, Ly = nano_size * 1E-9, Lz = nano_size * 1E-9 * 3; //meters
double delta_r = a0 * 0.5; // minimal distance between particles // order of magnitude of the oleic acid molecule width
double delta_r_init = 20 * 1E-9; // delta_r;

// periodic boundary conditions
int is_periodic = 0;

//double kExtra = 0.27 * 20 * 1; // change only this coef. instead Lz
//double Lx = 10 * kExtra * 1E-6, Ly = kExtra * 1E-6, Lz = kExtra * 1E-6;

// Basic physical model parameters
double dt_neel = 1E-9; // [s] // Neel relaxation time threashold for d ~ 10 nm, (Fertman-p-62)
const double dt0 = 1E-9; //1E-1 * dt_neel; // 1E2 * dt_neel; // 100 * dt_neel; // s // old approach 1.5E-5
double d_neel = 10 * 1E-9; // [m] // Neel-to-Brown relaxation diameter threashold, (Fertman-p-62)
double d_min = 2 * 1E-9; // [m] // minimal diameter of particle, wher Ms and T_curie is identical to ones in the bulk material, (Fertman-p-43)
long k_bm_inst_max = 100; // coefficient of a brownian motion instantiation: dt_bm_inst = dt * k_bm_inst_max
long k_bm_inst = 1;
double k_force_adapt_0 = 1.00; // 1.05 is a regular value for the adaptive force model // 1.00 means an overdamped model at the dt >> m / gamma
long slow_steps = 0;
//double smooth_v = 10; // disabled in code
double smooth_r = 0.4;
//double m_h_eff_tol = 1; // max. angle [rad] between m and B

double T = 273.15; // K
double sigma_sf_nano = 5E-4; //1E-4;

double kr = gl_scale; // particle size parameter []
//double R00 = 0.5 * 15E-9; // Radius of the nanoparticle [m]
double delta = 2.0E-9;
//double R0 = R00 + delta; // Radius including the acid sphere [m]
//double Vself = (4 * pi / 3.0) * pow(R00, 3); // [m^3]

int is_uniform_field_test = 0;

// Parameters of oleic acid drop
int is_large_mode = 1; // largest particles mode
double large_fraction = 7.5E-2; // 0.1 (Ivanov, Phase separation in bidisperse ferrocolloids) // 6.92E-02 (my sim); // fraction of the largest particles which form the primary aggregate (circle) in case of mode is_large_mode == 0
double k_large = 1.0; // 0.91; // 0.925; // correction of the large particles size
int is_oleic = 0;
double R_oleic_0 = (Lx / 2.0);

// microdrop mode parameters
int isMicroDrop = 0; // enable the microdrop mode
double phi_v = 0.12; // volume concentration of the dispersed phase [Padalka, exp]
double alpha = 0.5; // default saturation of the microdrop

int isPGasMode = 0;
int isPGasPrevail = 0;
double P_pgas = 0; // particles gas pressure inside the oleic drop
double P_sf_oleic = 0; // oleic droplet surface tension pressure
//double eta_oleic = 25.6 * 1E-3; // [Pa * s]
//double eta_oleic = 14.285 * 1E-3; // [Pa * s]
double eta_oleic = 0;
double a3_eta_oleic = - 1E-07; // my approximation of DOI: 10.1007/s11746-000-0197-z
double a2_eta_oleic = 2E-05;
double a1_eta_oleic = - 0.0018;
double a0_eta_oleic = 0.0559;
//double sigma_sf = 32.5 * 1E-3; // [N / m]
double sigma_sf = 0;
double a_sigma_sf = 34.060119 * 1E-3; // [N / m] // linear coefficient: sigma_sf = a + b * t
double b_sigma_sf = - 0.061298 * 1E-3; // [N / m] // linear coefficient
double ro_oleic = 0.853;// density of oleic acid g/sm^3
double mol_mass_oleic = 282; //mol mass of oleic acid C18H34O2
double v_oleic =  1* 1E-4; //volume of oleic acid
double mass_oleic = ro_oleic * v_oleic;

//double sigma_sf_nano = sigma_sf * 1 * 5.017559E-04; // [N / m]
//double sigma_sf_nano = 5E-2 * sigma_sf;

double rop = 5240; // magnetite mass density [kg / m^3]
//double M0 = Vself * rop;  // mass [kg]
double A_H = 4E-20; // Hamaker constant [J]

double N_oa = 1E18; // Surface density of oleic acid at the 50% coating [m-2] //[Fertman]
double k_o = 0.5; // 5E-4 - same result //0.1; // Level of coverage of surface by the oleic acid
double N_o; // Surface density of oleic acid

//double G_barrier = pow(kr, 2) * 25 * kb * T; // [TEMP] barrier which prevents particles aggregation
double K1 = 1.35 * 1E4 ; // [J/m] // First constant of the magnetite crystallographic anisotropy // (Goya2003)
double Ms_mass = 80 /* emu / g */ * (1E3) /* emu / kg */ * (1 / (9.274009 * (1E-21))); /* Bohr magnetons / kg */
//double m0 = Ms_mass * M0 /* Bohr magnetons */* 927.400915 * (1E-26); // Magnetic moment [J / T]
double Ms = 478 * 1E3; // [A / m]

double Ch = 0.01; // [DEPRECATED] adhesion / magnetic relation
double Ch_ss = 0; // 3E4; // 50 * 1E5; // soft-sphere repulsion parameter (parameter of the numerical model)

int load_at_start = 0;
int auto_save = 1;
int manual_field_control = 1; // 0-1-2-3 keys control to skip, Bx+, By+, Bz+ control 
int ext_field_is_homo = 1;
int auto_reversal = 0;

int setting_plot = 1; // cluster creation plot m_tot / m0 and I.

double start_t = 30 /* [micro_s] */ * (1E-6); // [s]
double T_ext = 0.2 * (1E6); // [micro_s] // external field period
double nu_ext = (1 / T_ext) * (1E6); // [Hz] // frequency of the external field (sin(w*t) dependence)

double start_ideal = 1; // start chaos (ideal superparam. gas)
double start_sediment = 0;
double ro0 = 0.5 * (0.78 + 0.85) * 1E3; // kerosene density
double eta_car = 0.00164; //Pa * s //kerosene (carrier liquid)
double mol_mass_car = 170; //molar mass of kerosene
double v_car = 3*1E-3; //volume of kerosine
double mass_car = ro0*v_car;//mas of kerosine
//int __deprecated__brownian_shifts = 0;
//int __deprecated__brownian_force = 1;

//default order of magnitude of the external field but exact function is hardcoded 
double B0 = 2500 /*Oe*/ * 79.577 * mu0; // Tesla
double gradPerc = 5E-1;
double gradL = Lx / 2.0;

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

double alpha_damp = 0.05; //magnetization dynamic damping