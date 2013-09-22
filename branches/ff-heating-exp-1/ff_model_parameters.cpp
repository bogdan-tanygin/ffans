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
double R = 8.31;
double Na = 6.02214179 * 1E23;
double g = 9.81;
double kb = 1.38 * (1E-23);

// Space

double scale = 0.4;
double Lx = 60 * 1E-6 * scale, Ly = 30 * 1E-6 * scale, Lz = 5 * 1E-6 * scale; //meters

//double kExtra = 0.27 * 20 * 1; // change only this coef. instead Lz
//double Lx = 10 * kExtra * 1E-6, Ly = kExtra * 1E-6, Lz = kExtra * 1E-6;

// Basic physical model parameters
double dt = 1E-1; //s
long slow_steps;
//double smooth_v = 10; // disabled in code
double smooth_r = 0.1;
//double m_h_eff_tol = 1; // max. angle [rad] between m and B

double R00 = 10 * (1E-9); // Radius of the nanoparticle [m]
double R0 = 12 * (1E-9); // Radius including the acid sphere [m]
double Vself = (4 * pi / 3.0) * pow(R00, 3); // [m^3]

double rop = 0.5 * (4.9 + 5.2) * (1E+3); // mass density [kg / m^3]
double M0 = Vself * rop;  // mass [kg]

double Ms_mass = 80 /* emu / g */ * (1E3) /* emu / kg */ * (1 / (9.274009 * (1E-21))); /* Bohr magnetons / kg */
double m0 = Ms_mass * M0 /* Bohr magnetons */* 927.400915 * (1E-26); // Magnetic moment [J / T]

double Ch = 1E-8; // adhesion / magnetic relation

int load_at_start = 0; // the best method of loading
int auto_save = 1;
int manual_field_control = 1; // 0-1-2-3 keys control to skip, Bx+, By+, Bz+ control 
int ext_field_is_homo = 1;
int auto_reversal = 0;

int setting_plot = 0; // cluster creation plot m_tot / m0 and I.

double start_t =30 /* [micro_s] */ * (1E-6); // [s]
double T_ext = 0.2 * (1E6); // [micro_s] // external field period
double nu_ext = (1 / T_ext) * (1E6); // [Hz] // frequency of the external field (sin(w*t) dependence)

double start_ideal = 1; // start chaos (ideal superparam. gas)
double ro0 = 0.5 * (0.78 + 0.85) * 1E3; // kerosene density
double eta = 0.00164; //Pa * s //kerosene

int brownian_shifts = 0;
int brownian_force = 1;

double T = 300; // K

//default order of magnitude of the external field but exact function is hardcoded 
double B0 = 100 /*Oe*/ * 79.577 * mu0; // Tesla

// Derived parameters
double C1 = 3 * mu0 / (4 * pi);
double C2 = 6 * pi * eta * R0;
double C3 = M0 * g;
double D = R * T / (6 * Na * pi * R0 * eta);
//double r0 = sqrt(3 * 2 * D * dt);
double r0mod = sqrt(3 * 2 * D); // needs extra * dt^(1/2.0)
//double C4 = r0;
double C5 = mu0 / (4 * pi);

//double Mself = m0 / Vself;

//double Hself = - (1 / 3.0) * Mself;
//double Bself = mu0 * (Hself + Mself);

// 0.740 is atomic packing factor of the fcc lattice
// total mean total vol. of the clusters
double Vtot = (4 / 3.0) * pi * pow(R0, 3) * (1 / 0.740) * pN;

double C6 = ro0 * Vself * g;

//double gamma_e = 1.7609 * (1E11); //s^-1 T^-1
//double alpha   = 0.05; //magnetization dynamic damping