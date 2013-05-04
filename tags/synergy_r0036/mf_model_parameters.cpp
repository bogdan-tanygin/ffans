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

// SYSTEM SI

#include <math.h>
#include "mf_model_parameters.h"

// Math constants
long double pi = acos(-1);

// Physics constants
long double mu0 = 4 * pi * 1E-7;
long double R = 8.31;
long double Na = 6.02214179 * 1E23;
long double g = 9.81;
long double kb = 1.38 * (1E-23);

// Space
long double kExtra = 0.27; // please, change only this coef.
long double kL = 1 * kExtra;
long double Lx = 1 * kL * 1E-6, Ly = 1 * kL * 1E-6, Lz = kL * 1E-6; //meters

// Basic physical model parameters
long double dt = 1E-5; //s
long slow_steps;
//long double smooth_v = 10; // disabled in code
long double smooth_r = 0.2;
long double m_h_eff_tol = 0.2; // max. angle [rad] between m and B

long double R00 = 10 * (1E-9); // Radius of the nanoparticle [m]
long double R0 = 12 * (1E-9); // Radius including the acid sphere [m]
long double Vself = (4 * pi / 3.0) * pow(R00, 3); // [m^3]

long double rop = 0.5 * (4.9 + 5.2) * (1E-3) * (1E+6); // mass density [kg / m^3]
long double M0 = Vself * rop;  // mass [kg]

long double Ms_mass = 80 /* emu / g */ * (1E3) /* emu / kg */ * (1 / (9.274009 * (1E-21))); /* Bohr magnetons / kg */
long double m0 = Ms_mass * M0 /* Bohr magnetons */* 927.400915 * (1E-26); // Magnetic moment [J / T]

long double Ch = 0.25; // adhesion / magnetic relation
int auto_reversal = 0;
int auto_save = 1;
int setting_plot = 1; // cluster creation plot m_tot / m0 and I.
long double start_t =30 /* [micro_s] */ * (1E-6); // [s]
long double T_ext = 0.2 * (1E6); // [micro_s] // external field period
long double nu_ext = (1 / T_ext) * (1E6); // [Hz] // frequency of the external field (sin(w*t) dependence)

//long double EPS = 8 * (1E-3) /* eV */ * 1.6 * (1E-19); // attraction and repulsive steric layer forces [J]
//long double ro1 = 0.25 * (1E-9); // [m]
//long double ro2 = 0.5 * (1E-9); // [m]

long double start_ideal = 1; // start chaos (ideal superparam. gas)
//long double start_sediment = 0;

int brownian_enable = 1;

//long double gap = 4; // initial gap in R0 units; only for the start_sediment
// gap = 3 - connected protoclusters // needs a good classification

//long double gap = 0.0; // initial gap in R0 units

long double ro0 = 1 * 1E3; // water density

long double T = 293; // K
long double eta = 8.8 * (1E-4); //Pa * s

//order of magnitude of the external field but exact function is hardcoded 
long double B0 = 100 /*Oe*/ * 79.577 * mu0; // Tesla

// Derived parameters
long double C1 = 3 * mu0 / (4 * pi);
long double C2 = 6 * pi * eta * R0;
long double C3 = M0 * g;
long double D = R * T / (6 * Na * pi * R0 * eta);
//long double r0 = sqrt(3 * 2 * D * dt);
long double r0mod = sqrt(3 * 2 * D); // needs extra * dt^(1/2.0)
//long double C4 = r0;
long double C5 = mu0 / (4 * pi);

long double Mself = m0 / Vself;

long double Hself = - (1 / 3.0) * Mself;
long double Bself = mu0 * (Hself + Mself);

// 0.740 is atomic packing factor of the fcc lattice
// total mean total vol. of the clusters
long double Vtot = (4 / 3.0) * pi * pow(R0, 3) * (1 / 0.740) * pN;

long double C6 = ro0 * Vself * g;

//long double gamma_e = 1.7609 * (1E11); //s^-1 T^-1
//long double alpha   = 0.05; //magnetization dynamic damping