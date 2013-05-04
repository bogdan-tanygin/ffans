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

// Math constants
extern long double pi;

// Physics constants
extern long double mu0;
extern long double R;
extern long double Na;
extern long double g;
extern long double kb;

extern long double C1;

extern long double kExtra;

// Geometry
extern long double Lx, Ly, Lz; //meters

// Basic physical model parameters
extern long double dt;
extern long slow_steps;
extern long double smooth_v;
extern long double smooth_r;
extern long double m_h_eff_tol;

extern long double Ch;
//extern long double EPS;
//extern long double ro1;
//extern long double ro2;
extern int auto_reversal;
extern int auto_save;
extern int setting_plot;
extern long double start_t; // [s]
extern long double nu_ext;

#define pN 230 // number of the particles

extern long double start_ideal;
extern long double start_sediment;

extern int brownian_enable;

//extern long double gap;

extern long double R00; // Radius of the nanoparticle [m]
extern long double R0; // Radius including the acid sphere [m]
extern long double M0; // mass [kg]
extern long double m0; // Magnetic moment [J / T]

extern long double T;
extern long double eta;

//order of magnitude of the external field but exact function is hardcoded 
extern long double B0; // Tesla

// Derived parameters
extern long double C2;
extern long double C3;
extern long double D;
//extern long double r0;
extern long double r0mod;
//extern long double C4;
extern long double C5;

extern long double Vself;
extern long double Mself;

extern long double Hself;
extern long double Bself;

extern long double Vtot;

extern long double C6;

extern long double gamma_e;
extern long double alpha;