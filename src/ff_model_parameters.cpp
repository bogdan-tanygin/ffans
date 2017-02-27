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

// SYSTEM SI

#include <math.h>
#include "ff_model_parameters.h"
#include "ff_iniParse.h"
#include <fstream>
#include <iostream>

using namespace std;

double gl_scale = 1.0;//1.0; // particle diameter scale parameter
int isAutoSetPosition = 0;
// Math constants
double pi = acos(-1.0);
int isShowInfo = 0;

// Physics constants
const double mu0 = 4 * pi * 1E-7;
double muB =1.0;// 9.27400968 * 1E-24; // Bohr magneton
double R =1.0;// 8.31;
double Na =1.0;// 6.02214179 * 1E23;
double g =1.0;// 9.81;
double kb =1.0;// 1.3806488 * 1E-23; // [m2 kg s-2 K-1]
double ta0 =1.0;// -273.15; // [°C]
double gamma_e =1.0;// 1.760859708 * 1E11; // [s^-1 T^-1] // electron gyromagnetic ratio

// Material
double a0 =1.0;// 0.8397E-9; // 

// Space
//double volume_reduce =iniGet("Space","volume_reduce");// 0.15; //0.25; // 0.09; //0.085; // 0.1862; // initial density set
//double scale = 2 * volume_reduce * gl_scale / pow(500.0 / pN, 1 / 3.0); // / pow(50.0, 1 / 3.0);
//double Lx = 1E-6 * scale, Ly = 1E-6 * scale, Lz = 25E-6 * scale; //meters
double nano_size =1.0;// 950;
double Lx = 1.0, Ly = 1.0, Lz = 1.0; //meters
double delta_r = 1; // minimal distance between particles // order of magnitude of the oleic acid molecule width
double delta_r_init =1.0;// 20 * 1E-9; // delta_r;

// periodic boundary conditions
int is_periodic =0;// 0;

// Basic physical model parameters
double dt_neel =1.0;// 1E-9; // [s] // Neel relaxation time threashold for d ~ 10 nm, (Fertman-p-62)
double dt0 =1.0;// 1.088* 1E-8; //1E-1 * dt_neel; // 1E2 * dt_neel; // 100 * dt_neel; // s // old approach 1.5E-5
double d_neel =1.0;// 10 * 1E-9; // [m] // Neel-to-Brown relaxation diameter threashold, (Fertman-p-62)
double d_min =1.0;// 2 * 1E-9; // [m] // minimal diameter of particle, wher Ms and T_curie is identical to ones in the bulk material, (Fertman-p-43)
long slow_steps =1.0;// 0;
double smooth_r =1.0;// 0.4;

double T =1.0;// 273.15; // K

double kr = 1.0; // particle size parameter []
double delta =1.0;// 2.0E-9;

int is_large_mode =1;// 1; // largest particles mode
double large_fraction =1.0;// 7.5E-2; // 0.1 (Ivanov, Phase separation in bidisperse ferrocolloids) // 6.92E-02 (my sim); // fraction of the largest particles which form the primary aggregate (circle) in case of mode is_large_mode == 0
double k_large =1.0;// 1.0; // 0.91; // 0.925; // correction of the large particles size

// microdrop mode parameters

double eta_oleic =1.0;// 0;
double a3_eta_oleic =1.0;// - 1E-07; // my approximation of DOI: 10.1007/s11746-000-0197-z
double a2_eta_oleic =1.0;// 2E-05;
double a1_eta_oleic =1.0;// - 0.0018;
double a0_eta_oleic =1.0;// 0.0559;
double ro_oleic =1.0;// 853;// density of oleic acid kg/m^3
double mol_mass_oleic =1.0;// 282*1E-3; //mol mass of oleic acid C18H34O2
double v_oleic =1.0;//  2* 1E-7; //m^3 volume of oleic acid
double mass_oleic = 1.0;

double rop =1.0;// 5240; // magnetite mass density [kg / m^3]
double A_H =1.0;// 4E-20; // Hamaker constant [J]

double N_oa =1.0;// 1E18; // Surface density of oleic acid at the 50% coating [m-2] //[Fertman]
double k_o =1.0;// 0.5; // 5E-4 - same result //0.1; // Level of coverage of surface by the oleic acid
double N_o; // Surface density of oleic acid

double K1 =1.0;// 1.35 * 1E4 ; // [J/m] // First constant of the magnetite crystallographic anisotropy // (Goya2003)
double Ms_mass = 80 /* emu / g */ * (1E3) /* emu / kg */ * (1 / (9.274009 * (1E-21))); /* Bohr magnetons / kg */
double Ms =1.0;// 478 * 1E3; // [A / m]

int load_at_start =0;// 0;
int auto_save =1;// 1;
int manual_field_control =1;// 1; // 0-1-2-3 keys control to skip, Bx+, By+, Bz+ control 
int ext_field_is_homo =1;// 1;
int auto_reversal =0;// 0;

int setting_plot =1;// 1; // cluster creation plot m_tot / m0 and I.

double start_ideal =1.0;// 1; // start chaos (ideal superparam. gas)
double ro0 =1.0;// 0.5 * (0.78 + 0.85) * 1E3; // kerosene density
double eta_car0 =1.0;// 0.00164; //Pa * s //kerosene (carrier liquid)
double eta_car =1.0;// 1; //Pa * s //Mix eta
double mol_mass_car =1.0;// 170*1E-3; //molar mass of kerosene
double v_car =1.0;// 3*1E-6; //volume of kerosine
double mass_car = ro0*v_car;//mas of kerosine

//default order of magnitude of the external field but exact function is hardcoded 
double B0 = 1.0; // Tesla
double gradPerc =1.0;// 5E-1;
double gradL = 1.0;

// Derived parameters
double C1 = 1.0;
double C5 = 1.0;

double alpha_damp =1.0;// 0.05; //magnetization dynamic damping
int ScreenCaptureStep =1;// 10000; //every ScreenCaptureStep`s steps will make screen  shot


void ParamInfo()
{
	fstream file1 = fstream("info2.txt");
	file1<<"infoTest"<<endl;
	file1<<"gl_scale"<<"=="<<gl_scale<<endl;
	file1<<"pi"<<"=="<<pi<<endl;
	file1<<"mu0"<<"=="<<mu0<<endl;
	file1<<"muB"<<"=="<<muB<<endl;
	file1<<"R"<<"=="<<R<<endl;
	file1<<"Na"<<"=="<<Na<<endl;
	file1<<"g"<<"=="<<g<<endl;
	file1<<"kb"<<"=="<<kb<<endl;
	file1<<"ta0"<<"=="<<ta0<<endl;
	file1<<"gamma_e"<<"=="<<gamma_e<<endl;
	file1<<"a0"<<"=="<<a0<<endl;
	file1<<"nano_size"<<"=="<<nano_size<<endl;
	file1<<"Lx"<<"=="<<Lx<<endl;
	file1<<"Ly"<<"=="<<Ly<<endl;
	file1<<"Lz"<<"=="<<Lz<<endl;
	file1<<"delta_r"<<"=="<<delta_r<<endl;
	file1<<"delta_r_init"<<"=="<<delta_r_init<<endl;
	file1<<"is_periodic"<<"=="<<is_periodic<<endl;
	file1<<"dt_neel"<<"=="<<dt_neel<<endl;
	file1<<"dt0"<<"=="<<dt0<<endl;
	file1<<"d_neel"<<"=="<<d_neel<<endl;
	file1<<"d_min"<<"=="<<d_min<<endl;
	file1<<"slow_steps"<<"=="<<slow_steps<<endl;
	file1<<"smooth_r"<<"=="<<smooth_r<<endl;
	file1<<"T"<<"=="<<T<<endl;
	file1<<"kr"<<"=="<<kr<<endl;
	file1<<"delta"<<"=="<<delta<<endl;
	file1<<"is_large_mode"<<"=="<<is_large_mode<<endl;
	file1<<"large_fraction"<<"=="<<large_fraction<<endl;
	file1<<"k_large"<<"=="<<k_large<<endl;
	file1<<"eta_oleic"<<"=="<<eta_oleic<<endl;
	file1<<"a3_eta_oleic"<<"=="<<a3_eta_oleic<<endl;
	file1<<"a2_eta_oleic"<<"=="<<a2_eta_oleic<<endl;
	file1<<"a1_eta_oleic"<<"=="<<a1_eta_oleic<<endl;
	file1<<"a0_eta_oleic"<<"=="<<a0_eta_oleic<<endl;
	file1<<"ro_oleic"<<"=="<<ro_oleic<<endl;
	file1<<"mol_mass_oleic"<<"=="<<mol_mass_oleic<<endl;
	file1<<"v_oleic"<<"=="<<v_oleic<<endl;
	file1<<"mass_oleic"<<"=="<<mass_oleic<<endl;
	file1<<"rop"<<"=="<<rop<<endl;
	file1<<"A_H"<<"=="<<A_H<<endl;
	file1<<"N_oa"<<"=="<<N_oa<<endl;
	file1<<"k_o"<<"=="<<k_o<<endl;
	file1<<"N_o"<<"=="<<N_o<<endl;
	file1<<"K1"<<"=="<<K1<<endl;
	file1<<"Ms_mass"<<"=="<<Ms_mass<<endl;
	file1<<"Ms"<<"=="<<Ms<<endl;
	file1<<"load_at_start"<<"=="<<load_at_start<<endl;
	file1<<"auto_save"<<"=="<<auto_save<<endl;
	file1<<"manual_field_control"<<"=="<<manual_field_control<<endl;
	file1<<"ext_field_is_homo"<<"=="<<ext_field_is_homo<<endl;
	file1<<"auto_reversal"<<"=="<<auto_reversal<<endl;
	file1<<"setting_plot"<<"=="<<setting_plot<<endl;
	file1<<"start_ideal"<<"=="<<start_ideal<<endl;
	file1<<"ro0"<<"=="<<ro0<<endl;
	file1<<"eta_car0"<<"=="<<eta_car0<<endl;
	file1<<"eta_car"<<"=="<<eta_car<<endl;
	file1<<"mol_mass_car"<<"=="<<mol_mass_car<<endl;
	file1<<"v_car"<<"=="<<v_car<<endl;
	file1<<"mass_car"<<"=="<<mass_car<<endl;
	file1<<"B0"<<"=="<<B0<<endl;
	file1<<"gradPerc"<<"=="<<gradPerc<<endl;
	file1<<"gradL"<<"=="<<gradL<<endl;
	file1<<"C1"<<"=="<<C1<<endl;
	file1<<"C5"<<"=="<<C5<<endl;
	file1<<"alpha_damp"<<"=="<<alpha_damp<<endl;
	file1<<"ScreenCaptureStep"<<"=="<<ScreenCaptureStep<<endl;

}
