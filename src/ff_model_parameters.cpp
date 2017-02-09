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
#include "ff_iniParam.h"
#include <fstream>
#include <iostream>

using namespace std;
int isAnalysis = ff_paramAnalysis();
double gl_scale = iniGet("ParticleDiam","gl_scale");//1.0; // particle diameter scale parameter
	//

// Math constants
double pi = acos(-1.0);

// Physics constants
const double mu0 = 4 * pi * 1E-7;
const double muB =iniGet("PhysConst","muB");// 9.27400968 * 1E-24; // Bohr magneton
const double R =iniGet("PhysConst","R");// 8.31;
const double Na =iniGet("PhysConst","Na");// 6.02214179 * 1E23;
double g =iniGet("PhysConst","g");// 9.81;
const double kb =iniGet("PhysConst","kb");// 1.3806488 * 1E-23; // [m2 kg s-2 K-1]
const double ta0 =iniGet("PhysConst","ta0");// -273.15; // [°C]
const double gamma_e =iniGet("PhysConst","gamma_e");// 1.760859708 * 1E11; // [s^-1 T^-1] // electron gyromagnetic ratio

// Material
double a0 =iniGet("Material","a0");// 0.8397E-9; // [m] // magnetite unit cell size - a cubic spinel structure with space group Fd3m (above the Verwey temperature) // \cite{Fleet1981}

// Space
double volume_reduce =iniGet("Space","volume_reduce");// 0.15; //0.25; // 0.09; //0.085; // 0.1862; // initial density set
double scale = 2 * 7 * 0.15 * volume_reduce * gl_scale / pow(500.0 / pN, 1 / 3.0); // / pow(50.0, 1 / 3.0);
//double Lx = 1E-6 * scale, Ly = 1E-6 * scale, Lz = 25E-6 * scale; //meters
double nano_size =iniGet("Space","nano_size");// 950;
double Lx = nano_size * 1E-9 / 3.0, Ly = nano_size * 1E-9, Lz = nano_size * 1E-9 * 3; //meters
double delta_r = a0 * 0.5; // minimal distance between particles // order of magnitude of the oleic acid molecule width
double delta_r_init =iniGet("Space","delta_r_init");// 20 * 1E-9; // delta_r;

// periodic boundary conditions
int is_periodic =(int)iniGet("PeriodBoundary","is_periodic");// 0;

//double kExtra = 0.27 * 20 * 1; // change only this coef. instead Lz
//double Lx = 10 * kExtra * 1E-6, Ly = kExtra * 1E-6, Lz = kExtra * 1E-6;

// Basic physical model parameters
double dt_neel =iniGet("BasicPhysModelPar","dt_neel");// 1E-9; // [s] // Neel relaxation time threashold for d ~ 10 nm, (Fertman-p-62)
const double dt0 =iniGet("BasicPhysModelPar","dt0");// 1.088* 1E-8; //1E-1 * dt_neel; // 1E2 * dt_neel; // 100 * dt_neel; // s // old approach 1.5E-5
double d_neel =iniGet("BasicPhysModelPar","d_neel");// 10 * 1E-9; // [m] // Neel-to-Brown relaxation diameter threashold, (Fertman-p-62)
double d_min =iniGet("BasicPhysModelPar","d_min");// 2 * 1E-9; // [m] // minimal diameter of particle, wher Ms and T_curie is identical to ones in the bulk material, (Fertman-p-43)
long k_bm_inst_max =(long)iniGet("BasicPhysModelPar","k_bm_inst_max");// 100; // coefficient of a brownian motion instantiation: dt_bm_inst = dt * k_bm_inst_max
long k_bm_inst =(long)iniGet("BasicPhysModelPar","k_bm_inst");// 1;
double k_force_adapt_0 =iniGet("BasicPhysModelPar","k_force_adapt_0");// 1.00; // 1.05 is a regular value for the adaptive force model // 1.00 means an overdamped model at the dt >> m / gamma
long slow_steps =(long)iniGet("BasicPhysModelPar","slow_steps");// 0;
//double smooth_v = 10; // disabled in code
double smooth_r =iniGet("BasicPhysModelPar","smooth_r");// 0.4;
//double m_h_eff_tol = 1; // max. angle [rad] between m and B

double T =iniGet("BasicPhysModelPar","T");// 273.15; // K
double sigma_sf_nano =iniGet("BasicPhysModelPar","sigma_sf_nano");// 5E-4; //1E-4;

double kr = gl_scale; // particle size parameter []
//double R00 = 0.5 * 15E-9; // Radius of the nanoparticle [m]
double delta =iniGet("BasicPhysModelPar","delta");// 2.0E-9;
//double R0 = R00 + delta; // Radius including the acid sphere [m]
//double Vself = (4 * pi / 3.0) * pow(R00, 3); // [m^3]

int is_uniform_field_test =(int)iniGet("BasicPhysModelPar","is_uniform_field_test");// 0;

// Parameters of oleic acid drop
int is_large_mode =(int)iniGet("OleicDrop","is_large_mode");// 1; // largest particles mode
double large_fraction =iniGet("OleicDrop","large_fraction");// 7.5E-2; // 0.1 (Ivanov, Phase separation in bidisperse ferrocolloids) // 6.92E-02 (my sim); // fraction of the largest particles which form the primary aggregate (circle) in case of mode is_large_mode == 0
double k_large =iniGet("OleicDrop","k_large");// 1.0; // 0.91; // 0.925; // correction of the large particles size
int is_oleic =(int)iniGet("OleicDrop","is_oleic");// 0;
double R_oleic_0 = (Lx / 2.0);

// microdrop mode parameters
int isMicroDrop =(int)iniGet("MicroDropModePar","isMicroDrop");// 0; // enable the microdrop mode
double phi_v =iniGet("MicroDropModePar","phi_v");// 0.12; // volume concentration of the dispersed phase [Padalka, exp]
double alpha =iniGet("MicroDropModePar","alpha");// 0.5; // default saturation of the microdrop

int isPGasMode =(int)iniGet("MicroDropModePar","isPGasMode");// 0;
int isPGasPrevail =(int)iniGet("MicroDropModePar","isPGasPrevail");// 0;
double P_pgas =iniGet("MicroDropModePar","P_pgas");// 0; // particles gas pressure inside the oleic drop
double P_sf_oleic =iniGet("MicroDropModePar","P_sf_oleic");// 0; // oleic droplet surface tension pressure
//double eta_oleic = 25.6 * 1E-3; // [Pa * s]
//double eta_oleic = 14.285 * 1E-3; // [Pa * s]
double eta_oleic =iniGet("MicroDropModePar","eta_oleic");// 0;
double a3_eta_oleic =iniGet("MicroDropModePar","a3_eta_oleic");// - 1E-07; // my approximation of DOI: 10.1007/s11746-000-0197-z
double a2_eta_oleic =iniGet("MicroDropModePar","a2_eta_oleic");// 2E-05;
double a1_eta_oleic =iniGet("MicroDropModePar","a1_eta_oleic");// - 0.0018;
double a0_eta_oleic =iniGet("MicroDropModePar","a0_eta_oleic");// 0.0559;
//double sigma_sf = 32.5 * 1E-3; // [N / m]
double sigma_sf =iniGet("MicroDropModePar","sigma_sf");// 0;
double a_sigma_sf =iniGet("MicroDropModePar","a_sigma_sf");// 34.060119 * 1E-3; // [N / m] // linear coefficient: sigma_sf = a + b * t
double b_sigma_sf =iniGet("MicroDropModePar","b_sigma_sf");// - 0.061298 * 1E-3; // [N / m] // linear coefficient
double ro_oleic =iniGet("MicroDropModePar","ro_oleic");// 853;// density of oleic acid kg/m^3
double mol_mass_oleic =iniGet("MicroDropModePar","mol_mass_oleic");// 282*1E-3; //mol mass of oleic acid C18H34O2
double v_oleic =iniGet("MicroDropModePar","v_oleic");//  2* 1E-7; //m^3 volume of oleic acid
double mass_oleic = ro_oleic * v_oleic;

//double sigma_sf_nano = sigma_sf * 1 * 5.017559E-04; // [N / m]
//double sigma_sf_nano = 5E-2 * sigma_sf;

double rop =iniGet("MicroDropModePar","rop");// 5240; // magnetite mass density [kg / m^3]
//double M0 = Vself * rop;  // mass [kg]
double A_H =iniGet("MicroDropModePar","A_H");// 4E-20; // Hamaker constant [J]

double N_oa =iniGet("MicroDropModePar","N_oa");// 1E18; // Surface density of oleic acid at the 50% coating [m-2] //[Fertman]
double k_o =iniGet("MicroDropModePar","k_o");// 0.5; // 5E-4 - same result //0.1; // Level of coverage of surface by the oleic acid
double N_o; // Surface density of oleic acid

//double G_barrier = pow(kr, 2) * 25 * kb * T; // [TEMP] barrier which prevents particles aggregation
double K1 =iniGet("MicroDropModePar","K1");// 1.35 * 1E4 ; // [J/m] // First constant of the magnetite crystallographic anisotropy // (Goya2003)
double Ms_mass = 80 /* emu / g */ * (1E3) /* emu / kg */ * (1 / (9.274009 * (1E-21))); /* Bohr magnetons / kg */
//double m0 = Ms_mass * M0 /* Bohr magnetons */* 927.400915 * (1E-26); // Magnetic moment [J / T]
double Ms =iniGet("MicroDropModePar","Ms");// 478 * 1E3; // [A / m]

double Ch =iniGet("MicroDropModePar","Ch");// 0.01; // [DEPRECATED] adhesion / magnetic relation
double Ch_ss =iniGet("MicroDropModePar","Ch_ss");// 0; // 3E4; // 50 * 1E5; // soft-sphere repulsion parameter (parameter of the numerical model)

int load_at_start =(int)iniGet("MicroDropModePar","load_at_start");// 0;
int auto_save =(int)iniGet("MicroDropModePar","auto_save");// 1;
int manual_field_control =(int)iniGet("MicroDropModePar","manual_field_control");// 1; // 0-1-2-3 keys control to skip, Bx+, By+, Bz+ control 
int ext_field_is_homo =(int)iniGet("MicroDropModePar","ext_field_is_homo");// 1;
int auto_reversal =(int)iniGet("MicroDropModePar","auto_reversal");// 0;

int setting_plot =(int)iniGet("MicroDropModePar","setting_plot");// 1; // cluster creation plot m_tot / m0 and I.

double start_t =iniGet("MicroDropModePar","start_t");// 30 /* [micro_s] */ * (1E-6); // [s]
double T_ext =iniGet("MicroDropModePar","T_ext");// 0.2 * (1E6); // [micro_s] // external field period
double nu_ext = (1 / T_ext) * (1E6); // [Hz] // frequency of the external field (sin(w*t) dependence)

double start_ideal =iniGet("MicroDropModePar","start_ideal");// 1; // start chaos (ideal superparam. gas)
double start_sediment =iniGet("MicroDropModePar","start_sediment");// 0;
double ro0 =iniGet("MicroDropModePar","ro0");// 0.5 * (0.78 + 0.85) * 1E3; // kerosene density
double eta_car0 =iniGet("MicroDropModePar","eta_car0");// 0.00164; //Pa * s //kerosene (carrier liquid)
double eta_car =iniGet("MicroDropModePar","eta_car");// 1; //Pa * s //Mix eta
double mol_mass_car =iniGet("MicroDropModePar","mol_mass_car");// 170*1E-3; //molar mass of kerosene
double v_car =iniGet("MicroDropModePar","v_car");// 3*1E-6; //volume of kerosine
double mass_car = ro0*v_car;//mas of kerosine
//int __deprecated__brownian_shifts = 0;
//int __deprecated__brownian_force = 1;

//default order of magnitude of the external field but exact function is hardcoded 
double B0 = 2500 /*Oe*/ * 79.577 * mu0; // Tesla
double gradPerc =iniGet("DefaultOrder","gradPerc");// 5E-1;
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

double alpha_damp =iniGet("DefaultOrder","alpha_damp");// 0.05; //magnetization dynamic damping
int ScreenCaptureStep =(int)iniGet("DefaultOrder","ScreenCaptureStep");// 10000; //every ScreenCaptureStep`s steps will make screen  shot


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
	file1<<"volume_reduce"<<"=="<<volume_reduce<<endl;
	file1<<"scale"<<"=="<<scale<<endl;
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
	file1<<"k_bm_inst_max"<<"=="<<k_bm_inst_max<<endl;
	file1<<"k_bm_inst"<<"=="<<k_bm_inst<<endl;
	file1<<"k_force_adapt_0"<<"=="<<k_force_adapt_0<<endl;
	file1<<"slow_steps"<<"=="<<slow_steps<<endl;
	file1<<"smooth_r"<<"=="<<smooth_r<<endl;
	file1<<"T"<<"=="<<T<<endl;
	file1<<"sigma_sf_nano"<<"=="<<sigma_sf_nano<<endl;
	file1<<"kr"<<"=="<<kr<<endl;
	file1<<"delta"<<"=="<<delta<<endl;
	file1<<"is_uniform_field_test"<<"=="<<is_uniform_field_test<<endl;
	file1<<"is_large_mode"<<"=="<<is_large_mode<<endl;
	file1<<"large_fraction"<<"=="<<large_fraction<<endl;
	file1<<"k_large"<<"=="<<k_large<<endl;
	file1<<"is_oleic"<<"=="<<is_oleic<<endl;
	file1<<"R_oleic_0"<<"=="<<R_oleic_0<<endl;
	file1<<"isMicroDrop"<<"=="<<isMicroDrop<<endl;
	file1<<"phi_v"<<"=="<<phi_v<<endl;
	file1<<"alpha"<<"=="<<alpha<<endl;
	file1<<"isPGasMode"<<"=="<<isPGasMode<<endl;
	file1<<"isPGasPrevail"<<"=="<<isPGasPrevail<<endl;
	file1<<"P_pgas"<<"=="<<P_pgas<<endl;
	file1<<"P_sf_oleic"<<"=="<<P_sf_oleic<<endl;
	file1<<"eta_oleic"<<"=="<<eta_oleic<<endl;
	file1<<"a3_eta_oleic"<<"=="<<a3_eta_oleic<<endl;
	file1<<"a2_eta_oleic"<<"=="<<a2_eta_oleic<<endl;
	file1<<"a1_eta_oleic"<<"=="<<a1_eta_oleic<<endl;
	file1<<"a0_eta_oleic"<<"=="<<a0_eta_oleic<<endl;
	file1<<"sigma_sf"<<"=="<<sigma_sf<<endl;
	file1<<"a_sigma_sf"<<"=="<<a_sigma_sf<<endl;
	file1<<"b_sigma_sf"<<"=="<<b_sigma_sf<<endl;
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
	file1<<"Ch"<<"=="<<Ch<<endl;
	file1<<"Ch_ss"<<"=="<<Ch_ss<<endl;
	file1<<"load_at_start"<<"=="<<load_at_start<<endl;
	file1<<"auto_save"<<"=="<<auto_save<<endl;
	file1<<"manual_field_control"<<"=="<<manual_field_control<<endl;
	file1<<"ext_field_is_homo"<<"=="<<ext_field_is_homo<<endl;
	file1<<"auto_reversal"<<"=="<<auto_reversal<<endl;
	file1<<"setting_plot"<<"=="<<setting_plot<<endl;
	file1<<"start_t"<<"=="<<start_t<<endl;
	file1<<"T_ext"<<"=="<<T_ext<<endl;
	file1<<"nu_ext"<<"=="<<nu_ext<<endl;
	file1<<"start_ideal"<<"=="<<start_ideal<<endl;
	file1<<"start_sediment"<<"=="<<start_sediment<<endl;
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
