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

#define MUL(A,B) (##A.x * ##B.x + ##A.y * ##B.y + ##A.z * ##B.z)

#include "ff_model_parameters.h"

typedef struct {
    double x,y,z;
} ff_vect_t;

void ff_model_init(void);
void ff_model_next_step(void);
void ff_model_check_collisions(long);

void ff_model_upgrade_ext_field(void);

ff_vect_t Bext(double x, double y, double z);
ff_vect_t B_het_bem(ff_vect_t r);
ff_vect_t dBz_het_bem(ff_vect_t r);

// working variables 
////////////////////
extern ff_vect_t r[];  //particles positions
extern ff_vect_t m[];  //particles magneitc moment direction
extern ff_vect_t v[];
extern ff_vect_t w[];

extern double Rp0[];
extern double Rp[];

extern double G_dd[];

//extern ff_vect_t P[];  // effective instantiated random force

//extern ff_vect_t F0[];

//extern ff_vect_t B[]; //distribution of the magnetic field between the particles

extern double B_hyst[];
extern double Mz_hyst[];
extern double B_hyst_n[];
extern double Mz_hyst_n[];

extern long step;
extern double dt;
extern double t;
extern double dT;
extern double T_mean;
extern long k_mean;
extern double k_force_adapt_mean_print;
extern double T_mean_loc_prev;

extern double Ek;
extern double Ek_tr;
extern double Ek_rot;

extern double kB;

extern int hyst_mode;

extern double mz_tot;

extern long glob_start_step;
extern long glob_start_step_susc;
extern double glob_start_t;
extern double mz_glob;
extern ff_vect_t m_tot_glob;

extern int time_go;

extern int g_hyst_start_line;
extern int g_hyst_up_line;
extern int g_hyst_bottom_line;

extern double g_Bz_prev;

extern double BmanX;
extern double BmanY;
extern double BmanZ;

extern ff_vect_t m_tot;

extern double V0_tot;
extern double V0_largest_EV; // mathematical expected value of largest particles total volume // see is_large_mode variable
extern double V0_tot_EV; // mathematical expected value of particles total volume
extern double I_glob;

extern int exist_p[];
extern double Rp[];

extern long pN_oleic_drop;
extern long pN_oleic_drop_I;
extern long pN_oleic_drop_II;
extern long pN_oleic_drop_III;
extern double phi_vol_fract_oleic;

extern double k_delta_force_rel_tot;
extern double k_delta_force_rel_p;
extern double k_delta_torque_rel_tot;
extern double k_delta_torque_rel_p;
extern double k_force_adapt;

extern double R_oleic;

// Update of the effective instantiated random force
void ff_model_effective_random_motion_update(long);
void ff_model_update_dT(void);
void ff_model_update_dT_p(long);

void ff_model_size_dispersion_init(void);
void ff_model_size_dispersion_param_calc(double,long);

int ff_model_check_walls(long);

//void ff_model_check_overlapp(long);

void ff_model_brownian_validation(long);
void ff_model_update_conc_in_oleic(long);
void ff_model_update_mdrop_parameters(void);

double ff_model_G_steric(long p, long ps);
double ff_model_G_london(long p, long ps);
double ff_model_G_zeeman(long p);