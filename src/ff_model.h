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

#define MUL(A,B) (##A.x * ##B.x + ##A.y * ##B.y + ##A.z * ##B.z)

#include "ff_model_parameters.h"

typedef struct {
    double x,y,z;
} ff_vect_t;

void ff_model_init(void);
void ff_model_next_step(void);

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

extern long step;
extern double dt;
extern double t;
extern double T_mean;
extern long k_mean;

extern double Ek;
extern double Ek_tr;
extern double Ek_rot;

extern double kB;

extern double mz_tot;

extern double mz_glob;
extern ff_vect_t m_tot_glob;

extern int time_go;

extern double BmanX;
extern double BmanY;
extern double BmanZ;

extern ff_vect_t m_tot;

extern double V0_tot;
extern double V0_largest_EV; // mathematical expected value of largest particles total volume // see is_large_mode variable
extern double V0_tot_EV; // mathematical expected value of particles total volume
extern double I_glob;
extern double phi_vol_fract;

extern int exist_p[];
extern double Rp[];

// Update of the effective instantiated random force
void ff_model_effective_random_motion_update(long);
void ff_model_update_dT(void);

void ff_model_size_dispersion_init(void);
void ff_model_size_dispersion_param_calc(double,long);

int ff_model_check_walls(long);

double ff_model_G_steric(long p, long ps);
double ff_model_G_london(long p, long ps);
double ff_model_G_zeeman(long p);