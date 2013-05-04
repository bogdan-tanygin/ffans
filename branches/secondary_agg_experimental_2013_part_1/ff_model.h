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

// working variables 
////////////////////
extern ff_vect_t r[];  //particles positions
extern ff_vect_t m[];  //particles magneitc moment direction

extern ff_vect_t v[];

extern long number_of_mag_threads;
extern long mag_thread_limit_plus[]; // pN is a maximum number of threads. Index is id of thread.
extern long mag_thread_limit_minus[]; // pN is a maximum number of threads. Index is id of thread.
extern long mag_thread_p_id[];
extern long mag_thread_id_ordered_p[pN + 1][pN + 1]; // [thread_id][ordered position of particle inside thread] = p
extern long mag_thread_p_position[];
//extern ff_vect_t F0[];

//extern ff_vect_t B[]; //distribution of the magnetic field between the particles

extern double B_hyst[];
extern double Mz_hyst[];
extern double B_hyst_n[];
extern double Mz_hyst_n[];

extern long step;
extern double t;

extern double Ek;

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

extern int exist_p[];
extern double Rp[];