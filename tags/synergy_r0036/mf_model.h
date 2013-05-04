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

#define MUL(A,B) (##A.x * ##B.x + ##A.y * ##B.y + ##A.z * ##B.z)

typedef struct {
    long double x,y,z;
} mf_vect_t;

void mf_model_init(void);
void mf_model_next_step(void);
void mf_model_check_collisions(long);

void mf_model_upgrade_ext_field(void);

mf_vect_t Bext(void);

// working variables 
////////////////////
extern mf_vect_t r[];  //particles positions
extern mf_vect_t m[];  //particles magneitc moment direction

extern mf_vect_t v[];
//extern mf_vect_t F0[];

//extern mf_vect_t B[]; //distribution of the magnetic field between the particles

extern long double B_hyst[];
extern long double Mz_hyst[];
extern long double B_hyst_n[];
extern long double Mz_hyst_n[];

extern long step;
extern long double t;

extern long double Ek;

extern long double kB;

extern int hyst_mode;

extern long double mz_tot;

extern long glob_start_step;
extern long double glob_start_t;
extern long double mz_glob;

extern int time_go;



extern int g_hyst_start_line;
extern int g_hyst_up_line;
extern int g_hyst_bottom_line;


extern long double g_Bz_prev;