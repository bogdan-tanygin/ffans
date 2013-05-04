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

typedef struct {
    long double r,g,b,a;
} mf_color4f_t;

void mf_mgr_init(void);

void mf_mgr_print_info(void);

void mf_mgr_show_next_step(void);

extern int show_m, show_b, show_bext;

extern int transp;

extern int show_sphere;