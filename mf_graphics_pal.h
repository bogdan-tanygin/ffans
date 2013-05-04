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

#define PROGRAM_TITLE "Synergy 1.4.5"

void mf_gr_print(void *font, char *str);
void mf_gr_init(int argc, char **argv);
void mf_gr_loop(void);

// Cube position and rotation speed variables.
extern float x_rot;   
extern float y_rot;   
extern float x_speed; 
extern float y_speed; 
extern int rotating;  
extern float z_off;

extern int window_id;
