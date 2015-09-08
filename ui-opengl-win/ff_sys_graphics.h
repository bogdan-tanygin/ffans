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

extern int window_id;

extern float x_rot; 
extern float y_rot;   
extern float x_speed; 
extern float y_speed; 
extern int rotating;  
extern float z_off;

extern float frame_rate;
extern int frame_count;

extern int window_width;
extern int window_height;

extern int projection_type;

void ff_gr_init(int argc, char **argv);
void ff_gr_print(void *font, char *str);
void ff_gr_loop(void);
void cbResizeScene(int width, int height);








