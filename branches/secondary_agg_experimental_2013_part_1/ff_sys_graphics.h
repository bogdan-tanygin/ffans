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

class ff_gr_params_t
{   
public:
    int window_id;
    
    float x_rot; 
    float y_rot;   
    float x_speed; 
    float y_speed; 
    int rotating;  
    float z_off;
    
    float frame_rate;
    int frame_count;
};

extern int window_width;
extern int window_height;

void ff_gr_init(int argc, char **argv, ff_gr_params_t *ff_gr_params);
void ff_gr_print(void *font, char *str);
void ff_gr_loop(void);








