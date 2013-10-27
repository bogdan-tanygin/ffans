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

//void ff_io_save_hyst(void);
void ff_io_save_setting(ff_vect_t m_tot,double I);
void ff_io_autosave(void);
void ff_io_load(void);
/*void ff_io_save_susceptX(void);
void ff_io_save_susceptY(void);
void ff_io_save_susceptZ(void);*/

extern void cbMouseInput(int button, int state, int x, int y);
extern void cbMouseMove(int x, int y);
extern void cbSpecialKeyPressed(int key, int x, int y);
extern void cbKeyPressed(unsigned char key, int x, int y);
