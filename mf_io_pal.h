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

void mf_io_save_hyst(void);
void mf_io_save_setting(mf_vect_t m_tot,long double I);
void mf_io_autosave(void);
void mf_io_load(void);

extern void cbMouseInput(int button, int state, int x, int y);
extern void cbMouseMove(int x, int y);
extern void cbSpecialKeyPressed(int key, int x, int y);
extern void cbKeyPressed(unsigned char key, int x, int y);
