/**************************************************************************
* Copyright (C) 2016,2017 Dmytro Matskevych<dimqqqq@mail.ru>
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
#include <math.h>
#include <string>
using namespace std;
extern double ff_visousity_mix(
								double y1,
								double mu1,
								double y2,
								double mu2);
extern double ff_molar_part (double nu1, double nu2);
extern double ff_mol(double m, double mol_m);
extern void GetScreenShot(string name1);
extern void addPosition(float x, float y, float z, int pts);
extern void ActiveWindow();
extern void delPosition();
extern void ChangePosition();
extern int MaxPointOfPosition;
extern int counterOfPosition;
extern void auto_set_position(int isAutoSet_);
extern void ff_pieces_coord_info();


