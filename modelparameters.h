/******************************************************************************************
* Copyright (C) 2011,2013-2015 Dr. Bogdan Tanygin <b.m.tanygin@gmail.com> and contributors.
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
******************************************************************************************/

#ifndef MODELPARAMETERS_H
#define MODELPARAMETERS_H

#include <math.h>
#include <QString>

extern const double pi;
extern const double n_dense_packing;
extern const double mu0;
extern const double muB;
extern const double ta0;
extern const double Na;
extern const double kb;
extern const double R;

#define PWRITE(P) par = #P; qDebug() << par + " = " + QString::number(##P, 'E');

// reading the parameters from .ini file. "C" is a to-metric-system conversion factor.
#define PREAD(O,P,C) obj = #O; par = #P; double (##P) = parameters.value(obj+"/"+par, "0").toDouble() * ##C;\
    qDebug() << par + " = " + QString::number(##P, 'E');

class ModelParameters
{
public:
    ModelParameters(QString fileName);
private:
    QString obj; // physical object in the parameters .ini file
    QString par; // paramater of physical object in the parameters .ini file
};

#endif // MODELPARAMETERS_H
