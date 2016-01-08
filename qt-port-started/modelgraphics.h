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

#ifndef MODELGRAPH_H
#define MODELGRAPH_H

#include <QObject>
#include <QColor>
//#include <QGLWidget>

class Patch;
struct Geometry;

//! [0]
class ModelGraph : public QObject
{
public:
    explicit ModelGraph(QObject *parent, int d = 64, qreal s = 1.0);
    ~ModelGraph();
    void setColor(QColor c);
    void draw() const;
private:
    void buildGeometry(int d, qreal s);

    QList<Patch *> parts;
    Geometry *geom;

    QColor cartesianSystemColor;
};
//! [0]

#endif // MODELGRAPH_H
