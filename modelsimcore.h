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

#ifndef MODELSIMCORE_H
#define MODELSIMCORE_H

#include <QObject>
#include <QThread>
#include <QDebug>

class Worker : public QObject
{
    Q_OBJECT
    QThread workerThread;

public slots:
    void doWork(const QString &parameter);

signals:
    void resultReady(const QString &result);
};

class SimController : public QObject
{
    Q_OBJECT
    QThread workerThread;
public:
    SimController();
    //~SimController();
public slots:
    //void handleResults(const QString &);
signals:
    void operate(const QString &);
};

class ModelSimCore
{
public:
    ModelSimCore();
    void modelSimInit(void);
};

#endif // MODELSIMCORE_H
