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

#include "modelsimcore.h"
#include "modelparameters.h"

#include <QObject>
#include <QThread>
#include <QDebug>

void Worker::doWork(const QString &parameter){
    //double i = 0;
    //while (i++) qDebug() << "thread " << i;
}

SimController::SimController() {
        Worker *worker = new Worker;
        worker->moveToThread(&workerThread);
        connect(&workerThread, &QThread::finished, worker, &QObject::deleteLater);
        connect(this, &SimController::operate, worker, &Worker::doWork);
        //connect(worker, &Worker::resultReady, this, &SimController::handleResults);
        workerThread.start();
    }

ModelSimCore::ModelSimCore()
{
}

void ModelSimCore::modelSimInit(void)
{
    ModelParameters* modelParameters = new ModelParameters (QString("parameters.ini"));
}
