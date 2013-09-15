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
