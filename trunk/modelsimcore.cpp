#include "modelsimcore.h"
#include "modelparameters.h"

#include <QObject>
#include <QThread>
#include <QDebug>
#include <QSettings>

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
    qDebug() << "Fundamental constants init";
    qDebug() << "--------------------------------";
    qDebug() << "pi =" << pi;
    qDebug() << "mu0 =" << mu0;
    qDebug() << "Na =" << Na;
    qDebug() << "kb =" << kb;
    qDebug() << "R =" << R;

    qDebug() << "";
    qDebug() << "Experimental setup";
    qDebug() << "--------------------------------";
    QSettings parameters(QString("parameters.ini"), QSettings::IniFormat);

    double Lx = parameters.value("Vessel/Lx", "0").toDouble();
    double Ly = parameters.value("Vessel/Ly", "0").toDouble();
    double Lz = parameters.value("Vessel/Lz", "0").toDouble();
    qDebug() << "Lx =" << Lx;
    qDebug() << "Ly =" << Ly;
    qDebug() << "Lz =" << Lz;

    double T = parameters.value("Fluid/T", "0").toDouble();
    qDebug() << "T =" << T;

    double B1 = parameters.value("Field/B1", "0").toDouble();
    double B2 = parameters.value("Field/B2", "0").toDouble();
    qDebug() << "B1 =" << B1;
    qDebug() << "B2 =" << B2;
}
