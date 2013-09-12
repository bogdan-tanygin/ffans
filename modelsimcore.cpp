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

    double d_mean = parameters.value("Fluid/d_mean", "0").toDouble();
    qDebug() << "d_mean =" << d_mean;
    double m_mol_rel = parameters.value("Nanoparticle/m_mol_rel", "0").toDouble();
    qDebug() << "m_mol_rel =" << m_mol_rel;
    double s_mol = parameters.value("Nanoparticle/s_mol", "0").toDouble();
    qDebug() << "s_mol =" << s_mol;
    double rho = parameters.value("Nanoparticle/rho", "0").toDouble();
    qDebug() << "rho =" << rho;

    qDebug() << "";
    qDebug() << "Analysis";
    qDebug() << "--------------------------------";

    double V_particle = pi * pow(d_mean, 3) / 6.0;
    double m_mean_temp = rho * V_particle;
    double m_mol = m_mol_rel * 1E-3 / Na;

    double s_particle = muB * s_mol * m_mean_temp / m_mol;
    qDebug() << "s_particle =" << s_particle;
    double M_particle = s_particle / V_particle;
    qDebug() << "M_particle =" << M_particle;
    qDebug() << "relation of M #1 =" << M_particle / 46000;
}
