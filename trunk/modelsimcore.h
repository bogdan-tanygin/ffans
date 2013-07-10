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
};

#endif // MODELSIMCORE_H
