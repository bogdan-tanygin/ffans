#ifndef MODELPARAMETERS_H
#define MODELPARAMETERS_H

#include <math.h>
#include <QString>

extern const double pi;
extern const double mu0;
extern const double muB;
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
