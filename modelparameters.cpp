#include <QDebug>
#include <QSettings>

#include "modelparameters.h"

// Fundamental constants [Metric System]
// -------------------------------------

// 1. Math and geometry
const double pi = acos(-1.0); // []

// 2. Physics

// 2.1. Electrodynamics

// Vacuum permeability
const double mu0 = 4 * pi * 1E-7; // [V s(A m)-1]

// Bohr magneton
const double muB = 9.27400968 * 1E-24; // [J T−1]

// 2.2. Thermodynamics

// Avogadro constant
const double Na = 6.02214129 * 1E23; // [mol−1]

// Boltzmann constant
const double kb = 1.3806488 * 1E-23; // [m2 kg s-2 K-1]

// Gas constant
const double R = Na * kb; // [J K−1 mol−1]

ModelParameters::ModelParameters(QString fileName)
{
    qDebug() << "Fundamental constants init";
    qDebug() << "--------------------------------";
    PWRITE(pi)
    PWRITE(mu0)
    PWRITE(Na)
    PWRITE(kb)
    PWRITE(R)

    qDebug() << "";
    qDebug() << "Experimental setup parameters";
    qDebug() << "--------------------------------";
    QSettings parameters(fileName, QSettings::IniFormat);

    PREAD(Vessel, Lx, 1E-6) //[µm]
    PREAD(Vessel, Ly, 1E-6) //[µm]
    PREAD(Vessel, Lz, 1E-6) //[µm]

    PREAD(Carrier fluid, rho_f, 1) //[kg m-3]
    PREAD(Carrier fluid, t1, 1) // [°C]
    PREAD(Carrier fluid, t2, 1) // [°C]

    PREAD(Ferrofluid, d_mean, 1E-9) // [nm]
    PREAD(Ferrofluid, s_mean, 1) // [A m2]

    PREAD(Nanoparticle, m_mol_rel, 1) // []
    PREAD(Nanoparticle, rho_p, 1) // [kg m-3]
    PREAD(Nanoparticle, s_mol, muB) // [Bohr magnetons]
    PREAD(Nanoparticle, Ms_p, 1E3) // [G]

    qDebug() << "";
    qDebug() << "Analysis";
    qDebug() << "--------------------------------";

    double V_particle_mean = pi * pow(d_mean, 3) / 6.0;
    PWRITE(V_particle_mean)
    double m_mean = rho_p * V_particle_mean;
    PWRITE(m_mean)
    double m_mol = m_mol_rel * 1E-3 / Na;
    double N_mol = m_mean / m_mol;
    PWRITE(N_mol)
    double s_particle_mean = s_mol * N_mol;
    PWRITE(s_particle_mean)
    double M_particle_mean = s_particle_mean / V_particle_mean;
    PWRITE(M_particle_mean)

    qDebug() << "--------------------------------";
}
