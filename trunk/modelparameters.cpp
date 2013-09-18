#include <QDebug>
#include <QSettings>

#include "modelparameters.h"

// Fundamental constants [Metric System]
// -------------------------------------

// 1. Math and geometry
const double pi = acos(-1.0); // []
const double n_dense_packing = pi / sqrt(18.0); // []

// 2. Physics

// 2.1. Electrodynamics

// Vacuum permeability
const double mu0 = 4 * pi * 1E-7; // [V s(A m)-1]

// Bohr magneton
const double muB = 9.27400968 * 1E-24; // [J T−1]

// 2.2. Thermodynamics

// Absolute zero
const double ta0 = -273.15; // [°C]

// Avogadro constant
const double Na = 6.02214129 * 1E23; // [mol−1]

// Boltzmann constant
const double kb = 1.3806488 * 1E-23; // [m2 kg s-2 K-1]

// Gas constant
const double R = Na * kb; // [J K−1 mol−1]

ModelParameters::ModelParameters(QString fileName)
{
    qDebug() << "I. Fundamental constants init";
    qDebug() << "=================================";
    PWRITE(pi)
    PWRITE(n_dense_packing)
    PWRITE(mu0)
    PWRITE(Na)
    PWRITE(kb)
    PWRITE(R)

    qDebug() << "";
    qDebug() << "II. Experimental setup parameters";
    qDebug() << "=================================";
    QSettings parameters(fileName, QSettings::IniFormat);

    PREAD(Vessel, Lx, 1E-6) //[µm]
    PREAD(Vessel, Ly, 1E-6) //[µm]
    PREAD(Vessel, Lz, 1E-6) //[µm]

    PREAD(Carrier fluid, rho_f, 1) //[kg m-3]
    PREAD(Carrier fluid, t1, 1) // [°C]
    PREAD(Carrier fluid, t2, 1) // [°C]

    double T1 = t1 - ta0;
    PWRITE(T1)
    double T2 = t2 - ta0;
    PWRITE(T2)

    PREAD(Ferrofluid, d_mean, 1E-9) // [nm]
    PREAD(Ferrofluid, s_mean, 1) // [A m2]
    PREAD(Ferrofluid, n_p, 1) // [m-3]

    PREAD(Nanoparticle, m_mol_rel, 1) // []
    PREAD(Nanoparticle, rho_p, 1) // [kg m-3]
    PREAD(Nanoparticle, s_mol, muB) // [Bohr magnetons]
    PREAD(Nanoparticle, Ms_p, 1E3) // [G]

    PREAD(Field, H1, 1.0 / (4.0 * pi * 1E-3)) // [Oe]
    PREAD(Field, H_med, 1.0 / (4.0 * pi * 1E-3)) // [Oe]
    PREAD(Field, H2, 1.0 / (4.0 * pi * 1E-3)) // [Oe]

    qDebug() << "";
    qDebug() << "III. Analysis";
    qDebug() << "=================================";
    qDebug() << "Nanoparticle parameters validation";

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
    qDebug() << "Derived paramaters";
    qDebug() << "";

    qDebug() << "---";
    qDebug() << "Estimation of dipole-dipole interaction value";

    qDebug() << "Mean distance between nanoparticles. Simple packing. In case of cubic packing, we will have larger distance.";
    double l = 1 / pow(n_p, 1/3.0);
    PWRITE(l)
    double Wdd_order = (mu0 / (4.0 * pi)) * pow(s_mean, 2) / pow(l, 3);
    PWRITE(Wdd_order)
    double beta2 = 1 / (kb * T2);
    PWRITE(beta2)
    qDebug() << "Small parameters x";
    double x = beta2 * Wdd_order;
    PWRITE(x)

    qDebug() << "---";
    qDebug() << "Estimation of external field value";

    double Wz_order = mu0 * s_mean * H2;
    PWRITE(Wz_order)
    double y = beta2 * Wz_order;
    PWRITE(y)

    qDebug() << "---";
    qDebug() << "Model of Lennard-Jones fluid - critical  parameters";

    qDebug() << "Primary aggregate diameter estimation";

    double d_mean_primary = 285 * 1E-9;
    PWRITE(d_mean_primary)
    double N_p_primary = n_dense_packing * pow(d_mean_primary, 3) / pow(d_mean, 3);
    PWRITE(N_p_primary)
    double ks_primary = N_p_primary * 0.75 / 100.0; // [my original]
    if (ks_primary < 1.0) ks_primary = 1.0; // [Dzyan]
    PWRITE(ks_primary)

    double tc = ta0 + sqrt(1.3 / 3.0) * (mu0 / (4.0 * pi)) *
            pow(ks_primary * s_mean, 2) / (kb * pow(d_mean_primary, 3));
    PWRITE(tc)

    qDebug() << "Different particles diameter estimation";

    qDebug() << "--------------------------------";
}
