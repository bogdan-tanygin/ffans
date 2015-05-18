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
    PREAD(Ferrofluid, t1, 1) // [°C]
    PREAD(Ferrofluid, t2, 1) // [°C]

    double T1 = t1 - ta0;
    PWRITE(T1)
    double T2 = t2 - ta0;
    PWRITE(T2)

    PREAD(Ferrofluid, sigma_sf_oleic, 1) // [N / m]
    PREAD(Ferrofluid, eta_oleic, 1) // [Pa * s]
    PREAD(Ferrofluid, d_mean, 1E-9) // [nm]
    double R_mean = d_mean / 2.0;
    PREAD(Ferrofluid, s_mean, 1) // [A m2]
    PREAD(Ferrofluid, n_p, 1) // [m-3]
    PREAD(Ferrofluid, t_vap_ldc, 1) // [°C]

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
    double v_particle_rms = sqrt(3 * kb * T1 / m_mean);
    PWRITE(v_particle_rms)

    qDebug() << "--------------------------------";
    qDebug() << "Derived paramaters";
    qDebug() << "";

    qDebug() << "---";
    qDebug() << "Estimation of dipole-dipole interaction value";

    qDebug() << "Mean distance between nanoparticles. Simple packing. In case of cubic packing, we will have larger distance.";
    double l = 1 / pow(n_p, 1 / 3.0);
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

    qDebug() << "[Hypothesis] Primary aggregate diameter estimation";

    double d_mean_primary = 285 * 1E-9;
    //double d_mean_primary = 11.5 * 1E-9;
    PWRITE(d_mean_primary)
    double N_p_primary = n_dense_packing * pow(d_mean_primary, 3) / pow(d_mean, 3);
    PWRITE(N_p_primary)
    double ks_primary = N_p_primary * 0.75 / 100.0; // [my original]
    if (ks_primary < 1.0) ks_primary = 1.0; // [Dzyan]
    PWRITE(ks_primary)

    double tc = ta0 + sqrt(1.3 / 3.0) * (mu0 / (4.0 * pi)) *
            pow(ks_primary * s_mean, 2) / (kb * pow(d_mean_primary, 3));
    PWRITE(tc)

    qDebug() << "[Hypothesis] Different particles diameter estimation";

    //t_vap_ldc = 419.3663;
    double d_mean_estimate = pow( ((4.0 * pi) / mu0) * kb * (t_vap_ldc - ta0) *
            pow(m_mol_rel / (pi * rho_p * Na * s_mol), 2), 1 / 3.0) * pow(pow(6.0,5) / 2.6, 1 / 6.0) * 1E-2;
    PWRITE(d_mean_estimate)

    qDebug() << "--------------------------------";

    qDebug() << "[Hypothesis] Viscosity of the oleic nanolayer estimation based on cluster evaporation temperature";

    double T_evap_ldc = t_vap_ldc - ta0;
    double Ek = (3 / 2.0) * kb * T_evap_ldc;
    double sigma_sf_nano_coef = Ek / (4 * pi * pow(R_mean, 2) * sigma_sf_oleic);
    PWRITE(sigma_sf_nano_coef)

    qDebug() << "--------------------------------";

    double gamma_tr_mean = 6 * pi * eta_oleic * R_mean;
    double time_damping = m_mean / gamma_tr_mean;
    PWRITE(time_damping)
    double v_heat_mean = sqrt(8 * kb * T_evap_ldc / (pi * m_mean));
    PWRITE(v_heat_mean)
    double F_heat_adaptive = gamma_tr_mean * v_heat_mean;
    double F_oleic_droplet = sigma_sf_nano_coef * sigma_sf_oleic * 2 * pi * R_mean;
    double F_delta_heat = sqrt(2 * kb * T_evap_ldc * gamma_tr_mean / time_damping);
    double forces_relation = F_heat_adaptive / F_oleic_droplet;
    double forces_relation_1 = F_heat_adaptive / F_delta_heat;
    PWRITE(forces_relation)
    PWRITE(forces_relation_1)
    double eta_k = 0.00164; //Pa * s //kerosene
    //R_mean *= 2;s_mean *= 8;m_mean *= 8;
    double Fdd_order = 3 * (mu0 / (4.0 * pi)) * pow(s_mean, 2) / pow(2 * R_mean, 4);
    double gamma_tr_mean_k = 6 * pi * eta_k * R_mean;
    v_heat_mean = sqrt(8 * kb * T_evap_ldc / (pi * m_mean));
    double F_heat_adaptive_k = gamma_tr_mean_k * v_heat_mean;
    double forces_relation_2 = F_heat_adaptive_k / Fdd_order;
    PWRITE(forces_relation_2)

    qDebug() << "--------------------------------";

    qDebug() << "[Validation] Value of T_ch in the GPU research [arXiv:1212.2934v1]";

    double T_ch = (R_mean / M_particle_mean) * sqrt(3 * rho_p / (4 * pi * mu0));
    PWRITE(T_ch)

    qDebug() << "--------------------------------";

    qDebug() << "[Validation] Value of oleic nano surface tension, which stabilizes the osmotic pressure";
    double n_p_oleic = 1.539e+023; // [m-3]
    double p_osm = kb * T_evap_ldc * n_p_oleic;
    double Lx1 = 2.514 * 1E-7; // [m]
    double R_oleic = Lx1 / 8.0;
    double sigma_estimate = p_osm * R_oleic / 2.0;
    PWRITE(sigma_estimate)
    double sigma_sim_empiric = 32.5 * 1E-3 * 5E-2;
    PWRITE(sigma_sim_empiric)
}
