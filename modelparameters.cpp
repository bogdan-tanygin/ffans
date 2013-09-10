#include "modelparameters.h"

// 1. Fundamental constants [Metric System]
// -------------------------------------

// 1.1. Math and geometry
double pi = acos(-1.0); // []

// 1.2. Physics

// 1.2.1. Electrodynamics

// Vacuum permeability
double mu0 = 4 * pi * 1E-7; // [V·s/(A·m)]

// 1.2.2. Thermodynamics

// Avogadro constant
double Na = 6.02214129 * 1E23; // [mol−1]

// Boltzmann constant
double kb = 1.3806488 * 1E-23; // [m2 kg s-2 K-1]

// Gas constant
double R = Na * kb; // [J K−1 mol−1]

ModelParameters::ModelParameters()
{
}
