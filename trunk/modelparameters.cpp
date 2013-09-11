#include "modelparameters.h"

// Fundamental constants [Metric System]
// -------------------------------------

// 1. Math and geometry
const double pi = acos(-1.0); // []

// 2. Physics

// 2.1. Electrodynamics

// Vacuum permeability
const double mu0 = 4 * pi * 1E-7; // [V·s/(A·m)]

// 2.2. Thermodynamics

// Avogadro constant
const double Na = 6.02214129 * 1E23; // [mol−1]

// Boltzmann constant
const double kb = 1.3806488 * 1E-23; // [m2 kg s-2 K-1]

// Gas constant
const double R = Na * kb; // [J K−1 mol−1]

ModelParameters::ModelParameters()
{
}
