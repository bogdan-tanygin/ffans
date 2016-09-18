#include "ff_analysis.h"
 double ff_visousity_mix(double y1,
						double mu1,
						double y2,
						double mu2)
{
	return pow(10,((y1*log(mu1))+(y2*log(mu2))));
}
 double ff_molar_part(double nu1, double nu2)
 {
	 return nu1/(nu1+nu2);
 }
 double ff_mol(double m, double mol_m)
 {
	 return m/mol_m;
 }