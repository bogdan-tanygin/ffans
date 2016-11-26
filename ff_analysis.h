#include <math.h>
#include <string>
using namespace std;
extern double ff_visousity_mix(
								double y1,
								double mu1,
								double y2,
								double mu2);
extern double ff_molar_part (double nu1, double nu2);
extern double ff_mol(double m, double mol_m);
extern void GetScreenShot(string name1);
extern void MultiPosition(float my_x_rot, float my_y_rot, float my_z_off);
extern void ActiveWindow();