/**************************************************************************
* Copyright (C) 2011,2013-2014,2017 Dr. Bogdan Tanygin<b.m.tanygin@gmail.com>
* Copyright (C) 2016,2017 Dmytro Matskevych<dimqqqq@mail.ru>
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
*************************************************************************/

#include <stdio.h>

#include <GL/glut.h>
#include <sstream>

#include "ff_model.h"
#include "ff_model_io.h"
#include "ff_sys_graphics.h"
#include "ff_model_graphics.h"
#include "ff_model_parameters.h"
#include "ff_analysis.h"
double mouse_x, mouse_y;
int mouse_rotate = 0;
int is_ctrl;
ostringstream out;
void cbMouseInput(int button, int state, int x, int y)
{
    switch (button){
    case GLUT_LEFT_BUTTON: 
        if (state == GLUT_DOWN) {mouse_rotate = 1; mouse_x = x; mouse_y = y; is_ctrl = glutGetModifiers();}
        if (state == GLUT_UP) {mouse_rotate = 0;}
        break;
    case GLUT_MIDDLE_BUTTON:
    case GLUT_RIGHT_BUTTON: break;
    default:
        printf ("MouseInput: No action assigned for %d and %d.\n", button, state);
        break;

    }
}

void cbMouseMove(int x, int y)
{
    if (mouse_rotate)
    {
        if (is_ctrl == GLUT_ACTIVE_CTRL)
        {
            z_off += (y - mouse_y)/100.0;
        }
        else
        {
            y_rot += (x - mouse_x)*180/600;
            x_rot += (y - mouse_y)*180/600;
        }
        mouse_x = x;
        mouse_y = y;
    }
}

void cbSpecialKeyPressed(int key, int x, int y)
{
    switch (key) {
    case GLUT_KEY_PAGE_UP:
        z_off -= 0.05f;
        break;

    case GLUT_KEY_PAGE_DOWN:
        z_off += 0.05f;
        break;

    case GLUT_KEY_UP:
        space_k *= 1.1;
        break;

    case GLUT_KEY_DOWN:
        space_k /= 1.1;
        break;

    case GLUT_KEY_LEFT:
        y_speed -= 0.01f;
        break;

    case GLUT_KEY_RIGHT:
        y_speed += 0.01f;
        break;

    default:
        printf ("SKP: No action assigned for %d.\n", key);
        break;
    }
}

void cbKeyPressed(unsigned char key, int x, int y)
{
    int tmp;
    double tmpd;
    switch (key) {
    
    case 'Q': case 'q': case 27:
        glutDestroyWindow(window_id);
        //        exit(1);
        break;

    // force entropy changes
	case 'a':
		ff_io_entropy_change();
		break;
	
	case 'd':
        Lx *= 1.05;
        Ly *= 1.05;
        Lz *= 1.05;
        break;

    case 'D':
        Lx /= 1.05;
        Ly /= 1.05;
        Lz /= 1.05;
        break;

	case 'V': case 'v':
		show_m = show_m ? 0 : 1;
		show_sphere = show_sphere ? 0 : 1;
		break;

	case 'M': case 'm':
        show_m = show_m?0:1;
        break;

    case 'X': case 'x':
        show_cube = show_cube ? 0 : 1;
        break;

    case '5':
        show_steric = show_steric?0:1;
        break;

    case 'O': case 'o':
        show_droplet = show_droplet?0:1;
        break;

    case 'B': case 'b':
        show_b = show_b?0:1;
        break;

    case 'R': case 'r':
        if (BmanY == 0) {BmanY = BmanX; BmanX = 0;}
        else if (BmanX == 0) {BmanX = -BmanY; BmanY = 0;}
        break;

    case 'l':
        ff_io_load(0);
        break;

    case 'L':
        ff_io_load(100);
        break;

    case 'S': case 's':
        show_sphere = show_sphere?0:1;
        break;

    case 'T': case 't':
        T += 1;
        break;

    case 'C': case 'c':
        T -= 1;
        break;

    case 'P': case 'p':
        projection_type = projection_type ? 0 : 1;
        cbResizeScene(window_width, window_height);
        break;

    case 'E': case 'e':
        show_bext = show_bext?0:1;
        break;
	
	case 'N': case 'n':
		ParamInfo();
		out<<"step ="<<step << " V_oleic = " << v_oleic<<" V_car = "<<v_car<<" Bmanz ="<<BmanZ <<".bmp";
		GetScreenShot(out.str());
		break;
	case 'z':
		addPosition(x_rot, y_rot, space_k, projection_type);
		break;
	case 'Z':
		delPosition();
		break;
	case 'f':
		ChangePosition();
		break;

	case '0':
        BmanX = BmanY = BmanZ = 0;
        glob_start_step_susc = 0;
        break;
    case '1':
        if (manual_field_control)
        {
            BmanX += 10;
        }
        break;
    case '2':
        if (manual_field_control)
        {
            BmanY += 10;
        }
        break;
    case '3':
        if (manual_field_control)
        {
            BmanZ += 10;
        }
        break;
    case 'i':case 'I':
        show_info = show_info?0:1;
        break;
    case '9':
        Lx /= 2;
        Ly /= 2;
        Lz /= 2;
        break;
    case ' ':
        time_go = time_go?0:1;
		ActiveWindow();
        break;
    
    case '[':
        scaling_cube /= 1.05;
        break;
    
    case ']':
        scaling_cube *= 1.05;
        break;

    case ',':
        gr_quality --;
        break;
    
    case '.':
        gr_quality ++;
        break;

    case 'g':
        gr_x0 -= 0.05;
        break;
    
    case 'h':
        gr_x0 += 0.05;
        break;
    
    case 'j':
        gr_y0 -= 0.05;
        break;
    
    case 'k':
        gr_y0 += 0.05;
        break;

    case ';':
        gr_z0 -= 0.05;
        break;
    
    case '\'':
        gr_z0 += 0.05;
        break;

    case 'u':
        dt *= 1.1;
        break;
    
    case 'y':
        dt /= 1.1;
        break;

    default:
        printf ("KP: No action assigned for %d.\n", key);
        break;
    }
}

void ff_io_save_setting(ff_vect_t m_tot,double I)
{
    FILE *file, *file1;
    double V_oleic = 0;

    file  = fopen("setting_M.dat", "a");
    file1 = fopen("setting_I.dat", "a");

    fprintf(file,  "%5.3e %5.3e \n", t, sqrt(MUL(m_tot,m_tot)));
    fprintf(file1, "%5.3e %5.3e \n", t, I);

    fclose(file);
    fclose(file1);
}

void ff_io_autosave(void)
{
    FILE* file;
    char str[50];

    if (auto_save)
    {

        sprintf(str, "%d.dat", step);

        file  = fopen(str, "w");

        fprintf(file,  "%5.3e ", t);

        for(long p = 1; p <= pN; p++)
        {
            fprintf(file, "p = %d ", p);

            fprintf(file,  "%5.3e ", Rp[p]);

            fprintf(file,  "%5.3e ", r[p].x);
            fprintf(file,  "%5.3e ", r[p].y);
            fprintf(file,  "%5.3e ", r[p].z);

            fprintf(file,  "%5.3e ", m[p].x);
            fprintf(file,  "%5.3e ", m[p].y);
            fprintf(file,  "%5.3e ", m[p].z);

            fprintf(file,  "%5.3e ", v[p].x);
            fprintf(file,  "%5.3e ", v[p].y);
            fprintf(file,  "%5.3e ", v[p].z);

            fprintf(file,  "%5.3e ", w[p].x);
            fprintf(file,  "%5.3e ", w[p].y);
            fprintf(file,  "%5.3e ", w[p].z);
        }

        fclose(file);
    } // end of if (auto_save)
}

void ff_io_load(long tstep)
{
    FILE* file;
    float tmp;
    char str[50];
    long p;

    glob_start_step = step;
    glob_start_step_susc = step;
    //k_bm_inst = 1;

    if (tstep == 0)
    {
        printf("Step to load (must be > 0) = ");
        scanf("%d", &tstep);
    }

    if (tstep > 0)
    {
        sprintf(str, "%d.dat", tstep);

        file  = fopen(str, "r");

        fscanf(file, "%f", &tmp);
        t = tmp;

        do
        {
            fscanf(file,"p = %d", &p);
            // -------------------

            fscanf(file, "%f", &tmp);
            Rp[p] = tmp;

            fscanf(file, "%f", &tmp);
            r[p].x = tmp;
            fscanf(file, "%f", &tmp);
            r[p].y = tmp;
            fscanf(file, "%f", &tmp);
            r[p].z = tmp;

            fscanf(file, "%f", &tmp);
            m[p].x = tmp;
            fscanf(file, "%f", &tmp);
            m[p].y = tmp;
            fscanf(file, "%f", &tmp);
            m[p].z = tmp;

            fscanf(file, "%f", &tmp);
            v[p].x = tmp;
            fscanf(file, "%f", &tmp);
            v[p].y = tmp;
            fscanf(file, "%f", &tmp);
            v[p].z = tmp;

            fscanf(file, "%f", &tmp);
            w[p].x = tmp;
            fscanf(file, "%f", &tmp);
            w[p].y = tmp;
            fscanf(file, "%f ", &tmp);
            w[p].z = tmp;

            Rp0[p] = Rp[p] - delta;
            ff_model_size_dispersion_param_calc(Rp0[p], p);
        }
        while(!feof(file));
		step = tstep;
        fclose(file);
    } //tstep > 0
}

void ff_io_entropy_change(void)
{
	double k_r = 1.01;
	long p, ps;
	double G = 0; // total energy per particles
	double S = 0; // entropy S = S - S0, where S0 is arbitrary value
	double V, a_free;
	double F = 0; // free energy
	FILE* file;

	// total energy calc
	G = 0;
	for (p = 1; p <= pN - 1; p++)
	{
		G += G_dd[p] / 2.0 + ff_model_G_zeeman(p);

		for (ps = p + 1; ps <= pN; ps++)
		{
			G += ff_model_G_london(p, ps) + ff_model_G_steric(p, ps);
		}
	}
	G /= pN;

	// entropy calculation (only for monodisperse ferrofluid!)
	a_free = pow(pi * 1.35 * pow(2 * Rp0[1], 3) / (6 * phi_vol_fract / 100), 1 / 3.0);
	V = (4 * pi / 3.0) * pow(a_free - 2 * Rp0[1], 3);
	S = kb * log(V);

	// free energy
	F = G - T * S;

	//printf("\n G = %e", G);
	printf("\n F / (k * T) = %e", F / (kb * T));

	file  = fopen("F(phi).dat", "a");
	fprintf(file, "%5.3e %5.3e \n", phi_vol_fract, F / (kb * T));
	fclose(file);

	for (p = 1; p <= pN - 1; p++)
	{
		r[p].x *= k_r;
		r[p].y *= k_r;
		r[p].z *= k_r;
	}

	Lx *= k_r;
	Ly *= k_r;
	Lz *= k_r;
}