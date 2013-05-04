/**************************************************************************
 * Copyright (C) 2011 Dr. Bogdan Tanygin<b.m.tanygin@gmail.com>
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
//#include <GL/gl.h>   // OpenGL.
//#include <GL/glu.h>  // GLU support library.
//#include <GL/glut.h>                 // GLUT support library.
#include "glut.h"    // GLUT support library.

#include "mf_model.h"
#include "mf_io_pal.h"
#include "mf_graphics_pal.h"
#include "mf_model_graphics.h"
#include "mf_model_parameters.h"

// var. for mouse control
long double mouse_x, mouse_y;
int mouse_rotate = 0;
int is_ctrl;

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
   case GLUT_KEY_PAGE_UP: // move the cube into the distance.
      z_off -= 0.05f;
      break;

   case GLUT_KEY_PAGE_DOWN: // move the cube closer.
      z_off += 0.05f;
      break;

   case GLUT_KEY_UP: // decrease x rotation speed;
      x_speed -= 0.01f;
      break;

   case GLUT_KEY_DOWN: // increase x rotation speed;
      x_speed += 0.01f;
      break;

   case GLUT_KEY_LEFT: // decrease y rotation speed;
      y_speed -= 0.01f;
      break;

   case GLUT_KEY_RIGHT: // increase y rotation speed;
      y_speed += 0.01f;
      break;

   default:
      printf ("SKP: No action assigned for %d.\n", key);
      break;
    }
}

///  keysetting
void cbKeyPressed(unsigned char key, int x, int y)
{
	int tmp;
    switch (key) {
    case 'Q': case 'q': case 27: // Q (Escape) - We're outta here.
        glutDestroyWindow(window_id);
        exit(1);
        break; // exit doesn't return, but anyway...

    case 'M': case 'm':
        show_m = show_m?0:1;
        break;
        
    case 'B': case 'b':
        show_b = show_b?0:1;
        break;

	case '[':
        smooth_r *= 0.9;
		m_h_eff_tol *= 0.9;
        break;

	case ']':
        smooth_r *= 1.1;
		m_h_eff_tol *= 1.1;
        break;

	case 'L': case 'l':
        mf_io_load();
        break;

	case 'T': case 't':
        transp = transp?0:1;
        break;

	case 'S': case 's':
        show_sphere = show_sphere?0:1;
        break;

    case 'E': case 'e':
        show_bext = show_bext?0:1;
        break;

	case 'h':
        hyst_mode++;
		tmp = 0;
		if (hyst_mode == 21) {hyst_mode = 5; tmp = 1;}
		mf_model_upgrade_ext_field();
		if (tmp) {g_hyst_up_line = 1;g_hyst_start_line = 0;}

		mf_io_save_hyst();

		mz_glob = 0;
		glob_start_step = step;

		break;
        
    case ' ':  // F (Space) - Freeze Rotation!
        time_go = time_go?0:1;
        break;
       
        break; 
    default:
        printf ("KP: No action assigned for %d.\n", key);
        break;
    }
}

void mf_io_save_hyst(void)
{
	FILE* file;
	char s[100];

	if (g_hyst_start_line)
	{
		file = fopen("start_hyst.dat", "a");
		B_hyst[hyst_mode] += 0.5 * (g_Bz_prev + Bext().z);
		B_hyst_n[hyst_mode] ++;

		Mz_hyst[hyst_mode] += mz_glob / (Vtot * (step - glob_start_step));
		Mz_hyst_n[hyst_mode] ++;

		fprintf(file, "%5.3e %5.3e \n", B_hyst[hyst_mode] / B_hyst_n[hyst_mode], Mz_hyst[hyst_mode] / Mz_hyst_n[hyst_mode]);
		fclose(file);
	}

	if (g_hyst_up_line)
	{
		if (hyst_mode == 5)
		{
			remove("up_hyst_old.dat");
			rename("up_hyst.dat", "up_hyst_old.dat");
			file = fopen("up_hyst.dat", "w");
		}
		else file = fopen("up_hyst.dat", "a");

		B_hyst[hyst_mode] += 0.5 * (g_Bz_prev + Bext().z);
		B_hyst_n[hyst_mode] ++;

		Mz_hyst[hyst_mode] += mz_glob / (Vtot * (step - glob_start_step));
		Mz_hyst_n[hyst_mode] ++;

		fprintf(file, "%5.3e %5.3e \n", B_hyst[hyst_mode] / B_hyst_n[hyst_mode], Mz_hyst[hyst_mode] / Mz_hyst_n[hyst_mode]);
		fclose(file);
	}

	if (g_hyst_bottom_line)
	{
		if (hyst_mode == 13)
		{
			remove("bottom_hyst_old.dat");
			rename("bottom_hyst.dat", "bottom_hyst_old.dat");
			file = fopen("bottom_hyst.dat", "w");
		}
		else file = fopen("bottom_hyst.dat", "a");

		B_hyst[hyst_mode] += 0.5 * (g_Bz_prev + Bext().z);
		B_hyst_n[hyst_mode] ++;

		Mz_hyst[hyst_mode] += mz_glob / (Vtot * (step - glob_start_step));
		Mz_hyst_n[hyst_mode] ++;

		fprintf(file, "%5.3e %5.3e \n", B_hyst[hyst_mode] / B_hyst_n[hyst_mode], Mz_hyst[hyst_mode] / Mz_hyst_n[hyst_mode]);
		fclose(file);
	}
}

void mf_io_save_setting(mf_vect_t m_tot,long double I)
{
	FILE* file, *file1;

	file  = fopen("setting_M.dat", "a");
	file1 = fopen("setting_I.dat", "a");

	fprintf(file,  "%5.3e %5.3e \n", t, sqrt(MUL(m_tot,m_tot)) / m0);
	fprintf(file1, "%5.3e %5.3e \n", t, I / (R00 * R00));

	fclose(file);
	fclose(file1);

}

void mf_io_autosave(void)
{
	FILE* file;
	char str[50];

	sprintf(str, "%d.dat", step);

	file  = fopen(str, "a");

	fprintf(file,  "%5.3e ", t);
	
	for(long p = 1; p <= pN; p++)
	{
		fprintf(file, "%d ", p);
		
		fprintf(file,  "%5.3e ", r[p].x);
		fprintf(file,  "%5.3e ", r[p].y);
		fprintf(file,  "%5.3e ", r[p].z);

		fprintf(file,  "%5.3e ", v[p].x);
		fprintf(file,  "%5.3e ", v[p].y);
		fprintf(file,  "%5.3e ", v[p].z);

		fprintf(file,  "%5.3e ", m[p].x);
		fprintf(file,  "%5.3e ", m[p].y);
		fprintf(file,  "%5.3e ", m[p].z);
	}

	fclose(file);
}

void mf_io_load(void)
{
	FILE* file;
	long tstep, tstep1;
	float tmp;
	char str[50];
	long p;

	printf("Step = ");
	scanf("%d", &tstep);
	sprintf(str, "%d.dat", tstep);

	file  = fopen(str, "r");

	fscanf(file, "%f", &tmp);
	t = tmp;

	do
	{
		fscanf(file,"%d", &p);

		fscanf(file, "%f", &tmp);
		r[p].x = tmp;
		fscanf(file, "%f", &tmp);
		r[p].y = tmp;
		fscanf(file, "%f", &tmp);
		r[p].z = tmp;

		fscanf(file, "%f", &tmp);
		v[p].x = tmp;
		fscanf(file, "%f", &tmp);
		v[p].y = tmp;
		fscanf(file, "%f", &tmp);
		v[p].z = tmp;

		fscanf(file, "%f", &tmp);
		m[p].x = tmp;
		fscanf(file, "%f", &tmp);
		m[p].y = tmp;
		fscanf(file, "%f", &tmp);
		m[p].z = tmp;
	}
	while(!feof(file));

	fclose(file);
}