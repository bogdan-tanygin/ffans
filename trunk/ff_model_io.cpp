/**************************************************************************
* Copyright (C) 2011,2013-2014 Dr. Bogdan Tanygin<b.m.tanygin@gmail.com>
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

#include "ff_model.h"
#include "ff_model_io.h"
#include "ff_sys_graphics.h"
#include "ff_model_graphics.h"
#include "ff_model_parameters.h"

double mouse_x, mouse_y;
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
    case GLUT_KEY_PAGE_UP:
        z_off -= 0.05f;
        break;

    case GLUT_KEY_PAGE_DOWN:
        z_off += 0.05f;
        break;

    case GLUT_KEY_UP:
        //x_speed -= 0.01f;
        space_k *= 1.1;
        break;

    case GLUT_KEY_DOWN:
        //x_speed += 0.01f;
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

    case 'D': case 'd':
        Lx *= 2;
        Ly *= 2;
        Lz *= 2;
        break;

    case 'M': case 'm':
        show_m = show_m?0:1;
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

        /*case '[':
        smooth_r *= 0.9;
        m_h_eff_tol *= 0.9;
        break;

        case ']':
        smooth_r *= 1.1;
        m_h_eff_tol *= 1.1;
        break;*/

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

    case 'h':
        hyst_mode++;
        tmp = 0;
        if (hyst_mode == 21) {hyst_mode = 5; tmp = 1;}
        ff_model_upgrade_ext_field();
        if (tmp) {g_hyst_up_line = 1;g_hyst_start_line = 0;}

        //        ff_io_save_hyst();

        mz_glob = 0;
        glob_start_step = step;

        break;

    case '0':
        BmanX = BmanY = BmanZ = 0;
        glob_start_step_susc = 0;
        break;
    case '1':
        if (manual_field_control)
        {
            /*ff_io_save_susceptX(); // save previous point
            m_tot_glob.x = 0;
            m_tot_glob.y = 0;
            m_tot_glob.z = 0;
            glob_start_step_susc = step;*/
            BmanX += 10;
        }
        break;
    case '2':
        if (manual_field_control)
        {
            /*ff_io_save_susceptY();
            m_tot_glob.x = 0;
            m_tot_glob.y = 0;
            m_tot_glob.z = 0;
            glob_start_step_susc = step; */
            BmanY += 10;
        }
        break;
    case '3':
        if (manual_field_control)
        {
            /*ff_io_save_susceptZ();
            m_tot_glob.x = 0;
            m_tot_glob.y = 0;
            m_tot_glob.z = 0;
            glob_start_step_susc = step; */
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
    case ' ':  // F (Space) - Freeze Rotation!
        time_go = time_go?0:1;
        break;

        break; 
    default:
        printf ("KP: No action assigned for %d.\n", key);
        break;
    }
}

/*void ff_io_save_hyst(void)
{
FILE* file;
char s[100];

if (g_hyst_start_line)
{
file = fopen("start_hyst.dat", "a");
B_hyst[hyst_mode] += 0.5 * (g_Bz_prev + Bext(0,0,0).z);
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

B_hyst[hyst_mode] += 0.5 * (g_Bz_prev + Bext(0,0,0).z);
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

B_hyst[hyst_mode] += 0.5 * (g_Bz_prev + Bext(0,0,0).z);
B_hyst_n[hyst_mode] ++;

Mz_hyst[hyst_mode] += mz_glob / (Vtot * (step - glob_start_step));
Mz_hyst_n[hyst_mode] ++;

fprintf(file, "%5.3e %5.3e \n", B_hyst[hyst_mode] / B_hyst_n[hyst_mode], Mz_hyst[hyst_mode] / Mz_hyst_n[hyst_mode]);
fclose(file);
}
}*/

void ff_io_save_setting(ff_vect_t m_tot,double I)
{
    FILE *file, *file1, *file2, *file2_phi, *file2_I, *file2_II, *file2_III;
    double V_oleic = 0;

    file  = fopen("setting_M.dat", "a");
    file1 = fopen("setting_I.dat", "a");
    file2 = fopen("setting_n_agg.dat", "a");
    file2_phi = fopen("setting_phi_agg.dat", "a");
    file2_I = fopen("setting_n_agg_I.dat", "a");
    file2_II = fopen("setting_n_agg_II.dat", "a");
    file2_III = fopen("setting_n_agg_III.dat", "a");

    //fprintf(file,  "%5.3e %5.3e \n", t, sqrt(MUL(m_tot,m_tot)) / m0);
    //fprintf(file1, "%5.3e %5.3e \n", t, I / (R00 * R00));
    fprintf(file,  "%5.3e %5.3e \n", t, sqrt(MUL(m_tot,m_tot)));
    fprintf(file1, "%5.3e %5.3e \n", t, I);

    V_oleic = (4 / 3.0) * pi * pow(R_oleic, 3);
    fprintf(file2, "%5.3e %5.3e \n", t, pN_oleic_drop / V_oleic);
    fprintf(file2_phi, "%5.3e %5.3e \n", t, phi_vol_fract_oleic);
    fprintf(file2_I, "%5.3e %5.3e \n", t, pN_oleic_drop_I / V_oleic);
    fprintf(file2_II, "%5.3e %5.3e \n", t, pN_oleic_drop_II / V_oleic);
    fprintf(file2_III, "%5.3e %5.3e \n", t, pN_oleic_drop_III / V_oleic);

    fclose(file);
    fclose(file1);
    fclose(file2);
    fclose(file2_phi);
    fclose(file2_I);
    fclose(file2_II);
    fclose(file2_III);
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

        fclose(file);
    } //tstep > 0
}

/*void ff_io_save_susceptX(void) // different susceptibility terms from the Bext.x
{
FILE* file;
double Mx, My, Mz; // time-average
double Hx; // external field

Hx = Bext(0,0,0).x / mu0;
Mx = m_tot_glob.x / (Vtot * (step - glob_start_step_susc));
My = m_tot_glob.y / (Vtot * (step - glob_start_step_susc));
Mz = m_tot_glob.z / (Vtot * (step - glob_start_step_susc));

if (BmanX == 0) remove("susceptXX.dat");
file  = fopen("susceptXX.dat", "a");
fprintf(file, "%5.3e %5.3e \n", Hx, Mx);
fclose(file);

if (BmanX == 0) remove("susceptXY.dat");
file  = fopen("susceptXY.dat", "a");
fprintf(file, "%5.3e %5.3e \n", Hx, My);
fclose(file);

if (BmanX == 0) remove("susceptXZ.dat");
file  = fopen("susceptXZ.dat", "a");
fprintf(file, "%5.3e %5.3e \n", Hx, Mz);
fclose(file);
}

void ff_io_save_susceptY(void) // different susceptibility terms from the Bext.x
{
FILE* file;
double Mx, My, Mz; // time-average
double Hy; // external field

Hy = Bext(0,0,0).y / mu0;
Mx = m_tot_glob.x / (Vtot * (step - glob_start_step_susc));
My = m_tot_glob.y / (Vtot * (step - glob_start_step_susc));
Mz = m_tot_glob.z / (Vtot * (step - glob_start_step_susc));

if (BmanY == 0) remove("susceptYX.dat");
file  = fopen("susceptYX.dat", "a");
fprintf(file, "%5.3e %5.3e \n", Hy, Mx);
fclose(file);

if (BmanY == 0) remove("susceptYY.dat");
file  = fopen("susceptYY.dat", "a");
fprintf(file, "%5.3e %5.3e \n", Hy, My);
fclose(file);

if (BmanY == 0) remove("susceptYZ.dat");
file  = fopen("susceptYZ.dat", "a");
fprintf(file, "%5.3e %5.3e \n", Hy, Mz);
fclose(file);
}

void ff_io_save_susceptZ(void) // different susceptibility terms from the Bext.x
{
FILE* file;
double Mx, My, Mz; // time-average
double Hz; // external field

Hz = Bext(0,0,0).z / mu0;
Mx = m_tot_glob.x / (Vtot * (step - glob_start_step_susc));
My = m_tot_glob.y / (Vtot * (step - glob_start_step_susc));
Mz = m_tot_glob.z / (Vtot * (step - glob_start_step_susc));

if (BmanZ == 0) remove("susceptZX.dat");
file  = fopen("susceptZX.dat", "a");
fprintf(file, "%5.3e %5.3e \n", Hz, Mx);
fclose(file);

if (BmanZ == 0) remove("susceptZY.dat");
file  = fopen("susceptZY.dat", "a");
fprintf(file, "%5.3e %5.3e \n", Hz, My);
fclose(file);

if (BmanZ == 0) remove("susceptZZ.dat");
file  = fopen("susceptZZ.dat", "a");
fprintf(file, "%5.3e %5.3e \n", Hz, Mz);
fclose(file);
}*/