/**************************************************************************
 * Copyright (C) 2011,2013 Dr. Bogdan Tanygin<b.m.tanygin@gmail.com>
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

#include <windows.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <gl/glut.h>
#include <gl/glaux.h>

#include "ff_sys_graphics.h"
#include "ff_model_graphics.h"
#include "ff_model.h"
#include "ff_model_parameters.h"

GLUquadric* g_quad;

int show_m, show_b, show_bext;
int show_droplet = 1;
int transp;
int show_sphere;


// vector parameters
double v_len;
//double v_diam;
ff_color4f_t v_col;

void ff_mgr_draw_vector(ff_vect_t r, ff_vect_t v)
{
    double x1,y1,z1;
    double xc,yc,zc; //center of the pyramid basis
    double xe,ye,ze; //top of pyramid
    double dx,dy,dz;
    double v_mod;
    ff_vect_t t, b, r2;
    double t_mod, b_mod;
    double v_diam;

    v_diam = v_len / 5.0;

    v_mod = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);

    t.x = - v.y;
    t.y =   v.x;
    t.z =   0;

    if ((v.x == 0)&&(v.y == 0)) {t.x = 1;}

    t_mod = sqrt(t.x * t.x + t.y * t.y + t.z * t.z);

    x1 = xc = r.x - 0.5 * v_len * v.x / v_mod;
    y1 = yc = r.y - 0.5 * v_len * v.y / v_mod;
    z1 = zc = r.z - 0.5 * v_len * v.z / v_mod;

    xe = r.x + 0.5 * v_len * v.x / v_mod;
    ye = r.y + 0.5 * v_len * v.y / v_mod;
    ze = r.z + 0.5 * v_len * v.z / v_mod;

    dx = 0.5 * v_diam * t.x / t_mod;
    dy = 0.5 * v_diam * t.y / t_mod;
    dz = 0.5 * v_diam * t.z / t_mod;

    x1 += dx;
    y1 += dy;
    z1 += dz;

    b.x = t.y * v.z - t.z * v.y;
    b.y = t.z * v.x - t.x * v.z;
    b.z = t.x * v.y - t.y * v.x;

    b_mod = sqrt(b.x * b.x + b.y * b.y + b.z * b.z);

    r2.x = - t.x * v_diam * 0.25 / t_mod + b.x * v_diam * 0.5 / b_mod;
    r2.y = - t.y * v_diam * 0.25 / t_mod + b.y * v_diam * 0.5 / b_mod;
    r2.z = - t.z * v_diam * 0.25 / t_mod + b.z * v_diam * 0.5 / b_mod;

    glBegin(GL_LINE_STRIP);
    //glBegin(GL_TRIANGLE_STRIP);

    glColor4f( v_col.r,  v_col.g,  v_col.b,  v_col.a);

    glVertex3d(x1, y1, z1);
    glVertex3d(x1 + r2.x, y1 + r2.y, z1 + r2.z);
    glVertex3d(xc + r2.x, yc + r2.y, zc + r2.z);
    glVertex3d(xc - dx, yc - dy, zc - dz);
    glVertex3d(xc - dx - r2.x, yc - dy - r2.y, zc - dz - r2.z);
    glVertex3d(xc - r2.x, yc - r2.y, zc - r2.z);
    glVertex3d(x1, y1, z1);

    glVertex3d(xe, ye, ze);
    glVertex3d(x1 + r2.x, y1 + r2.y, z1 + r2.z);
    glVertex3d(xe, ye, ze);
    glVertex3d(xc + r2.x, yc + r2.y, zc + r2.z);
    glVertex3d(xe, ye, ze);
    glVertex3d(xc - dx, yc - dy, zc - dz);
    glVertex3d(xe, ye, ze);
    glVertex3d(xc - dx - r2.x, yc - dy - r2.y, zc - dz - r2.z);
    glVertex3d(xe, ye, ze);
    glVertex3d(xc - r2.x, yc - r2.y, zc - r2.z);

    glEnd();
}

/*void ff_mgr_set_m_vector(void)
{
double k = 0.2;

v_len = k / 10.0;
//v_diam = k / 50.0;

v_col.r = 0;
v_col.g = 0;
v_col.b = 0;
v_col.a = 0;

glLineWidth(1);
}*/

void ff_mgr_set_bext_vector(void) //should be proportional to the field
{
    v_len  = (50 / 10.0) * fabs(kB);
    //v_diam = (50 / 50.0) * fabs(kB * B0 / Bself);

    /*v_len  = 1 / 10.0;
    v_diam = 1 / 50.0;*/

    v_col.r = 0;
    v_col.g = 1;
    v_col.b = 0;
    v_col.a = 0;

    glLineWidth(1);
}

/*void ff_mgr_set_b_vector(long p) //should be proportional to the field
{
double k = 0.2;

v_len  = (k / 10.0) * sqrt(MUL(B[p],B[p])) / Bself;
//v_diam = (k / 50.0) * sqrt(MUL(B[p],B[p])) / Bself;

v_col.r = 0;
v_col.g = 0;
v_col.b = 1;
v_col.a = 0;

glLineWidth(1);
}*/

/*void ff_mgr_draw_hysteresis() {
glBegin(GL_LINE_STRIP);

glColor3f(1.0, 1.0, 0.0);
glVertex3d(-0.5, -0.5, 0.0);


glEnd();
}*/




void ff_mgr_show_next_step()
{
    ff_vect_t r0, r1, v1;
    double theta, phi;
    double tmag;
    double kvec;
    double space_k = 1E6;
    double kSecondary = 1;
    

    //step++;
    //if (step%5 == 0) printf("\n !!! %e",r[15].x / Lx);

    glBegin(GL_LINE_STRIP);

    glColor3f(space_k * 0.0, space_k * 0.0, space_k * 0.0);
    glVertex3d(-space_k * 0.5 * Lx , -space_k * 0.5 * Ly , -space_k * 0.5 * Lz);
    glVertex3d( space_k * 0.5 * Lx , -space_k * 0.5 * Ly , -space_k * 0.5 * Lz);
    glVertex3d( space_k * 0.5 * Lx ,  space_k * 0.5 * Ly , -space_k * 0.5 * Lz);
    glVertex3d( space_k * 0.5 * Lx ,  space_k * 0.5 * Ly ,  space_k * 0.5 * Lz);
    glVertex3d(-space_k * 0.5 * Lx ,  space_k * 0.5 * Ly ,  space_k * 0.5 * Lz);
    glVertex3d(-space_k * 0.5 * Lx , -space_k * 0.5 * Ly ,  space_k * 0.5 * Lz);
    glVertex3d( space_k * 0.5 * Lx , -space_k * 0.5 * Ly ,  space_k * 0.5 * Lz);
    glVertex3d( space_k * 0.5 * Lx , -space_k * 0.5 * Ly ,  space_k * 0.5 * Lz);
    glVertex3d(-space_k * 0.5 * Lx , -space_k * 0.5 * Ly ,  space_k * 0.5 * Lz);
    glVertex3d(-space_k * 0.5 * Lx ,  space_k * 0.5 * Ly ,  space_k * 0.5 * Lz);
    glVertex3d(-space_k * 0.5 * Lx ,  space_k * 0.5 * Ly , -space_k * 0.5 * Lz);
    glVertex3d(-space_k * 0.5 * Lx , -space_k * 0.5 * Ly , -space_k * 0.5 * Lz);
    glVertex3d( space_k * 0.5 * Lx , -space_k * 0.5 * Ly , -space_k * 0.5 * Lz);
    glVertex3d( space_k * 0.5 * Lx , -space_k * 0.5 * Ly ,  space_k * 0.5 * Lz);
    glVertex3d( space_k * 0.5 * Lx ,  space_k * 0.5 * Ly ,  space_k * 0.5 * Lz);
    glVertex3d( space_k * 0.5 * Lx ,  space_k * 0.5 * Ly , -space_k * 0.5 * Lz);
    glVertex3d(-space_k * 0.5 * Lx ,  space_k * 0.5 * Ly , -space_k * 0.5 * Lz);
    glVertex3d(-space_k * 0.5 * Lx , -space_k * 0.5 * Ly , -space_k * 0.5 * Lz);
    glVertex3d(-space_k * 0.5 * Lx , -space_k * 0.5 * Ly ,  space_k * 0.5 * Lz);

    glEnd();    

    for (long p = 1; p <= pN; p++)
        if (exist_p[p])
        {
            //glEnable(GL_ALPHA_TEST);
            //glEnable(GL_BLEND);
            //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            /*if (transp) glColor4f(0.5, 0.5, 0.5, 0.8);
            else*/ 

            r1.x = space_k * r[p].x;
            r1.y = space_k * r[p].y;
            r1.z = space_k * r[p].z;

            glTranslatef(r1.x, r1.y, r1.z);
			glColor3f(1.0, 0.5, 0.5);
			//glColor3f(0.5, 0.5, 0.5);
			if (show_sphere) gluSphere(g_quad, space_k * Rp[p], 20, 20);
	        glTranslatef(-r1.x, -r1.y, -r1.z);
	
            //glDisable(GL_BLEND);
            //glDisable(GL_ALPHA_TEST);

			if (show_m)
            {
                kvec = 2;
				kSecondary = space_k * Rp[p] * 150;

                //r1.x = space_k * r[p].x / Lx;
                //r1.y = space_k * r[p].y / Ly;
                //r1.z = space_k * r[p].z / Lz;

                tmag = sqrt(m[p].x * m[p].x + m[p].y * m[p].y);

                theta = acos(m[p].z / sqrt(MUL(m[p],m[p])));
                /*if (m[p].x != 0) phi = atan(m[p].y / m[p].x);
                else phi = 0.5 * pi * m[p].y / (fabs(m[p].y) + 0.00000000001);
                if (m[p].x < 0) phi += pi;*/

                glTranslatef(r1.x, r1.y, r1.z);
                glRotatef(theta * 180 / pi, -m[p].y / tmag, m[p].x / tmag, 0);

                glutSolidCone(0.002 * kvec * kSecondary, 0.01 * kSecondary, 20, 20);
                glRotatef(180, 1, 0, 0);
                gluDisk(g_quad, 0, 0.002 * kvec * kSecondary, 20, 20);
                glRotatef(-180, 1, 0, 0);

                glTranslatef(0, 0, -0.01 * kSecondary);
                gluCylinder(g_quad, 0.0005 * kvec * 2 * kSecondary, 0.0005 * kvec * 2 * kSecondary, 0.01 * kSecondary, 20, 20);
                glRotatef(180, 1, 0, 0);
                gluDisk(g_quad, 0, 0.0005 * kvec * 2 * kSecondary, 20, 20);
                glRotatef(-180, 1, 0, 0);
                glTranslatef(0, 0,  0.01 * kSecondary);

                glRotatef(-theta * 180 / pi, -m[p].y / tmag, m[p].x / tmag, 0);
                glTranslatef(-r1.x, -r1.y, -r1.z);

                /*ff_mgr_set_m_vector();
                r1.x = space_k * r[p].x / Lx;
                r1.y = space_k * r[p].y / Ly;
                r1.z = space_k * r[p].z / Lz;
                ff_mgr_draw_vector(r1, m[p]);*/

            }

            /*if (show_b)
            {
            ff_mgr_set_b_vector(p);
            r1.x = space_k * r[p].x / Lx;
            r1.y = space_k * r[p].y / Ly;
            r1.z = space_k * r[p].z / Lz;
            ff_mgr_draw_vector(r1, B[p]);
            }*/
        } // particles loop
        if (show_bext)
        {
            r0.x = - 0.6;
            r0.y = 0;
            r0.z = 0;
            ff_mgr_set_bext_vector();
            v1 = Bext(0,0,0);
            ff_mgr_draw_vector(r0, v1);
        }

		if (show_droplet)
		{
			glEnable(GL_ALPHA_TEST);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			
			glColor4f(0.2, 0.1, 0.6, 0.5);
			gluSphere(g_quad, space_k * R_oleic, 20, 20);

			glDisable(GL_BLEND);
			glDisable(GL_ALPHA_TEST);
		}
}

void ff_mgr_print_info()
{
    char buf[80];
    int shift = -1;

	glDisable(GL_LIGHTING);
	//--------------------

    glColor4f(0,0,0,0.75);
    sprintf(buf,"step = %d", step);

    glRasterPos2i(6,0);
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

    /*glColor4f(0.9,0.2,0.2,0.75);
    sprintf(buf,"FPS: %f F: %2d", FrameRate, FrameCount);

    glRasterPos2i(6,0);
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);*/

    /*glColor4f(0.0,1.0,0.0,0.75);
    sprintf(buf,"Ek = %5.3e J (%5.3e K)", Ek, Ek * 2 / (3.0 * kb));

    glRasterPos2i(6, shift * 20);
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);*/

    glColor4f(0,0,0,0.75);
    sprintf(buf,"dt = %5.3e s", dt);

    glRasterPos2i(6, shift * 20);
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

    ///////////
    glColor4f(0,0,0,0.75);
    sprintf(buf,"t = %5.3e s", t);

    glRasterPos2i(6, shift * 20);
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

    ///////////
    glColor4f(0,0,0,0.75);
    sprintf(buf,"slow_steps = %d", slow_steps);

    glRasterPos2i(6, shift * 20);
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

    /*glColor4f(0,0,0,0.75);
    sprintf(buf,"Hysteresis mode = %d", hyst_mode);

    glRasterPos2i(6, shift * 20);
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);*/

    glRasterPos2i(6, shift * 20);
    shift--;

    ///////////
    /*glColor4f(0.0,1.0,1.0,0.75);
    sprintf(buf,"Average Concentration = %5.3e m^-3", pN / (Lx * Ly * Lz));

    glRasterPos2i(6, shift * 20);
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);*/

    ///////////
    /*glColor4f(0,0,0,0.75);
    sprintf(buf,"Mx = %5.3e A / m", m_tot.x / Vtot);

    glRasterPos2i(6, shift * 20);
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

    glColor4f(0,0,0,0.75);
    sprintf(buf,"My = %5.3e A / m", m_tot.y / Vtot);

    glRasterPos2i(6, shift * 20);
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

    glColor4f(0,0,0,0.75);
    sprintf(buf,"Mz = %5.3e A / m", m_tot.z / Vtot);

    glRasterPos2i(6, shift * 20);
    shift--;
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

	*/

    glColor4f(0,0,0,0.75);
    sprintf(buf,"Hx = %5.3e Oe", BmanX);

    glRasterPos2i(6, shift * 20);
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

    glColor4f(0,0,0,0.75);
    sprintf(buf,"Hy = %5.3e Oe", BmanY);

    glRasterPos2i(6, shift * 20);
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

    glColor4f(0,0,0,0.75);
    sprintf(buf,"Hz =%5.3e Oe", BmanZ);

    glRasterPos2i(6, shift * 20);
    shift--;
    shift--;
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

	glColor4f(0,0,0,0.75);
    sprintf(buf,"Lx = %5.3e m", Lx);

    glRasterPos2i(6, (shift--) * 20);
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

	glColor4f(0,0,0,0.75);
    sprintf(buf,"pseudo_Ek_tr = %5.3e J", Ek_tr);

    glRasterPos2i(6, (shift--) * 20);
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

    glColor4f(0,0,0,0.75);
    sprintf(buf,"pseudo_Ek_rot = %5.3e J", Ek_rot);

    glRasterPos2i(6, (shift--) * 20);
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

	glColor4f(0,0,0,0.75);
    sprintf(buf,"pseudo_dT = %5.3e K", dT);

    glRasterPos2i(6, (shift--) * 20);
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

	glColor4f(0,0,0,0.75);
    sprintf(buf,"T_mean = %5.3e K", T_mean / k_mean);

    glRasterPos2i(6, (shift--) * 20);
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

	glColor4f(0,0,0,0.75);
    sprintf(buf,"T_mean_loc = %5.3e K", T_mean_loc_prev / (k_bm_inst_max - 1));

    glRasterPos2i(6, (shift--) * 20);
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);

	glColor4f(0,0,0,0.75);
    sprintf(buf,"T = %5.3e K", T);

    glRasterPos2i(6, (shift--) * 20);
    ff_gr_print(GLUT_BITMAP_HELVETICA_12,buf);
	
	// ------------------
	glEnable(GL_LIGHTING);
}

void ff_mgr_init()
{
    glPointSize(4);

    transp = 0;
    show_sphere = 1;

    show_m    = 0;
    show_b    = 0;
    show_bext = 0;

    y_rot = 30.0f;
    x_rot = -60.0f;
    z_off = -1.0f;

    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);
    glEnable(GL_CULL_FACE);

    g_quad = gluNewQuadric();
}
