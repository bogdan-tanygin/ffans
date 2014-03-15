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

#include "ff_app.h"
#include "ff_sys_graphics.h"
#include "ff_model_graphics.h"
#include "ff_model.h"

#include <windows.h>
#include <time.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include "ff_model_io.h"

int window_id;

float x_rot; 
float y_rot;   
float x_speed; 
float y_speed; 
int rotating;  
float z_off;

float frame_rate;
int frame_count;

int window_width  = 300;
int window_height = 300;

GLdouble nearVal = 1;
GLdouble farVal = 100;

int projection_type = 1; // "0" is an orthographic type; "1" is a perspective one

//static void do_fps(void);

void cbRenderScene(void)
{
    //glDisable(GL_TEXTURE_2D);
    //glEnable( GL_TEXTURE_2D );

    //glEnable(GL_LIGHTING);

    //glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

    //glEnable(GL_DEPTH_TEST); 
    //glDisable(GL_DEPTH_TEST);

    //glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST_MIPMAP_NEAREST);
    //glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);

    //glMatrixMode(GL_MODELVIEW);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity(); 

    glTranslatef(0.0f,0.0f,z_off);

    glRotatef(x_rot,1.0f,0.0f,0.0f);
    glRotatef(y_rot,0.0f,0.0f,1.0f);

    //glTranslatef(0.0f, 0.0f, 0.0f);

    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();
    ff_model_next_step();
    //glPopMatrix();

    //glutSwapBuffers();

    // Output of information to the GL screen
    // --------------------------------------

    glLoadIdentity();

    glMatrixMode(GL_PROJECTION);

    glPushMatrix();

    glLoadIdentity();
    glOrtho(0,window_width,0,window_height,-1.0,1.0);

    glEnable( GL_TEXTURE_2D );

    glEnable(GL_LIGHTING);

    glDisable(GL_DEPTH_TEST);

    glColor4f(0.6,1.0,0.6,.75);

    glTranslatef(0.0f,window_height - 14,0.0f);

    if (show_info) ff_mgr_print_info();

    glPopMatrix();

    glutSwapBuffers();

    if(rotating)
    {
        x_rot += x_speed; 
        y_rot += y_speed; 
    }

    //do_fps();

    glEnable(GL_DEPTH_TEST);

    cbResizeScene(window_width, window_height);
}

void cbResizeScene(int width, int height)
{
    GLdouble clippingMagnitude = 0.5 * sqrt(2.0);
    int side;

    if (height == 0) height = 1;
    if (width < height) side = width;
    else side = height;

    glViewport((width - side) / 2, (height - side) / 2, side, side);
    //glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //gluPerspective(45.0f,(GLfloat)width/(GLfloat)height,0.1f,100.0f);
    if (projection_type == 1) glFrustum(-clippingMagnitude, +clippingMagnitude, -clippingMagnitude, +clippingMagnitude, nearVal, farVal);
    else					  glOrtho  (-clippingMagnitude, +clippingMagnitude, -clippingMagnitude, +clippingMagnitude, nearVal, farVal);

    glMatrixMode(GL_MODELVIEW);

    window_width  = width;
    window_height = height;
}

void graph_init(int width, int height) 
{
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

    glClearDepth(1.0);

    glDepthFunc(GL_LEQUAL);

    glShadeModel(GL_SMOOTH);

    cbResizeScene(width,height);

    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);    
}

/*static void do_fps(void) 
{
static clock_t last=0;
clock_t now;
float delta;
int frame_rate_samples = 50;

if (++(frame_count) >= frame_rate_samples) {
now  = clock();
delta= (now - last) / (float) CLOCKS_PER_SEC;
last = now;

frame_rate = frame_rate_samples / delta;
frame_count = 0;
}
}*/

void ff_gr_print(void *font, char *str)
{
    int i,l=strlen(str);
    for(i=0;i<l;i++) glutBitmapCharacter(font,*str++);
}

void ff_gr_init(int argc, char **argv)
{
    char title0[] = PROGRAM_TITLE;
    x_rot   = 0.0f;
    y_rot   = 0.0f;
    x_speed = 0.0000f;
    y_speed = 0.0000f;
    rotating = 1;

    frame_count = 0;
    frame_rate = 0;

    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(window_width, window_height);

    window_id = glutCreateWindow( PROGRAM_TITLE );

    GLfloat light_ambient[] =
    {0.2, 0.2, 0.2, 1.0};
    GLfloat light_diffuse[] =
    {1.0, 1.0, 1.0, 1.0};
    GLfloat light_specular[] =
    {1.0, 1.0, 1.0, 1.0};
    GLfloat light_position[] =
    {1.0, 1.0, 1.0, 0.0};

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    glEnable(GL_LIGHT0);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);


    glutDisplayFunc(&cbRenderScene);

    glutIdleFunc(&cbRenderScene);

    glutReshapeFunc(&cbResizeScene);

    glutKeyboardFunc(&cbKeyPressed);
    glutSpecialFunc(&cbSpecialKeyPressed);

    glutMouseFunc(&cbMouseInput);
    glutMotionFunc(&cbMouseMove); 

    graph_init(window_width, window_height);
}

void ff_gr_loop(void)
{
    glutMainLoop();
}