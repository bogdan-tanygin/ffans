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

#include "mf_graphics_pal.h"
#include "mf_model_graphics.h"
#include "mf_model.h"

#include <windows.h>
#include <time.h>
#include <GL/gl.h>   // OpenGL.
#include <GL/glu.h>  // GLU support library.
//#include <GL/glut.h>                 // GLUT support library.
#include "glut.h"    // GLUT support library.

#include "mf_io_pal.h"

// Cube position and rotation speed variables.
float x_rot   = 0.0f;
float y_rot   = 0.0f;
float x_speed = 0.0000f;
float y_speed = 0.0000f;
int rotating  = 1;
float z_off   =-2.0f;

int window_id;

// Window and texture IDs, window width and height.
int window_width  = 1500;
int window_height = 1000;

// Frames per second (FPS) statistic variables and routine.
#define FRAME_RATE_SAMPLES 50
int frame_count=0;
float frame_rate=0;


static void do_fps();

void cbRenderScene(void)
{
    //glDisable(GL_TEXTURE_2D);
	glEnable( GL_TEXTURE_2D );

    glEnable(GL_LIGHTING);

    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_DEPTH_TEST); 
    //glDisable(GL_DEPTH_TEST);

    glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST_MIPMAP_NEAREST);
    glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);

   // Need to manipulate the ModelView matrix to move our model around.
    glMatrixMode(GL_MODELVIEW);

    // Reset to 0,0,0; no rotation, no scaling.
    glLoadIdentity(); 

    // Move the object back from the screen.
    glTranslatef(0.0f,0.0f,z_off);

    // Rotate the calculated amount.
    glRotatef(x_rot,1.0f,0.0f,0.0f);
    glRotatef(y_rot,0.0f,0.0f,1.0f);

    // Clear the color and depth buffers.
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    mf_model_next_step();

    // Move back to the origin (for the text, below).
    glLoadIdentity();

    // We need to change the projection matrix for the text rendering.  
    glMatrixMode(GL_PROJECTION);

    // But we like our current view too; so we save it here.
    glPushMatrix();

    // Now we set up a new projection for the text.
    glLoadIdentity();
    glOrtho(0,window_width,0,window_height,-1.0,1.0);

    // Lit or textured text looks awful.
    //glDisable(GL_TEXTURE_2D);
	glEnable( GL_TEXTURE_2D );

    glEnable(GL_LIGHTING);

    // We don't want depth-testing either.
    glDisable(GL_DEPTH_TEST); 

    // But, for fun, let's make the text partially transparent too.
    glColor4f(0.6,1.0,0.6,.75);

    // Now we want to render the calulated FPS at the top.
    
    // To ease, simply translate up.  Note we're working in screen
    // pixels in this projection.
    
    glTranslatef(6.0f,window_height - 14,0.0f);

    mf_mgr_print_info();

    // Done with this special projection matrix.  Throw it away.
    glPopMatrix();

    // All done drawing.  Let's show it.
    glutSwapBuffers();

    // Now let's do the motion calculations.
    if(rotating)
    {
        x_rot += x_speed; 
        y_rot += y_speed; 
    }

    // And collect our statistics.
    do_fps();
}

// ------
// Callback routine executed whenever our window is resized.  Lets us
// request the newly appropriate perspective projection matrix for 
// our needs.  Try removing the gluPerspective() call to see what happens.

void cbResizeScene(int width, int height)
{
    // Let's not core dump, no matter what.
    if (height == 0) height = 1;

    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f,(GLfloat)width/(GLfloat)height,0.1f,100.0f);

    glMatrixMode(GL_MODELVIEW);

    window_width  = width;
    window_height = height;
}

// Does everything needed before losing control to the main
// OpenGL event loop.  
void graph_init(int width, int height) 
{
    // Color to clear color buffer to.
    //glClearColor(0.1f, 0.1f, 0.1f, 0.0f);
	// white layout
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

    // Depth to clear depth buffer to; type of test.
    glClearDepth(1.0);
    //glDepthFunc(GL_LESS); 
	glDepthFunc(GL_LEQUAL);

    // Enables Smooth Color Shading; try GL_FLAT for (lack of) fun.
    glShadeModel(GL_SMOOTH);

    // Load up the correct perspective matrix; using a callback directly.
    cbResizeScene(width,height);

    // A handy trick -- have surface material mirror the color.
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);    
}

static void do_fps() 
{
    static clock_t last=0;
    clock_t now;
    float delta;

    if (++frame_count >= FRAME_RATE_SAMPLES) {
        now  = clock();
        delta= (now - last) / (float) CLOCKS_PER_SEC;
        last = now;

        frame_rate = FRAME_RATE_SAMPLES / delta;
        frame_count = 0;
   }
}

void mf_gr_print(void *font, char *str)
{
    int i,l=strlen(str);
    for(i=0;i<l;i++) glutBitmapCharacter(font,*str++);
}

void mf_gr_init(int argc, char **argv)
{
    x_rot   = 0.0f;
	y_rot   = 0.0f;
	x_speed = 0.0000f;
	y_speed = 0.0000f;
	rotating  = 1;
	z_off   =-2.0f;

    glutInit(&argc, argv);
    SetConsoleTitle( PROGRAM_TITLE );

    // To see OpenGL drawing, take out the GLUT_DOUBLE request.
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(window_width, window_height);

    // Open a window 
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


    // Register the callback function to do the drawing. 
    glutDisplayFunc(&cbRenderScene);

    // If there's nothing to do, draw.
    glutIdleFunc(&cbRenderScene);

    // It's a good idea to know when our window's resized.
    glutReshapeFunc(&cbResizeScene);

    // And let's grab some keyboard input.
    glutKeyboardFunc(&cbKeyPressed);
    glutSpecialFunc(&cbSpecialKeyPressed);

    // Mouse input now also supported
    glutMouseFunc(&cbMouseInput);
    glutMotionFunc(&cbMouseMove); 

    // OK, OpenGL's ready to go.  Let's call our own init function.
    graph_init(window_width, window_height);
}

void mf_gr_loop(void)
{
    // Pass off control to OpenGL.
    // Above functions are called as appropriate.
    glutMainLoop();
}