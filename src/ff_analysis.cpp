/**************************************************************************
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
#include "ff_analysis.h"
#include "ff_sys_graphics.h"
#include "ff_model.h"
#include "ff_model_graphics.h"
#include "ff_iniParam.h"
#include "ff_image_module.h"

#include <Windows.h>
#include <GdiPlus.h>
#include <string>
#include <GdiPlusBitmap.h>
#include <Ole2.h>
#include <OleCtl.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdio.h>


using namespace Gdiplus;
using namespace std;

HWND hWnd;
bool Active = true;
struct CameraPosition
{
	float x;
	float y;
	float z;
	int projection_type_status;
};
int counterOfPosition = 0;
int MaxPointOfPosition =0;
vector<CameraPosition> position;
void ActiveWindow()
{
	if(Active){hWnd = GetForegroundWindow();Active = false;}
}
 double ff_visousity_mix(double y1,
						double mu1,
						double y2,
						double mu2)
{
	return pow(10,((y1*log10(mu1))+(y2*log10(mu2))));
}
 double ff_molar_part(double nu1, double nu2)
 {
	 return nu1/(nu1+nu2);
 }
 double ff_mol(double m, double mol_m)
 {
	 return m/mol_m;
 }

void GetScreenShot(string name1) //Make Screen Shot
{
    int x1, y1, x2, y2, w, h;
    // get screen dimensions
    x1  = GetSystemMetrics(SM_XVIRTUALSCREEN);
    y1  = GetSystemMetrics(SM_YVIRTUALSCREEN);
    x2  = GetSystemMetrics(SM_CXVIRTUALSCREEN);
    y2  = GetSystemMetrics(SM_CYVIRTUALSCREEN);
    w   = x2 - x1;
    h   = y2 - y1;

    // copy screen to bitmap 
	HDC     hScreen = GetDC(hWnd);
    HDC     hDC     = CreateCompatibleDC(hScreen);
    HBITMAP hBitmap = CreateCompatibleBitmap(hScreen, w, h);
    HGDIOBJ old_obj = SelectObject(hDC, hBitmap);
    BOOL    bRet    = BitBlt(hDC, 0, 0, w, h, hScreen, x1, y1, SRCCOPY);

    // save bitmap to clipboard
    OpenClipboard(NULL);
    EmptyClipboard ();
    SetClipboardData(CF_BITMAP, hBitmap);
    CloseClipboard();   


	//CImage image;
	//image.Attach(hBitmap);
	//image.Save("Test",ImageFormatBMP);
	//CloseClipboard();

	//var DataSet = 
	
	HBITMAP hb = (HBITMAP)SelectObject(hDC, old_obj);

	// Convert string to WCHAR*.
	const size_t cSize = strlen(name1.data()) + 1;
	wchar_t* name = new wchar_t[cSize];
	mbstowcs(name, name1.data(), cSize);

	int type = (int)iniGet("ImageSettings", "FileType");

	SaveImage(hBitmap, name, (FileTypes)type);
	
	// clean up
    SelectObject(hDC, old_obj);
    DeleteDC(hDC);
    ReleaseDC(NULL, hScreen);
    DeleteObject(hBitmap);
}

void addPosition(float x, float y, float z, int pts)
{
	position.push_back(CameraPosition());
	MaxPointOfPosition++;
	counterOfPosition = MaxPointOfPosition-1;
	position[counterOfPosition].x = x;
	position[counterOfPosition].y = y;
	position[counterOfPosition].z = z;
	position[counterOfPosition].projection_type_status = pts;
    if(isShowInfo!=0)
    {
	    cout<<"Position["<<counterOfPosition<<"] ("<<position[counterOfPosition].x<<","<<position[counterOfPosition].y<<","<<position[counterOfPosition].z<<") is saved"<<endl;
    }
    counterOfPosition++;
}

void delPosition()
{
	if(MaxPointOfPosition)
	{
        if(isShowInfo!=0)
        {
		    cout<<"Position["<<counterOfPosition<<"] ("<<position[counterOfPosition].x<<","<<position[counterOfPosition].y<<","<<position[counterOfPosition].z<<") is removed"<<endl;
        }
        position.erase(position.begin()+counterOfPosition);
		MaxPointOfPosition--;
		if(counterOfPosition<MaxPointOfPosition-1)
		{
			counterOfPosition++;
		}
		else
		{
			counterOfPosition=0;
		}
	}
	else
	{
		cout<<"U have not saved position yet"<<endl;
	}
}

void ChangePosition()
{
	if(MaxPointOfPosition)
	{
		if(counterOfPosition<MaxPointOfPosition-1)
		{
			counterOfPosition++;
		}
		else
		{
			counterOfPosition=0;
		}
        if(isShowInfo!=0)
        {
		    cout<<"Position["<<counterOfPosition<<"] ("<<position[counterOfPosition].x<<","<<position[counterOfPosition].y<<","<<position[counterOfPosition].z<<") is changed"<<endl;
        }
        projection_type = position[counterOfPosition].projection_type_status;
		x_rot = position[counterOfPosition].x;
		y_rot = position[counterOfPosition].y;
		space_k = position[counterOfPosition].z;
		cbResizeScene(window_width, window_height);
		
	}
	else
	{
		cout<<"U have not saved position yet"<<endl;
	}
}

void auto_set_position(int isAutoSet_)
{
    int iMax;
    if(isAutoSet_ != 0)
    {
    iMax=(int)iniGet("SetPosition","PositionMax");
        for(int i = 1;i<=iMax;i++)
        {
        ostringstream x;
        ostringstream y;
        ostringstream z;
        ostringstream pts;
        x<<"x"<<i;
        y<<"y"<<i;
        z<<"z"<<i;
        pts<<"pts"<<i;
        addPosition(iniGet("SetPosition",x.str()),iniGet("SetPosition",y.str()),iniGet("SetPosition",z.str()),iniGet("SetPosition",pts.str()));
        }
    }
}
