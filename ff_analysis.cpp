#include "ff_analysis.h"
#include "ff_sys_graphics.h"
#include "ff_model.h"

#include <Windows.h>
#include <GdiPlus.h>
#include <string>
#include <GdiPlusBitmap.h>
#include <Ole2.h>
#include <OleCtl.h>
#include <iostream>
#include <vector>

using namespace Gdiplus;
using namespace std;

HWND hWnd;
bool Active = true;
struct CameraPosition
{
	float x;
	float y;
	float z;
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

 bool saveBitmap(LPCSTR filename, HBITMAP bmp, HPALETTE pal)
{
    bool result = false;
    PICTDESC pd;

    pd.cbSizeofstruct   = sizeof(PICTDESC);
    pd.picType      = PICTYPE_BITMAP;
    pd.bmp.hbitmap  = bmp;
    pd.bmp.hpal     = pal;

    LPPICTURE picture;
    HRESULT res = OleCreatePictureIndirect(&pd, IID_IPicture, false,
                       reinterpret_cast<void**>(&picture));

    if (!SUCCEEDED(res))
    return false;

    LPSTREAM stream;
    res = CreateStreamOnHGlobal(0, true, &stream);

    if (!SUCCEEDED(res))
    {
    picture->Release();
    return false;
    }

    LONG bytes_streamed;
    res = picture->SaveAsFile(stream, true, &bytes_streamed);

    HANDLE file = CreateFile(filename, GENERIC_WRITE, FILE_SHARE_READ, 0,
                 CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, 0);

    if (!SUCCEEDED(res) || !file)
    {
    stream->Release();
    picture->Release();
    return false;
    }

    HGLOBAL mem = 0;
    GetHGlobalFromStream(stream, &mem);
    LPVOID data = GlobalLock(mem);

    DWORD bytes_written;

    result   = !!WriteFile(file, data, bytes_streamed, &bytes_written, 0);
    result  &= (bytes_written == static_cast<DWORD>(bytes_streamed));

    GlobalUnlock(mem);
    CloseHandle(file);

    stream->Release();
    picture->Release();

    return result;
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
	LPCSTR name = name1.c_str();
	bool error = saveBitmap(name,hBitmap,NULL);
	cout<<"Screen Shot ~~~ "<<name<<"  "<<error<<endl;
	
	
	// clean up
    SelectObject(hDC, old_obj);
    DeleteDC(hDC);
    ReleaseDC(NULL, hScreen);
    DeleteObject(hBitmap);
}

void addPosition()
{
	position.push_back(CameraPosition());
	MaxPointOfPosition++;
	position[counterOfPosition].x = x_rot;
	position[counterOfPosition].y = y_rot;
	position[counterOfPosition].z = z_off;
	counterOfPosition++;
}

void delPosition()
{
	if(MaxPointOfPosition)
	{
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
		x_rot = position[counterOfPosition].x;
		y_rot = position[counterOfPosition].y;
		z_off = position[counterOfPosition].z;
	}
	else
	{
		cout<<"U have not saved position yet"<<endl;
	}
}

