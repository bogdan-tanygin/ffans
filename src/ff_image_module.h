/**************************************************************************
* Copyright (C) 2017 Victor Voshchinskiy<voshchinskiyvitya@ukr.net>
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
#include <Windows.h>
#include <GdiPlus.h>

using namespace Gdiplus;

/*
Represents file types for images (bmp, png, etc.).
*/
enum FileTypes {
	Bmp = 0,
	Png = 1,
	Jpeg = 2,
	Gif = 3,
	Tiff = 4
};

/*
Saves bitmap to file with fileName and fileType (bmp, png, etc.)
Parameters:
- bitmap
- fileName - name of file for saving
- fileType - type of file (bmp, png, etc.)
*/
extern void SaveImage(HBITMAP bitmap, WCHAR* fileName, FileTypes fileType);

/*
Gets encoder for type (bmp, png, etc.)
Parameters:
- format - file tipe (bmp, png, etc.)
- pClsid - out parametr for image codec info
*/
int GetEncoderClsid(const WCHAR* format, CLSID* pClsid);

/*
Get string representation of image file type ("bmp", "png", etc).
Parameters:
- type (see image_module.h FileTypes enum)
*/
WCHAR* GetImageType(FileTypes type);