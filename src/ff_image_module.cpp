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
#include "ff_image_module.h"
#include <Windows.h>
#include <GdiPlus.h>
#include <GdiPlusFlat.h>
#include <stdio.h>

using namespace Gdiplus;
using namespace std;

#pragma comment (lib,"Gdiplus.lib")

// Indicates whether or not should we log file saving status.
bool shouldLog = true;

#pragma region Public methods
/*
Saves bitmap to file with fileName and fileType (bmp, png, etc.)
Parameters: 
 - bitmap
 - fileName - name of file for saving
 - fileType - type of file (bmp, png, etc.)
*/
void SaveImage(HBITMAP hBitmap, WCHAR* fileName, FileTypes fileType)
{
	WCHAR* fType = GetImageType(fileType);
	// Add "image/" to fileType for encoders searching (should be "image/png", "image/bmp", etc).
	WCHAR type[10];
	wcscpy(type, L"image/");
	wcsncat(type, fType, wcslen(fType));

	// Add file type to file name.
	wcsncat(fileName, L".", 1);
	wcsncat(fileName, fType, wcslen(fType));

	GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR gdiplusToken;
	GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

	CLSID encoderClsid;

	// Get the CLSID of the fileType encoder.
	GetEncoderClsid(type, &encoderClsid);

	GpBitmap* bitmap;
	GpStatus stat = DllExports::GdipCreateBitmapFromHBITMAP(hBitmap, NULL, &bitmap);

	if (shouldLog && stat != 0) 
		printf("Failure (image was not saved): stat = %d\n", stat);

	stat = DllExports::GdipSaveImageToFile(bitmap, fileName, &encoderClsid, NULL);
	DllExports::GdipDisposeImage(bitmap);

	if (shouldLog) {
		if (stat == GpStatus::Ok)
			printf("Image was saved successfully\n");
		else
			printf("Failure (image was not saved): stat = %d\n", stat);
	}

	GdiplusShutdown(gdiplusToken);
}
#pragma endregion

#pragma region Private methods
/*
Gets encoder for type (bmp, png, etc.)
Parameters:
 - format - file tipe (bmp, png, etc.)
 - pClsid - out parametr for image codec info
*/
int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
{
	UINT  num = 0;          // number of image encoders
	UINT  size = 0;         // size of the image encoder array in bytes

	ImageCodecInfo* pImageCodecInfo = NULL;

	GetImageEncodersSize(&num, &size);
	if (size == 0)
		return -1;  // Failure

	pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
	if (pImageCodecInfo == NULL)
		return -1;  // Failure

	GetImageEncoders(num, size, pImageCodecInfo);

	// Try to find appropriate encoder.
	for (UINT j = 0; j < num; ++j)
	{
		if (wcscmp(pImageCodecInfo[j].MimeType, format) == 0)
		{
			*pClsid = pImageCodecInfo[j].Clsid;
			free(pImageCodecInfo);
			return j;  // Success
		}
	}

	free(pImageCodecInfo);
	return -1;  // Failure
}

/*
Get string representation of image file type ("bmp", "png", etc).
Parameters: 
 - type (see image_module.h FileTypes enum)
*/
WCHAR* GetImageType(FileTypes type) {
	switch (type) {
		case FileTypes::Bmp: return L"bmp";
		case FileTypes::Png: return L"png";
		case FileTypes::Jpeg: return L"jpeg";
		case FileTypes::Gif: return L"gif";
		case FileTypes::Tiff: return L"tiff";
		default: return L"";
	}
}
#pragma endregion

