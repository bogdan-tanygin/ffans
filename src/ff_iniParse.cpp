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
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "ff_iniParam.h"
#include "ff_model_parameters.h"
using namespace std;
struct iniParam
{
	string fileName;
	string section;
	string name;
	double value;
};
vector <iniParam> Param;
void wordsInStruct(FILE *f1)
{
	int counter = 0;
	char ch;
	float value;
	string mySection="NoSection";
	string name = "noName";
	string svalue;
	string::size_type sz;
	while((ch=fgetc(f1))!=EOF)
	{
		if(ch==';'||ch=='#')
		{
			while((ch=fgetc(f1))!='\n')
			{
			}
			continue;
		}
		if(ch=='[')
		{
			mySection="";
			while((ch=fgetc(f1))!='\n')
			{
				if(ch!=']')
				{
					mySection+=ch;
				}
				else
				{
					break;
				}
			}
			continue;
		}
		if(ch>='A'&&ch<='z')
		{
			name=ch;
			while((ch=fgetc(f1))!='=')
			{
				if(ch!=' ')
				{
				name+=ch;
				}
			}
			value = 0;
			while((ch=fgetc(f1))==' '){}
			if((ch>='0'&&ch<='9')||ch=='-')
			{
				svalue = ch;
				while((ch=fgetc(f1))=='.'||(ch>='0'&&ch<='9')||ch=='e'||ch=='E'||ch=='+'||ch=='-')
				{
					svalue += ch;
				}		
				value = stof(svalue,&sz);
				Param.push_back(iniParam());
				Param[counter].fileName = "";
				Param[counter].name = name;
				Param[counter].section = mySection;
				Param[counter].value = value;
				counter++;
			}
		}
	}
}
void info()
{
	for(int i =0;i<Param.size();i++)
	{
		cout<<"Section ="<<Param[i].section<<"name ="<<Param[i].name<<"value ="<<Param[i].value<<endl;
	}
}
double iniGet(string section, string name)
{
	for(int i=0;i<Param.size();i++)
	{
		if(Param[i].section==section&&Param[i].name==name)
		{
            if(isShowInfo!=0)
            {
			  cout<<section<<"->"<<name<<"="<<Param[i].value<<"---Sucsses"<<endl;
            }
			return Param[i].value;
		}
	}
	cout<<section<<"->"<<name<<"---Fail"<<endl;
	return 0;
}
int ff_paramAnalysis( ) 
{
	string fileName = "ff_model_parameters.ini";
	FILE *f1 =  fopen ("samples/ff_model_parameters.ini", "r"); 
	wordsInStruct(f1);
//	info();
	cout<<"Import from:"<<fileName<<endl;
	fclose(f1);
	return 1;
}
