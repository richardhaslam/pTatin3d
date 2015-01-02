/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    pTatin3d
 **    filename:   main.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
/*
 Protoytpe for material properties
 
 
 1/ python mat_prop_class_generator.py
 2/ Edit files generated, remove #error() macro from each file
 3/ gcc -g -O0 -c MaterialConst_ViscosityArrhenius_def.c
 4/ gcc -g -O0 -o main.app main.c MaterialConst_ViscosityArrhenius_def.o
 5/ ./main.app
 
 */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "MaterialConst_ViscosityArrhenius_def.h"

void MatPropProcessOptions(const char   prefix[],
													 const int    regionid,
													 void         *dataarray,
													 const size_t class_size,
													 const char   classname[],
													 const int    nmembers,
													 const char   *membernames[],
													 const size_t member_sizes[],
													 const size_t member_offsets[] )
{
	void *data;
	int f;
	void *offset;
	double *mydata;
	
	printf("Region: %d \n", regionid );
	printf("Class: %s \n", classname );

	data = (void*)( (char*)dataarray + regionid * class_size );
	
	for (f=0; f<nmembers; f++) {
		//printf("member_offsets = %zd \n", member_offsets[f] );
		offset = (void*)( (char*)data + member_offsets[f] );
		
		if (prefix) {
			printf("option: -%s_%s_%s_%d\n",prefix,classname,membernames[f],regionid);
		} else {
			printf("option: -%s_%s_%d\n",classname,membernames[f],regionid);
		}
		mydata = (double*)offset;
		//printf("mydata %p \n", mydata);
		printf("  data value = %1.4e \n",*mydata);
	}
	
}


int main(void)
{
	int f;
	MaterialConst_ViscosityArrhenius data[5];
	
	
	for (f=0; f<5; f++) {
		MaterialConst_ViscosityArrhenius *me;
		
		data[f].enthalpy = 10.0 + (double)f;
		data[f].preexpA = (double)(20.0 + (double)f);
		data[f].nexp = 30.0 + (double)f;
		data[f].Vmol = 40.0 + (double)f;
		data[f].Tref = 50.0 + (double)f;
		data[f].Ascale = 60.0 + (double)f;

		me = &data[f];
		printf("&data[f] = %p, (%p %p %p %p %p) \n", me,((char*)me+0),((char*)me+8),((char*)me+16),((char*)me+24),((char*)me+32) );
		printf("e,A,n,V,T %p %p %p %p %p \n", &data[f].enthalpy,&data[f].preexpA,&data[f].nexp,&data[f].Vmol,&data[f].Tref );
	}
	
	
	MatPropProcessOptions(		"GENE",
											 3,
											 (void*)data,
											 sizeof(MaterialConst_ViscosityArrhenius),
											 MaterialConst_ViscosityArrhenius_classname_short,
											 MaterialConst_ViscosityArrhenius_nmembers,
											 MaterialConst_ViscosityArrhenius_member_names_short, /* or MaterialConst_ViscosityArrhenius_member_names for more descriptive option names */
											 MaterialConst_ViscosityArrhenius_member_sizes,
											 MaterialConst_ViscosityArrhenius_member_byte_offset);
											 
	
	return 0;
}