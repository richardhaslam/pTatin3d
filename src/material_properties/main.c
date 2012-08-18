/*
 Protoytpe for material properties
 
 
 1/ python mat_prop_class_generator.py
 2/ Edit files generated, remove #error() macro from each file
 3/ gcc -g -O0 -o main.app main.c MatProp_RheologyArrhenius_def.o
 4/ gcc -g -O0 -c MatProp_RheologyArrhenius_def.c
 5/ ./main.app
 
 */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "MatProp_RheologyArrhenius_def.h"

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
	MatProp_RheologyArrhenius data[5];
	
	
	for (f=0; f<5; f++) {
		MatProp_RheologyArrhenius *me;
		
		data[f].b = 10.0 + (double)f;
		data[f].n = (double)(20.0 + (double)f);
		data[f].E = 30.0 + (double)f;
		data[f].V = 40.0 + (double)f;
		data[f].R = 50.0 + (double)f;

		me = &data[f];
		printf("&data[f] = %p, (%p %p %p %p %p) \n", me,((char*)me+0),((char*)me+8),((char*)me+16),((char*)me+24),((char*)me+32) );
		printf("b,n,E,V,R %p %p %p %p %p \n", &data[f].b,&data[f].n,&data[f].E,&data[f].V,&data[f].R );
	}
	
	
	MatPropProcessOptions(		"GENE",
											 3,
											 (void*)data,
											 sizeof(MatProp_RheologyArrhenius),
											 MatProp_RheologyArrhenius_classname_short,
											 MatProp_RheologyArrhenius_nmembers,
											 MatProp_RheologyArrhenius_member_names_short, /* or MatProp_RheologyArrhenius_member_names for more descriptive option names */
											 MatProp_RheologyArrhenius_member_sizes,
											 MatProp_RheologyArrhenius_member_byte_offset);
											 
	
	return 0;
}