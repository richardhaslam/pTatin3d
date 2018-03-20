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
 **    filename:   phase_map.c
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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "ptatin3d.h"
#include "ptatin3d_defs.h"
#include "phase_map.h"

void PhaseMapCreate(PhaseMap *map)
{
	PhaseMap pm;
	pm = malloc(sizeof(struct _p_PhaseMap));
	memset(pm,0,sizeof(struct _p_PhaseMap));
	*map = pm;
}

void PhaseMapDestroy(PhaseMap *map)
{
	PhaseMap pm;

	if (map==NULL) { return; }
	pm = *map;

	if (pm->data!=NULL) {
		free(pm->data);
		pm->data = NULL;
	}
	*map = NULL;
}

void PhaseMapGetIndex(PhaseMap pm,const int i,const int j, int *index)
{
	if (i<0) { printf("ERROR(%s): i = %d  <0 \n", __func__, i ); exit(EXIT_FAILURE); }
	if (j<0) { printf("ERROR(%s): j = %d < 0 \n", __func__, j ); exit(EXIT_FAILURE); }
	if (i>=pm->mx) { printf("ERROR(%s): i = %d > %d\n", __func__, i, pm->mx ); exit(EXIT_FAILURE); }
	if (j>=pm->my) { printf("ERROR(%s): j = %d > %d\n", __func__, j, pm->my ); exit(EXIT_FAILURE); }


	*index = i + j * pm->mx;
}

void PhaseMapLoadFromFile_ASCII(const char filename[],PhaseMap *map)
{
	FILE *fp = NULL;
	PhaseMap phasemap;
  char dummy[1000];

	int i,j,phasemap_max;
  int index;

	/* open file to parse */
	fp = fopen(filename,"r");
	if (fp==NULL) {
		printf("Error(%s): Could not open file: %s \n",__func__, filename );
		exit(EXIT_FAILURE);
	}

	/* create data structure */
	PhaseMapCreate(&phasemap);

	/* read header information, mx,my,x0,y0,x1,y1 */
  //  fscanf(fp,"%s\n",dummy);
  if (!fgets(dummy,sizeof(dummy),fp)) {printf("fgets() failed. Exiting ungracefully.\n");exit(1);}
    if (fscanf(fp,"%d\n",&phasemap->mx) < 1) {printf("fscanf() failed. Exiting ungracefully.\n");exit(1);}
    if (fscanf(fp,"%d\n",&phasemap->my) < 1) {printf("fscanf() failed. Exiting ungracefully.\n");exit(1);}
    if (fscanf(fp,"%lf %lf %lf %lf\n",&phasemap->x0,&phasemap->y0,&phasemap->x1,&phasemap->y1) < 4) {printf("fscanf() failed. Exiting ungracefully.\n");exit(1);}
	//
	phasemap->dx = (phasemap->x1 - phasemap->x0)/(double)(phasemap->mx);
	phasemap->dy = (phasemap->y1 - phasemap->y0)/(double)(phasemap->my);


	/* allocate data */
	phasemap->data = malloc( sizeof(int)* phasemap->mx * phasemap->my );

	/* parse phase map from file */
  index = 0;
  for (j=0; j<phasemap->my; j++) {
    for (i=0; i<phasemap->mx; i++) {
      if (fscanf(fp,"%d",&phasemap->data[index]) < 1) {printf("fscanf() failed. Exiting ungracefully.\n");exit(1);}
      //printf("%d \n", phasemap->data[index]);
      index++;
    }
  }

	/* compute max number of phases */
	phasemap_max = -1;
	for (j=0; j<phasemap->my; j++) {
		for (i=0; i<phasemap->mx; i++) {
			int index;

			PhaseMapGetIndex(phasemap,i,j,&index);

			if (phasemap->data[index] > phasemap_max) {
				phasemap_max = phasemap->data[index];
			}
		}
	}
	if (phasemap_max==-1) {
		printf("Error(%s): Zero phases have been detected \n", __func__);
	}
	phasemap->nphases = phasemap_max+1;

	/* set pointer */
	*map = phasemap;
	fclose(fp);
}

int PhaseMapLoadFromFile_ASCII_ZIPPED(const char filename[],PhaseMap *map)
{
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Ascii zipped loading is not supported");
}


void PhaseMapLoadFromFile(const char filename[],PhaseMap *map)
{
	size_t len;
	int is_zipped;
	int matched_extension;

	is_zipped = 0;

	/* check extensions for common zipped file extensions */
	len = strlen(filename);
	matched_extension = strcmp(&filename[len-8],".tar.gz");
	if (matched_extension == 0) {
		printf("  Detected .tar.gz\n");
		is_zipped = 1;
	}
	matched_extension = strcmp(&filename[len-5],".tgz");
	if (matched_extension == 0) {
		printf("  Detected .tgz\n");
		is_zipped = 1;
	}
	matched_extension = strcmp(&filename[len-3],".Z");
	if (matched_extension == 0) {
		printf("  Detected .Z\n");
		is_zipped = 1;
	}

	if (is_zipped == 1) {
    PhaseMapLoadFromFile_ASCII_ZIPPED(filename,map);
	} else {
		PhaseMapLoadFromFile_ASCII(filename,map);
	}
}

void PhaseMapGetPhaseIndex(PhaseMap phasemap,double xp[],int *phase)
{
	int i,j,index;

	(*phase) = (int)PHASE_MAP_POINT_OUTSIDE;

	if (xp[0] < phasemap->x0) { return; }
	if (xp[0] > phasemap->x1) { return; }
	if (xp[1] < phasemap->y0) { return; }
	if (xp[1] > phasemap->y1) { return; }

	i = (xp[0] - phasemap->x0)/phasemap->dx;
	j = (xp[1] - phasemap->y0)/phasemap->dy;
	if (i==phasemap->mx) { i--; }
	if (j==phasemap->my) { j--; }

	PhaseMapGetIndex(phasemap,i,j,&index);

	*phase = phasemap->data[index];

}

void PhaseMapCheckValidity(PhaseMap phasemap,int phase,int *is_valid)
{
	*is_valid = 0;

	if ( (phase>=0) && (phase<phasemap->nphases) ) {
		*is_valid = 1;
	} else if ( phase==(int)PHASE_MAP_POINT_OUTSIDE ) {
		*is_valid = 0;
	}
}

/*

 gnuplot> set pm3d map
 gnuplot> splot "filename"

 */
void PhaseMapViewGnuplot(const char filename[],PhaseMap phasemap)
{
	FILE *fp = NULL;
	int i,j;

	/* open file to parse */
	fp = fopen(filename,"w");
	if (fp==NULL) {
		printf("Error(%s): Could not open file: %s \n",__func__, filename );
		exit(EXIT_FAILURE);
	}
	fprintf(fp,"# Phase map information \n");
	fprintf(fp,"# Phase map : (x0,y0) = (%lf,%lf) \n",phasemap->x0,phasemap->y0);
	fprintf(fp,"# Phase map : (x1,y1) = (%lf,%lf) \n",phasemap->x1,phasemap->y1);
	fprintf(fp,"# Phase map : (dx,dy) = (%lf,%lf) \n",phasemap->dx,phasemap->dy);
	fprintf(fp,"# Phase map : (mx,my) = (%d,%d) \n",phasemap->mx,phasemap->my);
	fprintf(fp,"# Phase map : nphases = %d \n",phasemap->nphases);

	for (j=0; j<phasemap->my; j++) {
		for (i=0; i<phasemap->mx; i++) {
			double x,y;
			int index;

			x = phasemap->x0 + phasemap->dx * 0.5 + i * phasemap->dx;
			y = phasemap->y0 + phasemap->dy * 0.5 + j * phasemap->dy;
			PhaseMapGetIndex(phasemap,i,j,&index);

			fprintf(fp,"%lf %lf %lf \n", x,y,(double)phasemap->data[index]);
		}fprintf(fp,"\n");
	}
	fclose(fp);

}

void PhaseMapGetMaxPhases(PhaseMap phasemap,int *maxphase)
{
	*maxphase = phasemap->nphases;
}

PetscErrorCode pTatinCtxAttachPhaseMap(pTatinCtx ctx,PhaseMap map)
{
	PetscErrorCode ierr;

    PetscFunctionBegin;

	ierr = pTatinCtxAttachModelData(ctx,"phasemap",(void*)map);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

PetscErrorCode pTatinCtxGetPhaseMap(pTatinCtx ctx,PhaseMap *map)
{
	void *mymap;
	PetscErrorCode ierr;

    PetscFunctionBegin;
	ierr = pTatinCtxGetModelData(ctx,"phasemap",&mymap);CHKERRQ(ierr);
	*map = (PhaseMap)mymap;

	PetscFunctionReturn(0);
}

/*
void test_PhaseMap(void)
{
	PhaseMap phasemap;
	double xp[2];
	int phase;

  //	PhaseMapLoadFromFile("test.bmp",&phasemap);
	PhaseMapLoadFromFile("model_geometry",&phasemap);

	xp[0] = 0.0;	xp[1] = 1.0;
	PhaseMapGetPhaseIndex(phasemap,xp,&phase);
	printf("x = ( %lf , %lf ) ==> phase = %d \n", xp[0],xp[1],phase);

	xp[0] = 5.0;	xp[1] = 3.2;
	PhaseMapGetPhaseIndex(phasemap,xp,&phase);
	printf("x = ( %lf , %lf ) ==> phase = %d \n", xp[0],xp[1],phase);

	xp[0] = -1.0;	xp[1] = 1.0;
	PhaseMapGetPhaseIndex(phasemap,xp,&phase);
	printf("x = ( %lf , %lf ) ==> phase = %d \n", xp[0],xp[1],phase);

	PhaseMapViewGnuplot("test.gp",phasemap);

	PhaseMapDestroy(&phasemap);

}
*/
