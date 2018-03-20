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
 **    filename:   spmA.c
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
	This is a do nothing surface process model (SPM) written in C.
  The purpose here is simply to demonstrate how to link against external codes/libraries.


*/
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "time.h"

#include "spmA.h"


int spmA_NID(int i,int j,int nx)
{
	return (i + j * nx );
}

void* spmA_malloc(size_t a)
{
	void *t;


	t = malloc(a);
	memset(t,0,a);
	return t;
}

int spmA_New(SPMAData **c)
{
	SPMAData *ctx;


	ctx = malloc( sizeof(SPMAData) );

	(*c) = ctx;

	return 1;
}

int spmA_Initialise(SPMAData *c,int nx,int ny,double Lx,double Ly,double k0,double k1,double e0,double dt,double final_time)
{
	int i,j,n;
	double dx,dy;

	c->output_frequency = 1;

	c->nx = nx;
	c->ny = ny;
	c->Lx  = Lx;
	c->Ly  = Ly;

	c->diffusion0  = k0;
	c->diffusion1  = k1;
	c->erodibility = e0;

	c->dt = dt;
	c->final_time = final_time;

	n = nx * ny;
	c->x = spmA_malloc(sizeof(double)*n);
	c->y = spmA_malloc(sizeof(double)*n);
	c->h = spmA_malloc(sizeof(double)*n);
	c->h_old = spmA_malloc(sizeof(double)*n);
	c->w = spmA_malloc(sizeof(double)*n);

	dx = Lx/(double)(nx-1);
	dy = Ly/(double)(ny-1);
	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) {
			int nid;

			nid = i + j * nx;
			c->x[nid] = 0.0 + i*dx;
			c->y[nid] = 0.0 + j*dy;
		}
	}

	return 1;
}


int spmA_InitialiseTopo_SlopeY(SPMAData *c,double h0,double h1)
{
	int i,j;
	double dy,dh;


	dy = c->Ly;
	dh = (h1 - h0)/dy;

	for (i=0; i<c->nx; i++) {
		for (j=0; j<c->ny; j++) {
			int nid;
			double ypos;

			nid = i + j * c->nx;
			ypos = c->y[nid];
			c->h_old[nid] = h0 + dh * ypos;
		}
	}

	return 1;
}


int spmA_InitialiseTopo_ApplyPositiveRandomNoise(SPMAData *c,double a0)
{
	int i,j;
	double r;

	srand(0);
	for (i=0; i<c->nx; i++) {
		for (j=0; j<c->ny; j++) {
			int nid;

			nid = i + j * c->nx;

			r = rand()/( (double)(RAND_MAX) );
			c->h_old[nid] += a0 * r;
		}
	}

	return 1;
}


int spmA_InitialiseUplift_Constant(SPMAData *c,double w0)
{
	int i,j;


	for (i=0; i<c->nx; i++) {
		for (j=0; j<c->ny; j++) {
			int nid;

			nid = i + j * c->nx;
			c->w[nid] = w0;
		}
	}

	return 1;
}


int spmA_InitialiseUplift_StepX(SPMAData *c,double x0,double w0,double w1)
{
	int i,j;


	for (i=0; i<c->nx; i++) {
		for (j=0; j<c->ny; j++) {
			int nid;
			double xpos;

			nid = i + j * c->nx;
			xpos = c->x[nid];

			c->w[nid] = w0;
			if (xpos >= x0) {
				c->w[nid] = w1;
			}
		}
	}

	return 1;
}


int spmA_UpdateTopo_basic(SPMAData *c,double dt)
{
	int i,j,nx;
	double dx,dy;
	double qR,qL,qU,qD;

	nx = c->nx;
	dx = c->Lx/(double)(c->nx-1);
	dy = c->Ly/(double)(c->ny-1);

	/* only uplift on interior plus edges */
	for (i=0; i<c->nx; i++) {
		for (j=0; j<c->ny; j++) {
			int nid;

			nid = i + j * nx;
			c->h[nid] += c->h_old[nid] + dt * c->w[nid];
		}
	}

	/* h_old = h */
	for (i=0; i<c->nx; i++) {
		for (j=0; j<c->ny; j++) {
			int nid;

			nid = i + j * nx;
			c->h_old[nid] = c->h[nid];
		}
	}

	return 1;
}

int spmA_ComputeTimestep_linear_diffusion(SPMAData *c,double *dt)
{
	double dx,dy,ds;

	dx = c->Lx/(double)(c->nx-1);
	dy = c->Ly/(double)(c->ny-1);

	ds = dx;
	if (dy < ds) { ds = dy; }

	*dt = 0.25 * ds*ds / c->diffusion0;

	return 1;
}

int spmA_UpdateTopo_linear_diffusion(SPMAData *c,double dt)
{
	int i,j,nx;
	double dx,dy;
	double qR,qL,qU,qD;

	nx = c->nx;
	dx = c->Lx/(double)(c->nx-1);
	dy = c->Ly/(double)(c->ny-1);

	/* diffusion on interior plus uplift */
	for (i=1; i<c->nx-1; i++) {
		for (j=1; j<c->ny-1; j++) {
			int nid;

			nid = i + j * nx;

			qR = c->diffusion0 * (c->h_old[spmA_NID(i+1,j,nx)] - c->h_old[spmA_NID(i,j,nx)])  /dx;
			qL = c->diffusion0 * (c->h_old[spmA_NID(i,j,nx)]   - c->h_old[spmA_NID(i-1,j,nx)])/dx;

			qU = c->diffusion0 * (c->h_old[spmA_NID(i,j+1,nx)] - c->h_old[spmA_NID(i,j,nx)])  /dy;
			qD = c->diffusion0 * (c->h_old[spmA_NID(i,j,nx)]   - c->h_old[spmA_NID(i,j-1,nx)])/dy;

			c->h[nid] = c->h_old[nid] + dt * ( (qR-qL)/dx + (qU-qD)/dy ) + dt * c->w[nid];
		}
	}


	/* fix h on edges */
	for (i=0; i<c->nx; i++) {
		for (j=0; j<c->ny; j++) {
			int nid;

			if ( (i>0) && (i<c->nx-1) ) {
				if ( (j>0) && (j<c->ny-1) ) {
					continue;
				}
			}

			nid = i + j * nx;
			c->h[nid] = c->h_old[nid] + dt * c->w[nid];
		}
	}


	/* h_old = h */
	for (i=0; i<c->nx; i++) {
		for (j=0; j<c->ny; j++) {
			int nid;

			nid = i + j * nx;
			c->h_old[nid] = c->h[nid];
		}
	}

	return 1;
}

int spmA_Apply(SPMAData *c)
{
	int    ierr;
	double dt,dt_mine,time;
	int    step,nt;
	char   name[100];
	SPMAEvolutionType evolution_type;

	nt = c->final_time/c->dt;
	nt++;
	dt_mine = c->final_time/((double)nt);

	printf("  SPMA: User requested dt=%1.4e, using dt=%1.4e \n",c->dt,dt_mine);
	evolution_type = SPMA_EVO_LINEAR_DIFFUSION;

	ierr = spmA_OutputIC(c,"spma_step_ic.out");

	time = 0.0;
	step = 0;
	while( time < c->final_time) {
		dt = c->dt;
		dt = dt_mine;

		/* evolution models */

		switch (evolution_type) {

			case SPMA_EVO_BASIC:
				ierr = spmA_UpdateTopo_basic(c,c->dt);
				break;

			case SPMA_EVO_LINEAR_DIFFUSION:

				ierr = spmA_ComputeTimestep_linear_diffusion(c,&dt);
				c->dt = dt;
				ierr = spmA_UpdateTopo_linear_diffusion(c,c->dt);
				break;

			case SPMA_EVO_NONLINEAR_DIFFUSION:
				printf("ERROR: SPMA_EVO_NONLINEAR_DIFFUSION unsupported!!\n");
				exit(0);
				break;

		}
		step++;
		time = time + c->dt;

		if (step%c->output_frequency == 0) {
			printf("  SPMA: performed step %d (time = %1.4e / %1.4e)\n",step,time,c->final_time);
			sprintf(name,"spma_step_%1.6d.out",step);
			ierr = spmA_Output(c,time,name);
		}
	}
	sprintf(name,"spma_step_%1.6d.out",step);
	ierr = spmA_Output(c,time,name);


	return 1;
}

int spmA_OutputIC(SPMAData *c,const char filename[])
{
	FILE *fp;
	int i,j;


	printf("  SPMA output: %s \n",filename);
	fp = fopen(filename,"w");
	if (fp == NULL) {
		printf("  SPMA: Cannot open file %s \n",filename);
		exit(0);
	}

	fprintf(fp,"# SPMA output \n");
	fprintf(fp,"# SPMA nx=%d ny=%d \n",c->nx,c->ny);
	fprintf(fp,"# SPMA Lx=%1.4e Ly=%1.4e \n",c->Lx,c->Ly);
	fprintf(fp,"# SPMA diffusion0=%1.4e diffusion1=%1.4e erodibility=%1.4e\n",c->diffusion0,c->diffusion1,c->erodibility);


	/* uplift on interior plus edges */
	for (i=0; i<c->nx; i++) {
		for (j=0; j<c->ny; j++) {
			int nid;

			nid = i + j * c->nx;
			fprintf(fp,"%1.4e %1.4e %1.4e %1.4e \n",c->x[nid],c->y[nid],c->h_old[nid],c->w[nid]);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

	return 1;
}

int spmA_Output(SPMAData *c,double time,const char filename[])
{
	FILE *fp;
	int i,j;


	printf("  SPMA output: %s \n",filename);
	fp = fopen(filename,"w");
	if (fp == NULL) {
		printf("  SPMA: Cannot open file %s \n",filename);
		exit(0);
	}

	fprintf(fp,"# SPMA output \n");
	fprintf(fp,"# SPMA nx=%d ny=%d \n",c->nx,c->ny);
	fprintf(fp,"# SPMA Lx=%1.4e Ly=%1.4e \n",c->Lx,c->Ly);
	fprintf(fp,"# SPMA diffusion0=%1.4e diffusion1=%1.4e erodibility=%1.4e\n",c->diffusion0,c->diffusion1,c->erodibility);
	fprintf(fp,"# SPMA time=%1.4e dt=%1.4e \n",time,c->dt);


	/* uplift on interior plus edges */
	for (i=0; i<c->nx; i++) {
		for (j=0; j<c->ny; j++) {
			int nid;

			nid = i + j * c->nx;
			fprintf(fp,"%1.4e %1.4e %1.4e \n",c->x[nid],c->y[nid],c->h[nid]);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

	return 1;
}

int spmA_Destroy(SPMAData **c)
{
	SPMAData *ctx;


	if (c) return(1);
	if (*c) return(1);


	ctx = (*c);
	free(ctx->x);
	free(ctx->y);
	free(ctx->h);
	free(ctx->h_old);
	free(ctx->w);

	free(ctx);
	(*c) = NULL;

	return 1;
}



