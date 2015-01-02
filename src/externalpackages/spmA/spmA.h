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
 **    filename:   spmA.h
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

#ifndef __spmA_dummy_h__
#define __spmA_dummy_h__

typedef struct {
	int    nx,ny;
	double Lx,Ly;
	double base_level;
	double *x,*y,*h,*h_old,*w;
	double diffusion0,diffusion1,erodibility;
	double dt,final_time;
	int output_frequency;
} SPMAData;

typedef enum { SPMA_EVO_BASIC=1, SPMA_EVO_LINEAR_DIFFUSION, SPMA_EVO_NONLINEAR_DIFFUSION } SPMAEvolutionType;

int spmA_NID(int i,int j,int nx);
void* spmA_malloc(size_t a);

int spmA_InitialiseTopo_SlopeY(SPMAData *c,double h0,double h1);
int spmA_InitialiseUplift_Constant(SPMAData *c,double w0);
int spmA_InitialiseUplift_StepX(SPMAData *c,double x0,double w0,double w1);
int spmA_InitialiseTopo_ApplyPositiveRandomNoise(SPMAData *c,double a0);

int spmA_UpdateTopo_basic(SPMAData *c,double dt);
int spmA_ComputeTimestep_linear_diffusion(SPMAData *c,double *dt);
int spmA_UpdateTopo_linear_diffusion(SPMAData *c,double dt);


int spmA_New(SPMAData **c);
int spmA_Initialise(SPMAData *c,int nx,int ny,double Lx,double Ly,double k0,double k1,double e0,double dt,double final_time);
int spmA_Apply(SPMAData *c);
int spmA_OutputIC(SPMAData *c,const char filename[]);
int spmA_Output(SPMAData *c,double time,const char filename[]);
int spmA_Destroy(SPMAData **c);

#endif
