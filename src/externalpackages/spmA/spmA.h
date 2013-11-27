
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
