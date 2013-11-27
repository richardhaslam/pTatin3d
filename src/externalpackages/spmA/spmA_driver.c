
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"

#include "spmA.h"





int main(void)
{
	int nx,ny;
	double Lx,Ly;
	SPMAData *c;
	int ierr;
	
	ierr = spmA_New(&c);
	
	nx = 61;
	ny = 61;
	Lx = 10.0;
	Ly = 10.0;
	ierr = spmA_Initialise(c,nx,ny,Lx,Ly,1.0,0.0,0.0,0.01,10.0);

	ierr = spmA_InitialiseTopo_SlopeY(c,1.0,2.0);
	ierr = spmA_InitialiseTopo_ApplyPositiveRandomNoise(c,0.3);

	//ierr = spmA_InitialiseUplift_Constant(c,0.0);
	ierr = spmA_InitialiseUplift_StepX(c,3.0,2.2,3.3);
	
	c->output_frequency = 100;
	ierr = spmA_Apply(c);
	
	ierr = spmA_Destroy(&c);

	return 0;
}