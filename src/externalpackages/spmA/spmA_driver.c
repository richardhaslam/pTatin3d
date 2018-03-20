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
 **    filename:   spmA_driver.c
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
