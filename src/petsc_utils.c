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
 **    filename:   petsc_utils.c
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

#include "petsc.h"
#include "petscvec.h"

/*
 Performs reduction on first rank in communicator
 Useful for KSP/SNES monitors which do not require result on all ranks (global reduction)
*/
#undef __FUNCT__
#define __FUNCT__ "VecNormLocal"
PetscErrorCode VecNormLocal(Vec x,NormType type,PetscReal *val)
{
    PetscScalar    *array;
    PetscInt       i,n;
    MPI_Comm       comm;
    PetscErrorCode ierr;
    
    ierr = PetscObjectGetComm((PetscObject)x,&comm);CHKERRQ(ierr);
    ierr = VecGetLocalSize(x,&n);CHKERRQ(ierr);
    ierr = VecGetArray(x,&array);CHKERRQ(ierr);
    
    switch (type) {
        
        case NORM_1:
        {
            PetscScalar sum;
            
            sum = 0.0;
            for (i=0; i<n; i++) {
                sum += PetscAbsScalar(array[i]);
            }
            ierr = MPI_Reduce(&sum,val,1,MPIU_REAL,MPI_SUM,0,comm);CHKERRQ(ierr);
            
        }
            break;
        
        case NORM_2:
        {
            PetscReal sum;

            sum = 0.0;
            for (i=0; i<n; i++) {
                sum = sum + array[i]*array[i];
            }
            ierr = MPI_Reduce(&sum,val,1,MPIU_REAL,MPI_SUM,0,comm);CHKERRQ(ierr);
            *val = PetscSqrtReal(*val);
        }
            break;
        
        case NORM_1_AND_2:
        {
            PetscReal sum2[] = { 0.0, 0.0 };
            
            for (i=0; i<n; i++) {
                sum2[0] += PetscAbsScalar(array[i]);
                sum2[1] += array[i]*array[i];
            }
            ierr = MPI_Reduce(sum2,val,2,MPIU_REAL,MPI_SUM,0,comm);CHKERRQ(ierr);
            val[1] = PetscSqrtReal(val[1]);
        }
            break;
            
        case NORM_INFINITY:
        {
            PetscReal max,abs;
            
            max = PETSC_MIN_REAL;
            for (i=0; i<n; i++) {
                abs = PetscAbsScalar(array[i]);
                max = PetscMax(max,abs);
            }
            ierr = MPI_Reduce(&max,val,1,MPIU_REAL,MPI_MAX,0,comm);CHKERRQ(ierr);
        }
            break;
    }
    ierr = VecRestoreArray(x,&array);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}
