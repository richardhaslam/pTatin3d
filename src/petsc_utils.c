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
    
    PetscFunctionBegin;
    ierr = PetscObjectGetComm((PetscObject)x,&comm);CHKERRQ(ierr);
    ierr = VecGetLocalSize(x,&n);CHKERRQ(ierr);
    ierr = VecGetArray(x,&array);CHKERRQ(ierr);
    
    switch (type) {
        
        case NORM_1:
        {
            PetscScalar sum = 0.0;
            
            for (i=0; i<n; i++) {
                sum += PetscAbsScalar(array[i]);
            }
            ierr = MPI_Reduce(&sum,val,1,MPIU_REAL,MPI_SUM,0,comm);CHKERRQ(ierr);
            
        }
            break;
        
        case NORM_2:
        {
            PetscScalar rsum,sum = 0.0;

            for (i=0; i<n; i++) {
                sum += array[i]*(PetscConj(array[i]));
            }
            rsum = PetscRealPart(sum);
            ierr = MPI_Reduce(&rsum,val,1,MPIU_REAL,MPI_SUM,0,comm);CHKERRQ(ierr);
            *val = PetscSqrtReal(*val);
        }
            break;
        
        case NORM_1_AND_2:
        {
            PetscScalar sum = 0.0;
            PetscReal   sum2[] = { 0.0, 0.0 };
            
            for (i=0; i<n; i++) {
                sum2[0] += PetscAbsScalar(array[i]);
                sum     += array[i]*(PetscConj(array[i]));
            }
            sum2[1] = PetscRealPart(sum);
            ierr = MPI_Reduce(sum2,val,2,MPIU_REAL,MPI_SUM,0,comm);CHKERRQ(ierr);
            val[1] = PetscSqrtReal(val[1]);
        }
            break;
            
        case NORM_INFINITY:
        {
            PetscReal max,abs;
            
            max = 0.0;
            for (i=0; i<n; i++) {
                abs = PetscAbsScalar(array[i]);
                if (abs > max) { max = abs; }
                /* check special case of abs == NaN */
                if (abs != abs) {
                    max = abs;
                    break;
                }
            }
            ierr = MPI_Reduce(&max,val,1,MPIU_REAL,MPI_MAX,0,comm);CHKERRQ(ierr);
        }
            break;

        default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown norm type");
            break;
    }
    ierr = VecRestoreArray(x,&array);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecStrideNormLocal"
/*@
 VecStrideNormLocal - Computes the norm of subvector of a vector defined
 by a starting point and a stride on the first rank in communicator
 
 Collective on Vec
 
 Input Parameter:
 +  v - the vector
 .  start - starting point of the subvector (defined by a stride)
 -  ntype - type of norm, one of NORM_1, NORM_2, NORM_INFINITY
 
 Output Parameter:
 .  norm - the norm
 
 Notes:
 One must call VecSetBlockSize() before this routine to set the stride
 information, or use a vector created from a multicomponent DMDA.
 
 If x is the array representing the vector x then this computes the norm
 of the array (x[start],x[start+stride],x[start+2*stride], ....)
 
 This is useful for computing, say the norm of the pressure variable when
 the pressure is stored (interlaced) with other variables, say density etc.
 
 This will only work if the desire subvector is a stride subvector
 
 Level: advanced
 
 Concepts: norm^on stride of vector
 Concepts: stride^norm
 
 .seealso: VecNorm(), VecStrideGather(), VecStrideScatter(), VecStrideMin(), VecStrideMax()
 @*/
PetscErrorCode VecStrideNormLocal(Vec v,PetscInt start,NormType ntype,PetscReal *nrm)
{
    PetscErrorCode ierr;
    PetscInt       i,n,bs;
    PetscScalar    *x;
    PetscReal      tnorm;
    MPI_Comm       comm;
    
    PetscFunctionBegin;
    ierr = VecGetLocalSize(v,&n);CHKERRQ(ierr);
    ierr = VecGetArray(v,&x);CHKERRQ(ierr);
    ierr = PetscObjectGetComm((PetscObject)v,&comm);CHKERRQ(ierr);
    
    ierr = VecGetBlockSize(v,&bs);CHKERRQ(ierr);
    if (start < 0) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Negative start %D",start);
    else if (start >= bs) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Start of stride subvector (%D) is too large for stride\n Have you set the vector blocksize (%D) correctly with VecSetBlockSize()?",start,bs);
    x += start;
    
    switch (ntype) {
        case NORM_2:
        {
            PetscScalar sum = 0.0;
            for (i=0; i<n; i+=bs) {
                sum += x[i]*(PetscConj(x[i]));
            }
            tnorm = PetscRealPart(sum);
            ierr  = MPI_Reduce(&tnorm,nrm,1,MPIU_REAL,MPIU_SUM,0,comm);CHKERRQ(ierr);
            *nrm  = PetscSqrtReal(*nrm);
            break;
        }
        
        case NORM_1:
        {
            tnorm = 0.0;
            for (i=0; i<n; i+=bs) {
                tnorm += PetscAbsScalar(x[i]);
            }
            ierr = MPI_Reduce(&tnorm,nrm,1,MPIU_REAL,MPIU_SUM,0,comm);CHKERRQ(ierr);
            break;
        }

        case NORM_INFINITY:
        {
            PetscReal tmp;
            tnorm = 0.0;
            
            for (i=0; i<n; i+=bs) {
                if ((tmp = PetscAbsScalar(x[i])) > tnorm) tnorm = tmp;
                /* check special case of tmp == NaN */
                if (tmp != tmp) {
                    tnorm = tmp;
                    break;
                }
            }
            ierr = MPI_Reduce(&tnorm,nrm,1,MPIU_REAL,MPIU_MAX,0,comm);CHKERRQ(ierr);
            break;
        }
            
        default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_UNKNOWN_TYPE,"Unknown norm type");
            break;
    }
    ierr = VecRestoreArray(v,&x);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
