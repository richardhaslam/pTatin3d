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
 **    filename:   test_mpiiob.c
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

#include <petsc.h>
#include <ptatin3d.h>
#include <ptatin_init.h>
#include <mpiio_blocking.h>


PetscErrorCode test_single_write(PetscInt length_ref,PetscBool skip_header)
{
  PetscErrorCode ierr;
  PetscMPIInt rank,root = 0;
  FILE *fp = NULL;
  double *data;
  int length,i;
  long int bytes,bytes_sum_g;
  PetscLogDouble t0,t1;
  
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);CHKERRQ(ierr);
  if (rank == root) {
    fp = fopen("test-ex1.bin","w");
  }
  
  length = length_ref + rank*3;
  data = malloc(length*sizeof(double));
  for (i=0; i<length; i++) {
    data[i] = 100.0*rank + i;
  }
  
  PetscTime(&t0);
  ierr = MPIWrite_Blocking(fp,data,length,sizeof(double),root,skip_header,PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscTime(&t1);
  
  bytes = length * sizeof(double);
  ierr = MPI_Allreduce(&bytes,&bytes_sum_g,1, MPI_LONG, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,".. [single][MPIWrite_Blocking] Wrote %1.2e (MB) in %1.2e (sec)\n",(double)bytes_sum_g/1.0e6,t1-t0);
  
  free(data);
  data = NULL;
  
  if (rank == root) {
    fclose(fp);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode test_single_read(PetscInt length_ref,PetscBool skip_header)
{
  PetscErrorCode ierr;
  PetscMPIInt commsize,rank,r,root = 0;
  FILE *fp = NULL;
  double *data = NULL;
  int length,i,start;
  long int bytes,bytes_sum_g;
  double diff,diff_g;
  PetscLogDouble t0,t1;
  
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);CHKERRQ(ierr);
  if (rank == root) {
    fp = fopen("test-ex1.bin","r");
  }
  
  length = length_ref + rank*3;
  start = 0;
  for (r=0; r<rank; r++) {
    start += (length_ref + r*3);
  }

  PetscTime(&t0);
  ierr = MPIRead_Blocking(fp,(void**)&data,length,sizeof(double),root,skip_header,PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscTime(&t1);
  
  bytes = length * sizeof(double);
  ierr = MPI_Allreduce(&bytes,&bytes_sum_g,1, MPI_LONG, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,".. [single][MPIRead_Blocking] Read %1.2e (MB) in %1.2e (sec)\n",(double)bytes_sum_g/1.0e6,t1-t0);
  
  /*
  if (rank == (commsize-1)) {
    for (i=0; i<length; i++) {
      printf("[rank %d][%d] %1.4e\n",(int)rank,start+i,data[i]);
    }
  }
  */
   
  diff = 0.0;
  for (i=0; i<length; i++) {
    diff += fabs( data[i] - (double)(100.0*rank + i) );
  }
  ierr = MPI_Allreduce(&diff,&diff_g,1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
  
  if (diff_g < 1.0e-12) {
    PetscPrintf(PETSC_COMM_WORLD,"++ [single][Pass] Max difference between expected and data read is < 1.0e-12\n");
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"-- [single][Fail] Max difference between expected and data read is %+1.4e\n",diff_g);
  }
  
  free(data);
  data = NULL;
  
  if (rank == root) {
    fclose(fp);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode test_multidata_write(PetscInt length_ref,PetscBool skip_header)
{
  PetscErrorCode ierr;
  PetscMPIInt rank,root = 0;
  FILE *fp = NULL;
  double *data1;
  int   *data2;
  short *data3;
  long int length[3],i;
  long int bytes,bytes_sum_g;
  PetscLogDouble t0,t1;
  
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);CHKERRQ(ierr);
  
  length[0] = length_ref + (rank+1)*3;
  length[1] = length_ref + (rank+2)*2;
  length[2] = length_ref + (rank+3)*1;
  
  data1 = malloc(length[0]*sizeof(double));
  data2 = malloc(length[1]*sizeof(int));
  data3 = malloc(length[2]*sizeof(short));
  
  for (i=0; i<length[0]; i++) {
    data1[i] = (double)(100.0*rank + i);
  }
  for (i=0; i<length[1]; i++) {
    data2[i] = (int)(10 * rank + i);
  }
  for (i=0; i<length[2]; i++) {
    data3[i] = (short)(2*rank + i);
  }
  
  if (rank == root) {
    fp = fopen("test_multi.bin","w");
  }

  PetscTime(&t0);
  ierr = MPIWrite_Blocking(fp,data1,length[0],sizeof(double),root,skip_header,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPIWrite_Blocking(fp,data2,length[1],sizeof(int),root,skip_header,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPIWrite_Blocking(fp,data3,length[2],sizeof(short),root,skip_header,PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscTime(&t1);
  
  bytes  = length[0] * sizeof(double);
  bytes += length[1] * sizeof(int);
  bytes += length[2] * sizeof(short);
  ierr = MPI_Allreduce(&bytes,&bytes_sum_g,1, MPI_LONG, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,".. [multi][MPIWrite_Blocking x 3] Wrote %1.2e (MB) in %1.2e (sec)\n",(double)bytes_sum_g/1.0e6,t1-t0);

  free(data3);
  free(data2);
  free(data1);
  
  if (rank == root) {
    fclose(fp);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode test_multidata_read(PetscInt length_ref,PetscBool skip_header)
{
  PetscErrorCode ierr;
  PetscMPIInt rank,root = 0;
  FILE *fp = NULL;
  double *data1 = NULL;
  int    *data2 = NULL;
  short  *data3 = NULL;
  long int i,length[3];
  double diff[3],diff_g[3],diff1,diff2,diff3;
  long int bytes,bytes_sum_g;
  PetscLogDouble t0,t1;
  
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);CHKERRQ(ierr);
  
  length[0] = length_ref + (rank+1)*3;
  length[1] = length_ref + (rank+2)*2;
  length[2] = length_ref + (rank+3)*1;
  
  if (rank == root) {
    fp = fopen("test_multi.bin","r");
  }
  
  PetscTime(&t0);
  ierr = MPIRead_Blocking(fp,(void**)&data1,length[0],sizeof(double),root,skip_header,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPIRead_Blocking(fp,(void**)&data2,length[1],sizeof(int),root,skip_header,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPIRead_Blocking(fp,(void**)&data3,length[2],sizeof(short),root,skip_header,PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscTime(&t1);
  
  bytes  = length[0] * sizeof(double);
  bytes += length[1] * sizeof(int);
  bytes += length[2] * sizeof(short);
  ierr = MPI_Allreduce(&bytes,&bytes_sum_g,1, MPI_LONG, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,".. [multi][MPIRead_Blocking x 3] Read %1.2e (MB) in %1.2e (sec)\n",(double)bytes_sum_g/1.0e6,t1-t0);
  
  
  diff1 = 0.0;
  for (i=0; i<length[0]; i++) {
    diff1 += fabs( data1[i] - (double)(100.0*rank + i) );
  }
  
  diff2 = 0.0;
  for (i=0; i<length[1]; i++) {
    diff2 += fabs( (double)(data2[i] - (int)(10 * rank + i)) );
  }
  
  diff3 = 0.0;
  for (i=0; i<length[2]; i++) {
    diff3 += fabs( (double)(data3[i] - (short)(2*rank + i)) );
  }

  diff[0] = diff1;
  diff[1] = diff2;
  diff[2] = diff3;
  ierr = MPI_Allreduce(diff,diff_g,3, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
  
  if (diff_g[0] < 1.0e-12) {
    PetscPrintf(PETSC_COMM_WORLD,"++ [multi][Pass] Max difference between expected and data read is < 1.0e-12 <field1>\n");
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"-- [multi][Fail] Max difference between expected and data read is %+1.4e <field1>\n",diff_g[0]);
  }
  if (diff_g[1] < 1.0e-12) {
    PetscPrintf(PETSC_COMM_WORLD,"++ [multi][Pass] Max difference between expected and data read is < 1.0e-12 <field2>\n");
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"-- [multi][Fail] Max difference between expected and data read is %+1.4e <field2>\n",diff_g[1]);
  }
  if (diff_g[2] < 1.0e-12) {
    PetscPrintf(PETSC_COMM_WORLD,"++ [multi][Pass] Max difference between expected and data read is < 1.0e-12 <field3>\n");
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"-- [multi][Fail] Max difference between expected and data read is %+1.4e <field3>\n",diff_g[2]);
  }
  
  
  free(data3);
  free(data2);
  free(data1);
  
  if (rank == root) {
    fclose(fp);
  }
  PetscFunctionReturn(0);
}


int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscInt       test_id,length;
  
  ierr = pTatinInitialize(&argc,&argv,(char *)0,NULL);CHKERRQ(ierr);
  
  test_id = 1;
  ierr = PetscOptionsGetInt(NULL,NULL,"-test_id",&test_id,NULL);CHKERRQ(ierr);
  length = 5;
  ierr = PetscOptionsGetInt(NULL,NULL,"-data_len_fac",&length,NULL);CHKERRQ(ierr);
  
  switch (test_id) {
      
    case 1:
      ierr = test_single_write(length,PETSC_TRUE);CHKERRQ(ierr);
      break;
      
    case 2:
      ierr = test_single_read(length,PETSC_TRUE);CHKERRQ(ierr);
      break;

    case 3:
      ierr = test_multidata_write(length,PETSC_FALSE);CHKERRQ(ierr);
      break;
    case 4:
      ierr = test_multidata_read(length,PETSC_FALSE);CHKERRQ(ierr);
      break;

    default:
      break;
  }
  
  ierr = pTatinFinalize();CHKERRQ(ierr);
  return 0;
}
