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
 **    filename:   stokes_operators_tensor.c
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
// -*- indent-tabs-mode:t c-basic-offset:8 -*-

#include <petscsys.h>
#include <petscfe.h>
#include <ptatin3d.h>
#include <ptatin3d_stokes.h>
#include <dmda_element_q2p1.h>
#include <stokes_operators.h>

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

extern PetscLogEvent MAT_MultMFA11_stp;
extern PetscLogEvent MAT_MultMFA11_cto;
extern PetscLogEvent MAT_MultMFA11_ker;
extern PetscLogEvent MAT_MultMFA11_cfr;

typedef struct _p_MFA11OpenCL *MFA11OpenCL;

struct _p_MFA11OpenCL {
  PetscObjectState state;
  cl_context       context;
  cl_command_queue queue;
  cl_program       program;
  cl_kernel        kernel0;
  cl_kernel        kernel1;

  cl_mem ufield;
  cl_mem LA_gcoords;
  cl_mem gaussdata_w;  // Data at Gauss points multiplied by respective quadrature weight
  PetscInt   element_colors;
  PetscInt  *elements_per_color;
  cl_mem    *el_ids_colored;
  cl_mem elnidx_u;
  cl_mem Yu;
  cl_mem D;
  cl_mem B;
};

/* OpenCL error checking */
#define ERROR_CHECKER_CASE(ERRORCODE)  case ERRORCODE: SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_LIB,"OpenCL Error: %d",ERRORCODE)
static PetscErrorCode checkError(cl_int err)
{
  PetscFunctionBeginUser;
  if (err != CL_SUCCESS)
  {
    switch (err)
    {
      ERROR_CHECKER_CASE(CL_DEVICE_NOT_FOUND);
      ERROR_CHECKER_CASE(CL_DEVICE_NOT_AVAILABLE);
      ERROR_CHECKER_CASE(CL_COMPILER_NOT_AVAILABLE);
      ERROR_CHECKER_CASE(CL_MEM_OBJECT_ALLOCATION_FAILURE);
      ERROR_CHECKER_CASE(CL_OUT_OF_RESOURCES);
      ERROR_CHECKER_CASE(CL_OUT_OF_HOST_MEMORY);
      ERROR_CHECKER_CASE(CL_PROFILING_INFO_NOT_AVAILABLE);
      ERROR_CHECKER_CASE(CL_MEM_COPY_OVERLAP);
      ERROR_CHECKER_CASE(CL_IMAGE_FORMAT_MISMATCH);
      ERROR_CHECKER_CASE(CL_IMAGE_FORMAT_NOT_SUPPORTED);
      ERROR_CHECKER_CASE(CL_BUILD_PROGRAM_FAILURE);
      ERROR_CHECKER_CASE(CL_MAP_FAILURE);

      ERROR_CHECKER_CASE(CL_INVALID_VALUE);
      ERROR_CHECKER_CASE(CL_INVALID_DEVICE_TYPE);
      ERROR_CHECKER_CASE(CL_INVALID_PLATFORM);
      ERROR_CHECKER_CASE(CL_INVALID_DEVICE);
      ERROR_CHECKER_CASE(CL_INVALID_CONTEXT);
      ERROR_CHECKER_CASE(CL_INVALID_QUEUE_PROPERTIES);
      ERROR_CHECKER_CASE(CL_INVALID_COMMAND_QUEUE);
      ERROR_CHECKER_CASE(CL_INVALID_HOST_PTR);
      ERROR_CHECKER_CASE(CL_INVALID_MEM_OBJECT);
      ERROR_CHECKER_CASE(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR);
      ERROR_CHECKER_CASE(CL_INVALID_IMAGE_SIZE);
      ERROR_CHECKER_CASE(CL_INVALID_SAMPLER);
      ERROR_CHECKER_CASE(CL_INVALID_BINARY);
      ERROR_CHECKER_CASE(CL_INVALID_BUILD_OPTIONS);
      ERROR_CHECKER_CASE(CL_INVALID_PROGRAM);
      ERROR_CHECKER_CASE(CL_INVALID_PROGRAM_EXECUTABLE);
      ERROR_CHECKER_CASE(CL_INVALID_KERNEL_NAME);
      ERROR_CHECKER_CASE(CL_INVALID_KERNEL_DEFINITION);
      ERROR_CHECKER_CASE(CL_INVALID_KERNEL);
      ERROR_CHECKER_CASE(CL_INVALID_ARG_INDEX);
      ERROR_CHECKER_CASE(CL_INVALID_ARG_VALUE);
      ERROR_CHECKER_CASE(CL_INVALID_ARG_SIZE);
      ERROR_CHECKER_CASE(CL_INVALID_KERNEL_ARGS);
      ERROR_CHECKER_CASE(CL_INVALID_WORK_DIMENSION);
      ERROR_CHECKER_CASE(CL_INVALID_WORK_GROUP_SIZE);
      ERROR_CHECKER_CASE(CL_INVALID_WORK_ITEM_SIZE);
      ERROR_CHECKER_CASE(CL_INVALID_GLOBAL_OFFSET);
      ERROR_CHECKER_CASE(CL_INVALID_EVENT_WAIT_LIST);
      ERROR_CHECKER_CASE(CL_INVALID_EVENT);
      ERROR_CHECKER_CASE(CL_INVALID_OPERATION);
      ERROR_CHECKER_CASE(CL_INVALID_GL_OBJECT);
      ERROR_CHECKER_CASE(CL_INVALID_BUFFER_SIZE);
      ERROR_CHECKER_CASE(CL_INVALID_MIP_LEVEL);
      ERROR_CHECKER_CASE(CL_INVALID_GLOBAL_WORK_SIZE);
        
      default: SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_LIB,"Unknown OpenCL error. Maybe OpenCL SDK not properly installed?");
    }
  }
  PetscFunctionReturn(0);
}

#define ERR_CHECK(err) checkError(err);



#define NQP 27			/* Number of quadrature points per element; must equal Q2_NODES_PER_EL_3D (27) */

#define WARPS_PER_BLOCK    4


static const char * opencl_spmv_kernel_sources =
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable \n"
"#define NQP                 27 \n"
"#define Q2_NODES_PER_EL_3D  27 \n"
"#define PetscReal           double \n"
"#define PetscScalar         double \n"
"#define PetscInt            int \n"
"#define WARPS_PER_BLOCK     4 \n"


"void TensorContract(PetscReal const *R, PetscReal const *S,PetscReal const *T,PetscReal const x[],PetscReal y[]) \n"
"{ \n"
"  __local PetscReal u[WARPS_PER_BLOCK][NQP],v[WARPS_PER_BLOCK][NQP]; \n"
"  PetscInt warp_in_block = get_local_id(0) / 32; \n"
"  PetscInt id_in_warp = get_local_id(0) % 32; \n"

"  PetscInt c = id_in_warp % 3; \n"
"  PetscInt kj = id_in_warp / 3; \n"
"  PetscInt k3 = (id_in_warp / 9) * 3; \n"
"  PetscInt ji = id_in_warp % 9; \n"

"  for (PetscInt l=0; l<3; l++) { \n"

	// u[l,k,j,c] = R[c,i] x[l,k,j,i]
"    PetscReal result = 0; \n"
"    v[warp_in_block][id_in_warp] = x[l]; \n"
"    for (PetscInt i=0; i<3; i++) result += R[i] * v[warp_in_block][kj*3+i]; \n"
"    u[warp_in_block][id_in_warp] = result; \n"

	// v[l,k,b,c] = S[b,j] u[l,k,j,c]
"    result = 0; \n"
"    for (PetscInt j=0; j<3; j++) result += S[j] * u[warp_in_block][(k3+j)*3+c]; \n"
"    v[warp_in_block][id_in_warp] = result; \n"

	// y[l,a,b,c] = T[a,k] v[l,k,b,c]
"	for (PetscInt k=0; k<3; k++) y[l] += T[k] * v[warp_in_block][k*9+ji]; \n"

"  } // for l \n"
"} \n"

" void JacobianInvert(PetscScalar dx[3][3],PetscScalar *dxdet) \n"
"{ \n"
"		PetscScalar a[3][3]; \n"
"		PetscScalar b0,b3,b6,idet; \n"
"		for (PetscInt j=0; j<3; j++) { \n"
"			for (PetscInt k=0; k<3; k++) { \n"
"				a[j][k] = dx[j][k]; \n"
"			} \n"
"		} \n"
"		b0 =  (a[1][1]*a[2][2] - a[2][1]*a[1][2]); \n"
"		b3 = -(a[1][0]*a[2][2] - a[2][0]*a[1][2]); \n"
"		b6 =  (a[1][0]*a[2][1] - a[2][0]*a[1][1]); \n"
"		*dxdet = a[0][0]*b0 + a[0][1]*b3 + a[0][2]*b6; \n"
"		idet = 1.0 / *dxdet; \n"
"		dx[0][0] =  idet*b0; \n"
"		dx[0][1] = -idet*(a[0][1]*a[2][2] - a[2][1]*a[0][2]); \n"
"		dx[0][2] =  idet*(a[0][1]*a[1][2] - a[1][1]*a[0][2]); \n"
"		dx[1][0] =  idet*b3; \n"
"		dx[1][1] =  idet*(a[0][0]*a[2][2] - a[2][0]*a[0][2]); \n"
"		dx[1][2] = -idet*(a[0][0]*a[1][2] - a[1][0]*a[0][2]); \n"
"		dx[2][0] =  idet*b6; \n"
"		dx[2][1] = -idet*(a[0][0]*a[2][1] - a[2][0]*a[0][1]); \n"
"		dx[2][2] =  idet*(a[0][0]*a[1][1] - a[1][0]*a[0][1]); \n"
"} \n"

"void QuadratureAction(PetscScalar gaussdata_eta_w_dxdet, \n"
"				       PetscScalar const dx[3][3], \n"
"				       PetscScalar const du[3][3], \n"
"				       PetscScalar dv[3][3]) \n"
"{ \n"
"		PetscScalar dux[3][3]; \n"
"		for (PetscInt l=0; l<3; l++) { // fields \n"
"			for (PetscInt k=0; k<3; k++) { // directions \n"
"				dux[k][l] = du[0][l] * dx[k][0] + du[1][l] * dx[k][1] + du[2][l] * dx[k][2]; \n"
"			} \n"
"		} \n"

"		PetscScalar dvx[3][3]; \n"
"		dvx[0][0] = 2 * gaussdata_eta_w_dxdet * dux[0][0]; \n"
"		dvx[0][1] =     gaussdata_eta_w_dxdet * (dux[0][1] + dux[1][0]); \n"
"		dvx[0][2] =     gaussdata_eta_w_dxdet * (dux[0][2] + dux[2][0]); \n"
"		dvx[1][0] =     dvx[0][1]; \n"
"		dvx[1][1] = 2 * gaussdata_eta_w_dxdet * dux[1][1]; \n"
"		dvx[1][2] =     gaussdata_eta_w_dxdet * (dux[1][2] + dux[2][1]); \n"
"		dvx[2][0] =     dvx[0][2]; \n"
"		dvx[2][1] =     dvx[1][2]; \n"
"		dvx[2][2] = 2 * gaussdata_eta_w_dxdet * dux[2][2]; \n"

"		for (PetscInt l=0; l<3; l++) { // fields \n"
"			for (PetscInt k=0; k<3; k++) { // directions \n"
"				dv[k][l] = (dvx[0][l] * dx[0][k] + dvx[1][l] * dx[1][k] + dvx[2][l] * dx[2][k]); \n"
"			} \n"
"		} \n"
"} \n"

"__kernel void set_zero(__global PetscScalar *Yu, PetscInt localsize) \n"
"{ \n"
"    for (int i=get_global_id(0); i<localsize; i+=get_global_size(0)) Yu[i] = 0; \n"
"} \n"

"__kernel void MFStokesWrapper_A11(PetscInt nel,PetscInt nen_u,__global PetscInt const *el_ids_colored, __global PetscInt const *elnidx_u,__global PetscReal const *LA_gcoords, \n"
"                                  __global PetscScalar const *ufield,__global PetscReal const *gaussdata_w,__global PetscScalar *Yu, \n"
"                                  __constant PetscReal const *D,__constant PetscReal const *B) \n"
"{ \n"
"	PetscScalar el_x[3]; \n"
"	PetscScalar el_uv[3]; // unifies elu, elv \n"
"	PetscScalar dx[3][3]={ {0} },du[3][3]={ {0} },dv[3][3]={ {0} }; \n"
"   PetscScalar dxdet = 0; \n"
"   PetscInt elidx = get_global_id(0) / 32;  // one warp per colored element. elidx is here the index within the same color. \n"
"   PetscInt id_in_warp = get_local_id(0) % 32; \n"
"   PetscInt E_times_3; \n"
"   PetscReal R[3],S[3],T[3]; \n"
"   PetscInt c = id_in_warp % 3; \n"
"   PetscInt b = (id_in_warp % 9) / 3; \n"
"   PetscInt a = id_in_warp / 9; \n"

"   if (elidx >= nel) \n"
"     return; \n"

"	if (id_in_warp < Q2_NODES_PER_EL_3D) { \n"

"      elidx = el_ids_colored[elidx]; // get global element index \n"
"      E_times_3 = 3 * elnidx_u[nen_u*elidx+id_in_warp]; \n"

"      for (PetscInt l=0; l<3; l++) { \n"
"        el_x[l] = LA_gcoords[E_times_3+l]; \n"
"        el_uv[l] = ufield[E_times_3+l]; \n"
"        R[l] = D[3*c+l]; \n"
"        S[l] = B[3*b+l]; \n"
"        T[l] = B[3*a+l]; \n"
"      } \n"
"	  TensorContract(R,S,T,el_x, dx[0]); //TensorContract(D,B,B,GRAD,el_uxv,dx[0]); \n"
"	  TensorContract(R,S,T,el_uv,du[0]); //TensorContract(D,B,B,GRAD,el_uxv,du[0]); \n"

"      for (PetscInt l=0; l<3; l++) { \n"
"        R[l] = B[3*c+l]; \n"
"        S[l] = D[3*b+l]; \n"
"      } \n"
"	  TensorContract(R,S,T,el_x, dx[1]); //TensorContract(B,D,B,GRAD,el_uxv,dx[1]); \n"
"	  TensorContract(R,S,T,el_uv,du[1]); //TensorContract(B,D,B,GRAD,el_uxv,du[1]); \n"

"      for (PetscInt l=0; l<3; l++) { \n"
"        S[l] = B[3*b+l]; \n"
"        T[l] = D[3*a+l]; \n"
"      } \n"
"	  TensorContract(R,S,T,el_x, dx[2]); //TensorContract(B,B,D,GRAD,el_uxv,dx[2]); \n"
"	  TensorContract(R,S,T,el_uv,du[2]); //TensorContract(B,B,D,GRAD,el_uxv,du[2]); \n"

"	  JacobianInvert(dx,&dxdet); \n"

"	  QuadratureAction(gaussdata_w[elidx*NQP + id_in_warp] * dxdet,dx,du,dv); \n"

"      for (PetscInt l=0; l<3; l++) { \n"
"        el_uv[l] = 0; \n"
"        R[l] = D[3*l + c]; \n"
"        S[l] = B[3*l + b]; \n"
"        T[l] = B[3*l + a]; \n"
"      } \n"
"	  TensorContract(R,S,T,dv[0],el_uv); //TensorContract(D,B,B,GRAD_TRANSPOSE,dv[0],el_uxv); \n"
"      for (PetscInt l=0; l<3; l++) { \n"
"        R[l] = B[3*l + c]; \n"
"        S[l] = D[3*l + b]; \n"
"      } \n"
"	  TensorContract(R,S,T,dv[1],el_uv); //TensorContract(B,D,B,GRAD_TRANSPOSE,dv[1],el_uxv); \n"
"      for (PetscInt l=0; l<3; l++) { \n"
"        S[l] = B[3*l + b]; \n"
"        T[l] = D[3*l + a]; \n"
"      } \n"
"	  TensorContract(R,S,T,dv[2],el_uv); //TensorContract(B,B,D,GRAD_TRANSPOSE,dv[2],el_uxv); \n"

"      for (PetscInt l=0; l<3; l++) { \n"
"        Yu[E_times_3+l] += el_uv[l]; \n"
"      } \n"
"    } \n"
"} \n";


PetscErrorCode MFA11SetUp_OpenCL(MatA11MF mf)
{
  PetscErrorCode ierr;
  MFA11OpenCL    ctx;
  PetscReal x1[3],w1[3],B[3][3],D[3][3];
  PetscInt i;
  /* OpenCL-related variables */
  cl_uint num_platforms, platform_index = 0;
  cl_platform_id platform_ids[42];   //no more than 42 platforms supported...
  char buffer[8192];
  cl_device_id device_ids[42];
  cl_uint num_devices, device_index = 0;
  cl_build_status status;

  PetscFunctionBegin;
  if (mf->ctx) PetscFunctionReturn(0);
  ierr = PetscMalloc1(1,&ctx);CHKERRQ(ierr);
  ctx->state = 0;
  mf->ctx = ctx;

  if (sizeof(PetscInt)    != 4) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"OpenCL SpMV kernels only support 32-bit integers!");
  if (sizeof(PetscScalar) != 8) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"OpenCL SpMV kernels only support double precision!");
  if (sizeof(PetscReal)   != 8) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"OpenCL SpMV kernels only support double precision!");

  /************** Initialize OpenCL **************/

  /* Query platform: */
  ierr = clGetPlatformIDs(42,platform_ids,&num_platforms);ERR_CHECK(ierr);

  for (i=0;i<num_platforms;++i) {
	ierr = clGetPlatformInfo(platform_ids[i],CL_PLATFORM_VENDOR,1024 * sizeof(char),buffer,NULL);ERR_CHECK(ierr);
  }

  /* Query devices: */
  ierr = clGetDeviceIDs(platform_ids[platform_index],CL_DEVICE_TYPE_ALL,42,device_ids,&num_devices);ERR_CHECK(ierr);
  for (i=0; i<num_devices; ++i) {
	ierr = clGetDeviceInfo(device_ids[i],CL_DEVICE_NAME,1024 * sizeof(char),buffer,NULL);ERR_CHECK(ierr);
  }

  /* now set up a context containing the selected device: */
  ctx->context = clCreateContext(0,1,&(device_ids[device_index]),NULL,NULL,&ierr);ERR_CHECK(ierr);

  /* create a command queue for the device: */
  ctx->queue = clCreateCommandQueue(ctx->context,device_ids[device_index],0,&ierr);ERR_CHECK(ierr);

  /* create OpenCL program and extract kernel */
  ctx->program = clCreateProgramWithSource(ctx->context,1,&opencl_spmv_kernel_sources,NULL,&ierr);ERR_CHECK(ierr);
  ierr = clBuildProgram(ctx->program,1,&device_ids[device_index],NULL,NULL,NULL);
  if (ierr != CL_SUCCESS) {
	PetscPrintf(PETSC_COMM_WORLD,"Build Scalar: Err = %d\n", ierr);
	ierr = clGetProgramBuildInfo(ctx->program,device_ids[device_index],CL_PROGRAM_BUILD_STATUS,sizeof(cl_build_status),&status,NULL);ERR_CHECK(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Build Status: %d\n", status);
	ierr = clGetProgramBuildInfo(ctx->program,device_ids[device_index],CL_PROGRAM_BUILD_LOG,   sizeof(char)*8192      ,buffer,NULL);ERR_CHECK(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Log: %s\n",buffer);
	PetscPrintf(PETSC_COMM_WORLD,"Sources: %s\n",opencl_spmv_kernel_sources);
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"OpenCL kernel compilation failed!");
  }
  ctx->kernel0 = clCreateKernel(ctx->program,"set_zero",&ierr);ERR_CHECK(ierr);
  ctx->kernel1 = clCreateKernel(ctx->program,"MFStokesWrapper_A11",&ierr);ERR_CHECK(ierr);

  ctx->ufield      = NULL;
  ctx->LA_gcoords  = NULL;
  ctx->gaussdata_w = NULL;
  ctx->elements_per_color = NULL;
  ctx->el_ids_colored     = NULL;
  ctx->elnidx_u    = NULL;
  ctx->Yu          = NULL;

  ierr = PetscDTGaussQuadrature(3,-1,1,x1,w1);CHKERRQ(ierr);
  for (i=0; i<3; i++) {
    B[i][0] = .5*(PetscSqr(x1[i]) - x1[i]);
    B[i][1] = 1 - PetscSqr(x1[i]);
    B[i][2] = .5*(PetscSqr(x1[i]) + x1[i]);
    D[i][0] = x1[i] - .5;
    D[i][1] = -2*x1[i];
    D[i][2] = x1[i] + .5;
  }

  ctx->D = clCreateBuffer(ctx->context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,3 * 3 * sizeof(PetscReal),D,&ierr);ERR_CHECK(ierr);
  ctx->B = clCreateBuffer(ctx->context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,3 * 3 * sizeof(PetscReal),B,&ierr);ERR_CHECK(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode MFA11Destroy_OpenCL(MatA11MF mf)
{
  PetscErrorCode ierr;
  PetscInt       i;
  MFA11OpenCL    ctx;
  
  PetscFunctionBegin;
  ctx = mf->ctx;
  if (!ctx) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"OpenCL MF-SpMV implementation should have a valid context");
  /* Free internal members */
  ierr = clReleaseKernel(ctx->kernel0);ERR_CHECK(ierr);
  ierr = clReleaseKernel(ctx->kernel1);ERR_CHECK(ierr);
  ierr = clReleaseProgram(ctx->program);ERR_CHECK(ierr);
  ierr = clReleaseCommandQueue(ctx->queue);ERR_CHECK(ierr);
  ierr = clReleaseMemObject(ctx->ufield);ERR_CHECK(ierr);
  ierr = clReleaseMemObject(ctx->LA_gcoords);ERR_CHECK(ierr);
  ierr = clReleaseMemObject(ctx->gaussdata_w);ERR_CHECK(ierr);
  for (i=0; i<ctx->element_colors; ++i) {
    ierr = clReleaseMemObject(ctx->el_ids_colored[i]);ERR_CHECK(ierr);
  }
  ierr = PetscFree(ctx->elements_per_color);ERR_CHECK(ierr);
  ierr = PetscFree(ctx->el_ids_colored);ERR_CHECK(ierr);
  ierr = clReleaseMemObject(ctx->elnidx_u);ERR_CHECK(ierr);
  ierr = clReleaseMemObject(ctx->Yu);ERR_CHECK(ierr);
  ierr = clReleaseMemObject(ctx->D);ERR_CHECK(ierr);
  ierr = clReleaseMemObject(ctx->B);ERR_CHECK(ierr);
  ierr = clReleaseContext(ctx->context);ERR_CHECK(ierr);
  /* Free context */
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  mf->ctx = NULL;

	PetscFunctionReturn(0);
}

PetscErrorCode MFStokesWrapper_A11_OpenCL(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[])
{
	PetscErrorCode ierr;
	Vec gcoords;
	const PetscReal *LA_gcoords;
	PetscInt nel,nen_u,e,i,j,k,localsize;
    PetscReal x1[3],w1[3],w[NQP];
	const PetscInt *elnidx_u;
	QPntVolCoefStokes *all_gausspoints;
	const QPntVolCoefStokes *cell_gausspoints;
    PetscReal *gaussdata_host;
    size_t lwsize = 128;
    size_t gwsize;
    MFA11OpenCL openclctx = mf->ctx;

	PetscFunctionBegin;
	/* setup for coords */
	ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);
    ierr = VecGetLocalSize(gcoords,&localsize);CHKERRQ(ierr);

	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);

	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);

    /* Set up OpenCL buffers */
    if (!openclctx->elnidx_u) {
      openclctx->elnidx_u = clCreateBuffer(openclctx->context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,nel * nen_u * sizeof(PetscInt),(void*)elnidx_u,&ierr);ERR_CHECK(ierr);

      /* Assign colors to elements such that there is no overlap in writes to Yu if elements are processed concurrently */
      PetscInt elements_colored = 0;
      PetscInt *element_color;
      PetscInt *Yu_color; // scratchpad

      ierr = PetscMalloc(nel * sizeof(PetscInt), &element_color);CHKERRQ(ierr);
      for (i=0; i<nel; ++i) element_color[i] = -1;
      ierr = PetscMalloc(nel * NQP * sizeof(PetscInt), &Yu_color);CHKERRQ(ierr);
      for (i=0; i<nel * NQP; ++i) Yu_color[i] = -1;

      openclctx->element_colors = 0;
      while (elements_colored < nel) {

        for (i=0; i<nel; ++i) {

          if (element_color[i] >= 0) continue;  /* element already has a color */

          /* Check if element can be colored: No element in Yu has current color */
          PetscInt can_be_colored = 1;
          for (j=0; j<nen_u; ++j) {
            if (Yu_color[elnidx_u[i*nen_u + j]] == openclctx->element_colors) {
              can_be_colored = 0;
              break;
            }
          }

          /* Color element if possible, update Yu indices to current color */
          if (can_be_colored) {
            element_color[i] = openclctx->element_colors;
            for (j=0; j<nen_u; ++j)
              Yu_color[elnidx_u[i*nen_u + j]] = openclctx->element_colors;

            ++elements_colored;
          }
        }

        ++openclctx->element_colors;
      }

      /* Generate OpenCL arrays with coloring information */
      ierr = PetscMalloc(openclctx->element_colors * sizeof(PetscInt),&openclctx->elements_per_color);CHKERRQ(ierr);
      ierr = PetscMalloc(openclctx->element_colors * sizeof(cl_mem),&openclctx->el_ids_colored);CHKERRQ(ierr);

      for (i=0; i<openclctx->element_colors; ++i) {
        /* count elements, collect element indices for this color and copy over to GPU: */
        openclctx->elements_per_color[i] = 0;
        for (j=0; j<nel; ++j) {
          if (element_color[j] == i) {
            Yu_color[openclctx->elements_per_color[i]] = j; /* Reusing Yu_color array here */
            openclctx->elements_per_color[i] += 1;
          }
        }

        openclctx->el_ids_colored[i] = clCreateBuffer(openclctx->context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,openclctx->elements_per_color[i] * sizeof(PetscInt),Yu_color,&ierr);ERR_CHECK(ierr);
      }

      /* clean up */
      ierr = PetscFree(element_color);CHKERRQ(ierr);
      ierr = PetscFree(Yu_color);CHKERRQ(ierr);

    }

    if (!openclctx->ufield) {
      openclctx->ufield = clCreateBuffer(openclctx->context,CL_MEM_READ_ONLY,localsize * sizeof(PetscScalar),NULL,&ierr);ERR_CHECK(ierr);
    }
    ierr = clEnqueueWriteBuffer(openclctx->queue,openclctx->ufield,CL_TRUE,0,localsize * sizeof(PetscScalar),(void*)ufield,0,NULL,NULL);ERR_CHECK(ierr);

    if (!openclctx->LA_gcoords) {
      openclctx->LA_gcoords = clCreateBuffer(openclctx->context,CL_MEM_READ_ONLY,localsize * sizeof(PetscReal),NULL,&ierr);ERR_CHECK(ierr);
    }

    if (!openclctx->gaussdata_w) {
      openclctx->gaussdata_w = clCreateBuffer(openclctx->context,CL_MEM_READ_ONLY,nel * NQP * sizeof(PetscReal),NULL,&ierr);ERR_CHECK(ierr);
    }
    
    if (mf->state != openclctx->state) {
      ierr = clEnqueueWriteBuffer(openclctx->queue,openclctx->LA_gcoords,CL_TRUE,0,localsize * sizeof(PetscReal),(void*)LA_gcoords,0,NULL,NULL);ERR_CHECK(ierr);

      ierr = PetscDTGaussQuadrature(3,-1,1,x1,w1);CHKERRQ(ierr);
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            w[(i*3+j)*3+k] = w1[i] * w1[j] * w1[k];

      ierr = PetscMalloc(nel * NQP * sizeof(PetscReal), &gaussdata_host);CHKERRQ(ierr);
      for (e=0; e<nel; e++) {
        ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,(QPntVolCoefStokes**)&cell_gausspoints);CHKERRQ(ierr);
        for (i=0; i<NQP; i++) gaussdata_host[e*NQP + i] = cell_gausspoints[i].eta * w[i];
      }
      ierr = clEnqueueWriteBuffer(openclctx->queue,openclctx->gaussdata_w,CL_TRUE,0,nel * NQP * sizeof(PetscReal),(void*)gaussdata_host,0,NULL,NULL);ERR_CHECK(ierr);
      ierr = PetscFree(gaussdata_host);CHKERRQ(ierr);

      /* Save new state to avoid unnecessary subsequent copies */
      openclctx->state = mf->state;
    }

    if (!openclctx->Yu) {
      openclctx->Yu = clCreateBuffer(openclctx->context,CL_MEM_READ_WRITE,localsize * sizeof(PetscScalar),NULL,&ierr);ERR_CHECK(ierr);
    }

	ierr = PetscLogEventBegin(MAT_MultMFA11_ker,0,0,0,0);CHKERRQ(ierr);

    /* Launch OpenCL kernel for zeroing Yu */
	ierr = clSetKernelArg(openclctx->kernel0,0,sizeof(cl_mem),  &openclctx->Yu);ERR_CHECK(ierr);
	ierr = clSetKernelArg(openclctx->kernel0,1,sizeof(PetscInt),&localsize);ERR_CHECK(ierr);
    lwsize = 256;
    gwsize = 256 * lwsize;
	ierr = clEnqueueNDRangeKernel(openclctx->queue,openclctx->kernel0,1,NULL,&gwsize,&lwsize,0,NULL,NULL);ERR_CHECK(ierr);

    /* Launch OpenCL kernel for matrix-free SpMV
     *  - inputs: elnidx_u, LA_gcoords, ufield, gaussdata
     *  - output: Yu
     */
    for (i=0; i<openclctx->element_colors; ++i) {
	  ierr = clSetKernelArg(openclctx->kernel1,0,sizeof(PetscInt),&openclctx->elements_per_color[i]);ERR_CHECK(ierr);
      ierr = clSetKernelArg(openclctx->kernel1,1,sizeof(PetscInt),&nen_u                           );ERR_CHECK(ierr);
      ierr = clSetKernelArg(openclctx->kernel1,2,sizeof(cl_mem),  &openclctx->el_ids_colored[i]    );ERR_CHECK(ierr);
      ierr = clSetKernelArg(openclctx->kernel1,3,sizeof(cl_mem),  &openclctx->elnidx_u             );ERR_CHECK(ierr);
      ierr = clSetKernelArg(openclctx->kernel1,4,sizeof(cl_mem),  &openclctx->LA_gcoords           );ERR_CHECK(ierr);
      ierr = clSetKernelArg(openclctx->kernel1,5,sizeof(cl_mem),  &openclctx->ufield               );ERR_CHECK(ierr);
      ierr = clSetKernelArg(openclctx->kernel1,6,sizeof(cl_mem),  &openclctx->gaussdata_w          );ERR_CHECK(ierr);
      ierr = clSetKernelArg(openclctx->kernel1,7,sizeof(cl_mem),  &openclctx->Yu                   );ERR_CHECK(ierr);
      ierr = clSetKernelArg(openclctx->kernel1,8,sizeof(cl_mem),  &openclctx->D                    );ERR_CHECK(ierr);
      ierr = clSetKernelArg(openclctx->kernel1,9,sizeof(cl_mem),  &openclctx->B                    );ERR_CHECK(ierr);
      lwsize = WARPS_PER_BLOCK*32;
      gwsize = ((openclctx->elements_per_color[i]-1)/WARPS_PER_BLOCK + 1) * lwsize;
	  ierr = clEnqueueNDRangeKernel(openclctx->queue,openclctx->kernel1,1,NULL,&gwsize,&lwsize,0,NULL,NULL);ERR_CHECK(ierr);
    }
    ierr = clFinish(openclctx->queue);ERR_CHECK(ierr);
    ierr = PetscLogEventEnd(MAT_MultMFA11_ker,0,0,0,0);CHKERRQ(ierr);

    PetscLogFlops((nel * 9) * 3*NQP*(6+6+6));           /* 9 tensor contractions per element */
    PetscLogFlops(nel*NQP*(14 + 1/* division */ + 27)); /* 1 Jacobi inversion per element */
    PetscLogFlops(nel*NQP*(5*9+6+6+6*9));               /* 1 quadrature action per element */

    /* Read back OpenCL data */
	ierr = PetscLogEventBegin(MAT_MultMFA11_cfr,0,0,0,0);CHKERRQ(ierr);
    ierr = clEnqueueReadBuffer(openclctx->queue,openclctx->Yu,CL_TRUE,0,localsize * sizeof(PetscScalar),&(Yu[0]),0,NULL,NULL);ERR_CHECK(ierr);
    ierr = PetscLogEventEnd(MAT_MultMFA11_cfr,0,0,0,0);CHKERRQ(ierr);

	ierr = VecRestoreArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

