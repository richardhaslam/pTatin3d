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

#include <assert.h>
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

typedef struct _p_MFA11OpenCL *MFA11OpenCL;

struct _p_MFA11OpenCL {
  PetscObjectState state;
};

/* OpenCL error checking */
#define ERROR_CHECKER_CASE(ERRORCODE)  case ERRORCODE: assert("#ERRORCODE");
static void checkError(cl_int err)
{
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
        
      default: assert("Unknown error. Maybe OpenCL SDK not properly installed?");
    }
  }
}

#define ERR_CHECK(err) checkError(err);



#define NQP 27			/* Number of quadrature points per element; must equal Q2_NODES_PER_EL_3D (27) */




static const char * opencl_spmv_kernel_sources =
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable \n"
"#define NQP                 27 \n"
"#define Q2_NODES_PER_EL_3D  27 \n"
"#define PetscReal           double \n"
"#define PetscScalar         double \n"
"#define PetscInt            int \n"

"typedef enum { \n"
"	GRAD, \n"
"	GRAD_TRANSPOSE \n"
"} GradMode; \n"


"void TensorContract(__global PetscReal const *Rf,__global PetscReal const *Sf,__global PetscReal const *Tf,GradMode gmode,PetscReal const x[][NQP],PetscReal y[][NQP]) \n"
"{ \n"
"	PetscReal R[3][3],S[3][3],T[3][3]; \n"
"	PetscReal u[3][NQP],v[3][NQP]; \n"
"    PetscInt i,j,k,l,kj,ji,a,b,c; \n"

"	for (j=0; j<3; j++) { \n"
"		for (i=0; i<3; i++) { \n"
"			R[i][j] = i<3 ? (gmode == GRAD ? Rf[3*i+j] : Rf[3*j + i]) : 0.; \n"
"			S[i][j] = i<3 ? (gmode == GRAD ? Sf[3*i+j] : Sf[3*j + i]) : 0.; \n"
"			T[i][j] = i<3 ? (gmode == GRAD ? Tf[3*i+j] : Tf[3*j + i]) : 0.; \n"
"		} \n"
"	} \n"

"	// u[l,k,j,c] = R[c,i] x[l,k,j,i] \n"
"   for (i=0; i<3; ++i) { \n"
"     for (j=0; j<NQP; ++j) { \n"
"         u[i][j] = 0; \n"
"     } \n"
"   } \n"
"	for (l=0; l<3; l++) { \n"
"		for (kj=0; kj<9; kj++) { \n"
"			for (i=0; i<3; i++) { \n"
"				for (c=0; c<3; c++) { \n"
"					u[l][kj*3+c] += R[c][i] * x[l][kj*3+i]; \n"
"				} \n"
"			} \n"
"		} \n"
"	} \n"

"	// v[l,k,b,c] = S[b,j] u[l,k,j,c] \n"
"   for (i=0; i<3; ++i) { \n"
"      for (j=0; j<NQP; ++j) { \n"
"        v[i][j] = 0; \n"
"      } \n"
"    } \n"
"	for (l=0; l<3; l++) { \n"
"		for (k=0; k<3; k++) { \n"
"			for (j=0; j<3; j++) { \n"
"				for (c=0; c<3; c++) { \n"
"					for (b=0; b<3; b++) { \n"
"						v[l][(k*3+b)*3+c] += S[b][j] * u[l][(k*3+j)*3+c]; \n"
"					} \n"
"				} \n"
"			} \n"
"		} \n"
"	} \n"

"	// y[l,a,b,c] = T[a,k] v[l,k,b,c] \n"
"	for (k=0; k<3; k++) { \n"
"		for (l=0; l<3; l++) { \n"
"			for (a=0; a<3; a++) { \n"
"				for (ji=0; ji<9; ji++) { \n"
"					y[l][a*9+ji] += T[a][k] * v[l][k*9+ji]; \n"
"				} \n"
"			} \n"
"		} \n"
"	} \n"
"} \n"

" void JacobianInvert(PetscScalar dx[3][3][NQP],PetscScalar dxdet[NQP]) \n"
"{ \n"
"	PetscInt i,j,k; \n"

"	for (i=0; i<NQP; i++) { \n"
"		PetscScalar a[3][3]; \n"
"		PetscScalar b0,b3,b6,det,idet; \n"
"		for (j=0; j<3; j++) { \n"
"			for (k=0; k<3; k++) { \n"
"				a[j][k] = dx[j][k][i]; \n"
"			} \n"
"		} \n"
"		b0 =  (a[1][1]*a[2][2] - a[2][1]*a[1][2]); \n"
"		b3 = -(a[1][0]*a[2][2] - a[2][0]*a[1][2]); \n"
"		b6 =  (a[1][0]*a[2][1] - a[2][0]*a[1][1]); \n"
"		det = a[0][0]*b0 + a[0][1]*b3 + a[0][2]*b6; \n"
"		idet = 1.0 / det; \n"
"		dx[0][0][i] =  idet*b0; \n"
"		dx[0][1][i] = -idet*(a[0][1]*a[2][2] - a[2][1]*a[0][2]); \n"
"		dx[0][2][i] =  idet*(a[0][1]*a[1][2] - a[1][1]*a[0][2]); \n"
"		dx[1][0][i] =  idet*b3; \n"
"		dx[1][1][i] =  idet*(a[0][0]*a[2][2] - a[2][0]*a[0][2]); \n"
"		dx[1][2][i] = -idet*(a[0][0]*a[1][2] - a[1][0]*a[0][2]); \n"
"		dx[2][0][i] =  idet*b6; \n"
"		dx[2][1][i] = -idet*(a[0][0]*a[2][1] - a[2][0]*a[0][1]); \n"
"		dx[2][2][i] =  idet*(a[0][0]*a[1][1] - a[1][0]*a[0][1]); \n"
"		dxdet[i] =  det; \n"
"	} \n"
"} \n"

"void QuadratureAction(__global const PetscScalar *gaussdata_eta, \n"
"				       PetscScalar const dx[3][3][Q2_NODES_PER_EL_3D], \n"
"				       PetscScalar const dxdet[Q2_NODES_PER_EL_3D], \n"
"				       __global PetscReal const *w, \n"
"				       PetscScalar const du[3][3][Q2_NODES_PER_EL_3D], \n"
"				       PetscScalar dv[3][3][Q2_NODES_PER_EL_3D]) \n"
"{ \n"
"	PetscInt i,l,k; \n"

"	for (i=0; i<NQP; i++) { \n"
"		PetscScalar Du[6],Dv[6]; // Symmetric gradient with respect to physical coordinates, xx, yy, zz, xy+yx, xz+zx, yz+zy \n"

"		PetscScalar dux[3][3]; \n"
"		for (l=0; l<3; l++) { // fields \n"
"			for (k=0; k<3; k++) { // directions \n"
"				dux[k][l] = du[0][l][i] * dx[k][0][i] + du[1][l][i] * dx[k][1][i] + du[2][l][i] * dx[k][2][i]; \n"
"			} \n"
"		} \n"
"		Du[0] = dux[0][0]; \n"
"		Du[1] = dux[1][1]; \n"
"		Du[2] = dux[2][2]; \n"
"		Du[3] = 0.5*(dux[0][1] + dux[1][0]); \n"
"		Du[4] = 0.5*(dux[0][2] + dux[2][0]); \n"
"		Du[5] = 0.5*(dux[1][2] + dux[2][1]); \n"

"		for (k=0; k<6; k++) { // Stress is coefficient of test function  \n"
"			Dv[k] = 2 * gaussdata_eta[i] * Du[k]; \n"
"		} \n"

"		PetscScalar dvx[3][3]; \n"
"		dvx[0][0] = Dv[0]; \n"
"		dvx[0][1] = Dv[3]; \n"
"		dvx[0][2] = Dv[4]; \n"
"		dvx[1][0] = Dv[3]; \n"
"		dvx[1][1] = Dv[1]; \n"
"		dvx[1][2] = Dv[5]; \n"
"		dvx[2][0] = Dv[4]; \n"
"		dvx[2][1] = Dv[5]; \n"
"		dvx[2][2] = Dv[2]; \n"

"		for (l=0; l<3; l++) { // fields \n"
"			for (k=0; k<3; k++) { // directions \n"
"				dv[k][l][i] = w[i] * dxdet[i] * (dvx[0][l] * dx[0][k][i] + dvx[1][l] * dx[1][k][i] + dvx[2][l] * dx[2][k][i]); \n"
"			} \n"
"		} \n"
"	} \n"
"} \n"

"__kernel void MFStokesWrapper_A11(PetscInt nel,PetscInt nen_u,__global PetscInt const *elnidx_u,__global PetscReal const *LA_gcoords, \n"
"                                  __global PetscScalar const *ufield,__global PetscReal const *gaussdata,__global PetscScalar *Yu_premerge, \n"
"                                  __global PetscReal const *D,__global PetscReal const *B,__global PetscReal const *w) \n"
"{ \n"
"	PetscScalar elu[3][Q2_NODES_PER_EL_3D]={},elx[3][Q2_NODES_PER_EL_3D]={},elv[3][Q2_NODES_PER_EL_3D]; \n"
"	PetscScalar dx[3][3][NQP],dxdet[NQP],du[3][3][NQP],dv[3][3][NQP]; \n"
"	PetscInt i,j,k,l; \n"
"	PetscInt elidx = get_global_id(0); \n"

"    if (elidx >= nel) \n"
"      return; \n"

"	for (i=0; i<Q2_NODES_PER_EL_3D; i++) { \n"
"		PetscInt E = elnidx_u[nen_u*elidx+i]; \n"
"		for (l=0; l<3; l++) { \n"
"			elx[l][i] = LA_gcoords[3*E+l]; \n"
"			elu[l][i] = ufield[3*E+l]; \n"
"            elv[l][i] = 0; \n"
"		} \n"
"	} \n"

"	for (i=0; i<3; i++) { \n"
"      for (j=0; j<3; j++) { \n"
"        for (k=0; k<NQP; k++) { \n"
"          dx[i][j][k] = 0; \n"
"          du[i][j][k] = 0; \n"
"          dv[i][j][k] = 0; \n"
"        } \n"
"      } \n"
"    } \n"

"	TensorContract(D,B,B,GRAD,elx,dx[0]); \n"
"	TensorContract(B,D,B,GRAD,elx,dx[1]); \n"
"	TensorContract(B,B,D,GRAD,elx,dx[2]); \n"

"	JacobianInvert(dx,dxdet); \n"

"	TensorContract(D,B,B,GRAD,elu,du[0]); \n"
"	TensorContract(B,D,B,GRAD,elu,du[1]); \n"
"	TensorContract(B,B,D,GRAD,elu,du[2]); \n"

"	QuadratureAction(gaussdata + elidx*NQP,dx,dxdet,w,du,dv); \n"

"	TensorContract(D,B,B,GRAD_TRANSPOSE,dv[0],elv); \n"
"	TensorContract(B,D,B,GRAD_TRANSPOSE,dv[1],elv); \n"
"	TensorContract(B,B,D,GRAD_TRANSPOSE,dv[2],elv); \n"

"    // Avoid worries about race conditions by applying the merge to Yu on the host! \n"
"	for (i=0; i<NQP; i++) { \n"
"		for (l=0; l<3; l++) { \n"
"			(Yu_premerge + 3*NQP*elidx)[3*i+l] = elv[l][i]; \n"
"		} \n"
"	} \n"
"} \n";


#undef __FUNCT__
#define __FUNCT__ "MFA11SetUp_OpenCL"
PetscErrorCode MFA11SetUp_OpenCL(MatA11MF mf)
{
  PetscErrorCode ierr;
  MFA11OpenCL    ctx;

	PetscFunctionBegin;
  
  if (mf->ctx) PetscFunctionReturn(0);
  
  ierr = PetscMalloc1(1,&ctx);CHKERRQ(ierr);
  ctx->state = 0;
	mf->ctx = ctx;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MFA11Destroy_OpenCL"
PetscErrorCode MFA11Destroy_OpenCL(MatA11MF mf)
{
  PetscErrorCode ierr;
  MFA11OpenCL    ctx;
  
  PetscFunctionBegin;
  ctx = mf->ctx;
  if (!ctx) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"OpenCL MF-SpMV implementation should have a valid context");
  /* Free internal members */
  /* Free context */
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  mf->ctx = NULL;
  
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MFStokesWrapper_A11_OpenCL"
PetscErrorCode MFStokesWrapper_A11_OpenCL(MatA11MF mf,Quadrature volQ,DM dau,PetscScalar ufield[],PetscScalar Yu[])
{
	PetscErrorCode ierr;
	DM cda;
	Vec gcoords;
	const PetscReal *LA_gcoords;
	PetscInt nel,nen_u,e,i,j,k,l,localsize;
	const PetscInt *elnidx_u;
	QPntVolCoefStokes *all_gausspoints;
	const QPntVolCoefStokes *cell_gausspoints;
	PetscReal x1[3],w1[3],B[3][3],D[3][3],w[NQP];
    /* OpenCL-related variables */
    static PetscBool opencl_init = PETSC_FALSE;
	static cl_uint num_platforms, platform_index = 0;
	static cl_platform_id platform_ids[42];   //no more than 42 platforms supported...
	char buffer[8192]; 
	static cl_device_id device_ids[42];
	static cl_uint num_devices, device_index = 0;
	static cl_context my_context;
	static cl_command_queue queue;
	static cl_program my_program;
	static cl_build_status status;
	static cl_kernel my_kernel;
    cl_mem elnidx_u_opencl,ufield_opencl,LA_gcoords_opencl,gaussdata_opencl,Yu_premerge_opencl,D_opencl,B_opencl,w_opencl;
    PetscReal *gaussdata_host;
    PetscScalar *Yu_premerge;
    size_t lwsize = 256;
    size_t gwsize = 256 * lwsize;


	PetscFunctionBegin;
    if (opencl_init == PETSC_FALSE) {

      if (sizeof(PetscInt)    != 4) assert("OpenCL SpMV kernels only support 32-bit integers!");
      if (sizeof(PetscScalar) != 8) assert("OpenCL SpMV kernels only support double precision!");
      if (sizeof(PetscReal)   != 8) assert("OpenCL SpMV kernels only support double precision!");

	  /************** Initialize OpenCL **************/
		
	  /*
	  * Query platform:
	  */
	  ierr = clGetPlatformIDs(42,platform_ids,&num_platforms);ERR_CHECK(ierr);

	  printf("# OpenCL Platforms found: %d\n",num_platforms);
	  for (i=0;i<num_platforms;++i)
	  {
		ierr = clGetPlatformInfo(platform_ids[i],CL_PLATFORM_VENDOR,1024 * sizeof(char),buffer,NULL);ERR_CHECK(ierr);
		printf("# (%d) %s\n", i, buffer);
	  }

	  /*
	   * Query devices:
	   */
	  ierr = clGetDeviceIDs(platform_ids[platform_index],CL_DEVICE_TYPE_ALL,42,device_ids,&num_devices);ERR_CHECK(ierr);
	  printf("# Devices found: %d\n", num_devices);
	  for (i=0; i<num_devices; ++i)
	  {
		ierr = clGetDeviceInfo(device_ids[i],CL_DEVICE_NAME,1024 * sizeof(char),buffer,NULL);ERR_CHECK(ierr);
		printf("# (%d) %s\n", i, buffer);
	  }

	  /* now set up a context containing the selected device: */
	  my_context = clCreateContext(0,1,&(device_ids[device_index]),NULL,NULL,&ierr);ERR_CHECK(ierr);
	   
	  /* create a command queue for the device: */
	  queue = clCreateCommandQueue(my_context,device_ids[device_index],0,&ierr);ERR_CHECK(ierr);

	  /* create OpenCL program and extract kernel */
	  my_program = clCreateProgramWithSource(my_context,1,&opencl_spmv_kernel_sources,NULL,&ierr);ERR_CHECK(ierr);
	  ierr = clBuildProgram(my_program,1,&device_ids[device_index],NULL,NULL,NULL);
	  if (ierr != CL_SUCCESS)
	  {
		printf("Build Scalar: Err = %d\n", ierr);
		ierr = clGetProgramBuildInfo(my_program,device_ids[device_index],CL_PROGRAM_BUILD_STATUS,sizeof(cl_build_status),&status,NULL);ERR_CHECK(ierr);
		printf("Build Status: %d\n", status);
		ierr = clGetProgramBuildInfo(my_program,device_ids[device_index],CL_PROGRAM_BUILD_LOG,   sizeof(char)*8192      ,buffer,NULL);ERR_CHECK(ierr);
		printf("Log: %s\n",buffer);
		printf("Sources: %s\n",opencl_spmv_kernel_sources);
	  }
	  my_kernel = clCreateKernel(my_program,"MFStokesWrapper_A11",&ierr);ERR_CHECK(ierr);

      /* OpenCL initialization completed */
      opencl_init = PETSC_TRUE;
    }

	ierr = PetscDTGaussQuadrature(3,-1,1,x1,w1);CHKERRQ(ierr);
	for (i=0; i<3; i++) {
		B[i][0] = .5*(PetscSqr(x1[i]) - x1[i]);
		B[i][1] = 1 - PetscSqr(x1[i]);
		B[i][2] = .5*(PetscSqr(x1[i]) + x1[i]);
		D[i][0] = x1[i] - .5;
		D[i][1] = -2*x1[i];
		D[i][2] = x1[i] + .5;
	}
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			for (k=0; k<3; k++) {
				w[(i*3+j)*3+k] = w1[i] * w1[j] * w1[k];}}}

	/* setup for coords */
	ierr = DMGetCoordinateDM(dau,&cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dau,&gcoords);CHKERRQ(ierr);
	ierr = VecGetArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);
    ierr = VecGetLocalSize(gcoords,&localsize);CHKERRQ(ierr);

	ierr = DMDAGetElements_pTatinQ2P1(dau,&nel,&nen_u,&elnidx_u);CHKERRQ(ierr);

	ierr = VolumeQuadratureGetAllCellData_Stokes(volQ,&all_gausspoints);CHKERRQ(ierr);

    /* Set up OpenCL buffers */
      elnidx_u_opencl = clCreateBuffer(my_context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nel * nen_u * sizeof(PetscInt),(void*)elnidx_u,  &ierr);ERR_CHECK(ierr);
        ufield_opencl = clCreateBuffer(my_context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,localsize * sizeof(PetscScalar),(void*)ufield,    &ierr);ERR_CHECK(ierr);
    LA_gcoords_opencl = clCreateBuffer(my_context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,localsize * sizeof(PetscReal  ),(void*)LA_gcoords,&ierr);ERR_CHECK(ierr);
    
    ierr = PetscMalloc(nel * NQP * sizeof(PetscReal), &gaussdata_host);CHKERRQ(ierr);
    for (e=0; e<nel; e++) {
      ierr = VolumeQuadratureGetCellData_Stokes(volQ,all_gausspoints,e,(QPntVolCoefStokes**)&cell_gausspoints);CHKERRQ(ierr);
      for (i=0; i<NQP; i++) gaussdata_host[e*NQP + i] = cell_gausspoints[i].eta;
    }
    gaussdata_opencl = clCreateBuffer(my_context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,nel * NQP * sizeof(PetscReal),&(gaussdata_host[0]),&ierr);ERR_CHECK(ierr);

    ierr = PetscMalloc(3 * nel * nen_u * sizeof(PetscScalar),&Yu_premerge);CHKERRQ(ierr);
    Yu_premerge_opencl = clCreateBuffer(my_context,CL_MEM_READ_WRITE, 3 * nel * nen_u * sizeof(PetscScalar), NULL, &ierr);ERR_CHECK(ierr);

    D_opencl = clCreateBuffer(my_context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,     3 * 3 * sizeof(PetscReal),&(D[0]),&ierr);ERR_CHECK(ierr);
    B_opencl = clCreateBuffer(my_context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,     3 * 3 * sizeof(PetscReal),&(B[0]),&ierr);ERR_CHECK(ierr);
    w_opencl = clCreateBuffer(my_context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 3 * 3 * 3 * sizeof(PetscReal),&(w[0]),&ierr);ERR_CHECK(ierr);

    /* Launch OpenCL kernel
     *  - inputs: elnidx_u, LA_gcoords, ufield, gaussdata
     *  - output: Yu_premerge
     */
	ierr = clSetKernelArg(my_kernel,0,sizeof(PetscInt),&nel               );ERR_CHECK(ierr);
	ierr = clSetKernelArg(my_kernel,1,sizeof(PetscInt),&nen_u             );ERR_CHECK(ierr);
	ierr = clSetKernelArg(my_kernel,2,sizeof(cl_mem),  &elnidx_u_opencl   );ERR_CHECK(ierr);
	ierr = clSetKernelArg(my_kernel,3,sizeof(cl_mem),  &LA_gcoords_opencl );ERR_CHECK(ierr);
	ierr = clSetKernelArg(my_kernel,4,sizeof(cl_mem),  &ufield_opencl     );ERR_CHECK(ierr);
	ierr = clSetKernelArg(my_kernel,5,sizeof(cl_mem),  &gaussdata_opencl  );ERR_CHECK(ierr);
	ierr = clSetKernelArg(my_kernel,6,sizeof(cl_mem),  &Yu_premerge_opencl);ERR_CHECK(ierr);
	ierr = clSetKernelArg(my_kernel,7,sizeof(cl_mem),  &D_opencl          );ERR_CHECK(ierr);
	ierr = clSetKernelArg(my_kernel,8,sizeof(cl_mem),  &B_opencl          );ERR_CHECK(ierr);
	ierr = clSetKernelArg(my_kernel,9,sizeof(cl_mem),  &w_opencl          );ERR_CHECK(ierr);
	ierr = clEnqueueNDRangeKernel(queue,my_kernel,1,NULL,&gwsize,&lwsize,0,NULL,NULL);ERR_CHECK(ierr);

    PetscLogFlops((nel * 9) * 3*NQP*(6+6+6));           /* 9 tensor contractions per element */
    PetscLogFlops(nel*NQP*(14 + 1/* division */ + 27)); /* 1 Jacobi inversion per element */
    PetscLogFlops((nel * 9) * 3*NQP*(6+6+6));           /* 1 quadrature action per element */

    /* Read back OpenCL data */
    ierr = clEnqueueReadBuffer(queue,Yu_premerge_opencl,CL_TRUE,0,3 * nel * NQP * sizeof(PetscScalar),&(Yu_premerge[0]),0,NULL,NULL);ERR_CHECK(ierr);

    for (e=0; e<nel; e++) {
      for (i=0; i<NQP; i++) {
		PetscInt E = elnidx_u[nen_u*e+i];
		for (l=0; l<3; l++) {
			Yu[3*E+l] += (Yu_premerge + 3*NQP*e)[3*i+l];
		}
      }
	}

	ierr = VecRestoreArrayRead(gcoords,&LA_gcoords);CHKERRQ(ierr);

    /* clean up */
    ierr = clReleaseMemObject(elnidx_u_opencl);ERR_CHECK(ierr);
    ierr = clReleaseMemObject(ufield_opencl);ERR_CHECK(ierr);
    ierr = clReleaseMemObject(LA_gcoords_opencl);ERR_CHECK(ierr);
    ierr = PetscFree(gaussdata_host);CHKERRQ(ierr);
    ierr = clReleaseMemObject(gaussdata_opencl);ERR_CHECK(ierr);
    ierr = PetscFree(Yu_premerge);CHKERRQ(ierr);
    ierr = clReleaseMemObject(Yu_premerge_opencl);ERR_CHECK(ierr);
    ierr = clReleaseMemObject(D_opencl);ERR_CHECK(ierr);
    ierr = clReleaseMemObject(B_opencl);ERR_CHECK(ierr);
    ierr = clReleaseMemObject(w_opencl);ERR_CHECK(ierr);

	PetscFunctionReturn(0);
}

