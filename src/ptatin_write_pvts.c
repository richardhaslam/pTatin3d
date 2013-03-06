/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012, 
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Geophysical Fluid Dynamics, 
 **        Department of Earth Sciences,
 **        ETH Zürich,
 **        Sonneggstrasse 5,
 **        CH-8092 Zurich,
 **        Switzerland
 **
 **    Project:       pTatin3d
 **    Filename:      ptatin_write_pvts.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published by
 **    the Free Software Foundation, either version 3 of the License, or
 **    (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **    GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d.  If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    $Id$
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@*/

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <petsc.h>
#include <petscvec.h>
#include <petscdm.h>

#include "ptatin3d.h"
#include "ptatin3d_defs.h"
#include "private/ptatin_impl.h"

#include "dmda_update_coords.h"
#include "dmda_view_petscvtk.h"
#include "dmda_checkpoint.h"
#include "mesh_deformation.h"
#include "ptatin3d_stokes.h"
#include "output_paraview.h"

//#define HAVE_STRLCAT

#undef __FUNCT__  
#define __FUNCT__ "_strlcat"
PetscErrorCode _strlcat(char orig[],char append[],size_t L)
{
#ifdef HAVE_STRLCAT
	strlcat(orig,append,L);
#else	
	char *new;
	size_t l1,l2;
	
	l1 = strlen(orig);
	l2 = strlen(append);
	
	if (l1+l2>=L) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Input string is not long enough");
	}
	
	asprintf(&new,"%s%s",orig,append);
	strcpy(orig,new);
	free(new);
#endif	
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "pTatinVecFieldWrite"
PetscErrorCode pTatinVecFieldWrite(Vec x,const char name[],PetscBool zip_file)
{
	char fieldname[PETSC_MAX_PATH_LEN];
	PetscViewer viewer;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	if (zip_file) {
		sprintf(fieldname,"%s.gz",name);
	} else {
		sprintf(fieldname,"%s",name);
	}

  ierr = PetscViewerCreate(((PetscObject)x)->comm,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
	ierr = PetscViewerBinarySetMPIIO(viewer);CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(viewer,fieldname);CHKERRQ(ierr);
	
	ierr = VecView(x,viewer);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "pTatinVecFieldRead"
PetscErrorCode pTatinVecFieldRead(const char name[],PetscBool zip_file,Vec x)
{
	char fieldname[PETSC_MAX_PATH_LEN];
	PetscViewer viewer;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	if (zip_file) {
		sprintf(fieldname,"%s.gz",name);
	} else {
		sprintf(fieldname,"%s",name);
	}

  ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_READ);CHKERRQ(ierr);
	ierr = PetscViewerBinarySetMPIIO(viewer);CHKERRQ(ierr);
	ierr = PetscViewerFileSetName(viewer,fieldname);CHKERRQ(ierr);
	
//	ierr = VecCreate(((PetscObject)viewer)->comm,x);CHKERRQ(ierr);
	ierr = VecLoad(x,viewer);CHKERRQ(ierr);
//	ierr = VecSetBlockSize(*x,6);CHKERRQ(ierr);
	ierr = PetscOptionsClearValue("-vecload_block_size"); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "test_pTatinVecFieldWrite"
PetscErrorCode test_pTatinVecFieldWrite(void) 
{
	DM  da;
	PetscInt nx,ny,nz;
	Vec x,coords;
	PetscViewer v;
	PetscReal val;
	PetscScalar max;
	PetscReal x0,y0,z0,x1,y1,z1;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	/* create the da */
	nx = ny = nz = 10;
	PetscOptionsGetInt( PETSC_NULL, "-mx", &nx, 0 );
	PetscOptionsGetInt( PETSC_NULL, "-my", &ny, 0 );
	PetscOptionsGetInt( PETSC_NULL, "-mz", &nz, 0 );
	
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,nx,ny,nz, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, 6,1, 0,0,0,&da);CHKERRQ(ierr);
	
	x0 = y0 = z0 = -1.0;
	x1 = y1 = z1 = 1.0;
	ierr = DMDASetUniformCoordinates(da, x0,x1, y0,y1, z0,z1);CHKERRQ(ierr);
	
	ierr = MeshDeformation_GaussianBump_YMAX(da);CHKERRQ(ierr);
	
	/* create a field */
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
	ierr = VecSetRandom( x, PETSC_NULL );CHKERRQ(ierr);
	ierr = VecNorm( x, NORM_1, &val );	PetscPrintf( PETSC_COMM_WORLD, "|x| = %1.5e \n", val );CHKERRQ(ierr);
	ierr = VecNorm( x, NORM_2, &val ); PetscPrintf( PETSC_COMM_WORLD, "|x|_2 = %1.5e \n", val );CHKERRQ(ierr);
	ierr = VecMin( x, 0, &max ); PetscPrintf( PETSC_COMM_WORLD, "min(x) = %1.5e \n", max );CHKERRQ(ierr);
	ierr = VecMax( x, 0, &max ); PetscPrintf( PETSC_COMM_WORLD, "max(x) = %1.5e \n", max );CHKERRQ(ierr);
	
	/* dump field to vtk */
	ierr = DMDAViewPetscVTK(da, x, "dmda_write.vtk");CHKERRQ(ierr);
	
	/* dump field to disk */
	ierr = pTatinVecFieldWrite(x,"mesh_ufield.ptatinvec",PETSC_TRUE);CHKERRQ(ierr);
	 
	/* dump dm to disk */
	ierr = DMDAPackDataToFile( da, "mesh.ptatindmda" );CHKERRQ(ierr);
	
	ierr = DMDestroy(&da);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "test_pTatinVecFieldLoad"
PetscErrorCode test_pTatinVecFieldLoad( void ) 
{
	DM  da;
	Vec x,coords;
	PetscReal val;
	PetscScalar max;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = DMDACreateFromPackDataToFile(PETSC_COMM_WORLD,"mesh.ptatindmda",&da);CHKERRQ(ierr);
	ierr = DMView(da,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
	
	ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
//	ierr = pTatinVecFieldRead("mesh_ufield.ptatinvec",PETSC_FALSE,x);CHKERRQ(ierr);
	ierr = pTatinVecFieldRead("mesh_ufield.ptatinvec",PETSC_TRUE,x);CHKERRQ(ierr);

	ierr = VecNorm( x, NORM_1, &val );	PetscPrintf( PETSC_COMM_WORLD, "|x| = %1.5e \n", val );CHKERRQ(ierr);
	ierr = VecNorm( x, NORM_2, &val ); PetscPrintf( PETSC_COMM_WORLD, "|x|_2 = %1.5e \n", val );CHKERRQ(ierr);
	ierr = VecMin( x, 0, &max ); PetscPrintf( PETSC_COMM_WORLD, "min(x) = %1.5e \n", max );CHKERRQ(ierr);
	ierr = VecMax( x, 0, &max ); PetscPrintf( PETSC_COMM_WORLD, "max(x) = %1.5e \n", max );CHKERRQ(ierr);
	
	/* dump field to vtk */
	ierr = DMDAViewPetscVTK(da,x,"dmda_write_2.vtk");CHKERRQ(ierr);
	
	ierr = DMDestroy(&da);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "test_DMDACheckPoint"
PetscErrorCode test_DMDACheckPoint(void)
{
	PetscErrorCode ierr;
	PetscBool restart, checkpoint,flg;
	
	PetscFunctionBegin;
	checkpoint = PETSC_FALSE;
	PetscOptionsGetBool( PETSC_NULL, "-checkpoint", &checkpoint, &flg );
	if( checkpoint == PETSC_TRUE ) {
		ierr = test_pTatinVecFieldWrite();CHKERRQ(ierr);
	}
	
	restart = PETSC_FALSE;
	PetscOptionsGetBool( PETSC_NULL, "-restart", &restart, &flg );
	if( restart == PETSC_TRUE ) {
		ierr = test_pTatinVecFieldLoad();CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompOutput_StokesRawVelocityPressure"
PetscErrorCode PhysCompOutput_StokesRawVelocityPressure(PhysCompStokes ctx,Vec X,PetscBool zip_file,const char outputpath[],const char prefix[])
{
	PetscLogDouble t0,t1,tl,tg;
	char start[PETSC_MAX_PATH_LEN];
	char f1[PETSC_MAX_PATH_LEN];
	char f3[PETSC_MAX_PATH_LEN];
	Vec Xu,Xp;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	PetscGetTime(&t0);
	/* dav,dap */
	if (prefix) {
		sprintf(f1,"%s/pt3d_stokes.dmda-velocity%s",outputpath,prefix);
		sprintf(f3,"%s/pt3d_stokes.dmda-pressure%s",outputpath,prefix);
	} else {
		sprintf(f1,"%s/pt3d_stokes.dmda-velocity",outputpath);
		sprintf(f3,"%s/pt3d_stokes.dmda-pressure",outputpath);
	}
	PetscPrintf(PETSC_COMM_WORLD,"  writing %s \n", f1 );
	PetscPrintf(PETSC_COMM_WORLD,"  writing %s \n", f3 );
	
	ierr = PhysCompSaveMesh_Stokes3d(ctx,f1,f3,PETSC_NULL);CHKERRQ(ierr);
	
	
	ierr = DMCompositeGetAccess(ctx->stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);

	if (prefix) { sprintf(f1,"%s/pt3d_stokes.dmda-Xu%s",outputpath,prefix); } 
	else {        sprintf(f1,"%s/pt3d_stokes.dmda-Xu",outputpath);           }
	ierr = pTatinVecFieldWrite(Xu,f1,zip_file);CHKERRQ(ierr);

	if (prefix) { sprintf(f1,"%s/pt3d_stokes.dmda-Xp%s",outputpath,prefix); } 
	else {        sprintf(f1,"%s/pt3d_stokes.dmda-Xp",outputpath);           }
	ierr = pTatinVecFieldWrite(Xp,f1,zip_file);CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(ctx->stokes_pack,X,&Xu,&Xp);CHKERRQ(ierr);
	PetscGetTime(&t1);
	tl = t1-t0;
	ierr = MPI_Allreduce(&tl,&tg,1,MPIU_REAL,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);	

	PetscPrintf(PETSC_COMM_WORLD,"%s() -> pt3d_stokes.{dmda.v,dmda,p,dmda.Xu,dmda.Xp}: CPU time %1.2e (sec) \n", __FUNCT__,prefix,tg);
	
	PetscFunctionReturn(0);
}

/* minimal loading from checkpoint file for viz of u,v,w,p fields */
#undef __FUNCT__  
#define __FUNCT__ "PhysCompStokesLoad_DM"
PetscErrorCode PhysCompStokesLoad_DM(const char vname[],const char pname[],PhysCompStokes *ctx)
{
	PhysCompStokes stokes;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = PhysCompCreate_Stokes(&stokes);CHKERRQ(ierr);	
	ierr = PhysCompLoadMesh_Stokes3d(stokes,vname,pname);CHKERRQ(ierr);
	*ctx = stokes;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompStokesLoad_VP"
PetscErrorCode PhysCompStokesLoad_VP(PhysCompStokes ctx,const char vname[],const char pname[],PetscBool zip_file,Vec *VP)
{
	Vec            X,Xv,Xp;
	DM             dav,dap;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = DMCreateGlobalVector(ctx->stokes_pack,&X);CHKERRQ(ierr);
	
	ierr = DMCompositeGetAccess(ctx->stokes_pack,X,&Xv,&Xp);CHKERRQ(ierr);
	
	/* load from file */
	ierr = pTatinVecFieldRead(vname,zip_file,Xv);CHKERRQ(ierr);
	ierr = pTatinVecFieldRead(pname,zip_file,Xp);CHKERRQ(ierr);
	
	ierr = DMCompositeRestoreAccess(ctx->stokes_pack,X,&Xv,&Xp);CHKERRQ(ierr);
	
	*VP = X;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompStokesLoad_X"
PetscErrorCode PhysCompStokesLoad_X(PhysCompStokes ctx,const char xname[],PetscBool zip_file,Vec *VP)
{
	Vec            X;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = DMCreateGlobalVector(ctx->stokes_pack,&X);CHKERRQ(ierr);
	
	/* load from file */
	ierr = pTatinVecFieldRead(xname,zip_file,X);CHKERRQ(ierr);
	
	*VP = X;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PhysCompStokesWrite_DM_X"
PetscErrorCode PhysCompStokesWrite_DM_X(PhysCompStokes ctx,Vec VP,const char outputpath[],const char name[])
{
	PetscErrorCode ierr;

	PetscFunctionBegin;
	
	ierr = pTatinOutputParaViewMeshVelocityPressure(ctx->stokes_pack,VP,outputpath,name);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/*
 I haven't figured out the reason why, but when I write the entire vector X into a single file and read it back,
 the ordering is completely wrong. Serial write - serial read works, however any parallel write/read's, even with
 the same number of processors during write/read produces output incorrectly ordered when plotted on the DMDA.
 For this reason, I write out the u,p vectors seperately, each in their own file.
*/
#undef __FUNCT__  
#define __FUNCT__ "_test_load_and_writevts_from_checkpoint_file"
PetscErrorCode _test_load_and_writevts_from_checkpoint_file(void)
{
	PhysCompStokes ctx;
	Vec            X;
	char           vname[PETSC_MAX_PATH_LEN],pname[PETSC_MAX_PATH_LEN],xname[PETSC_MAX_PATH_LEN];
	char           outputpath[PETSC_MAX_PATH_LEN];
	char           prefix[PETSC_MAX_PATH_LEN];
	char           suffix[PETSC_MAX_PATH_LEN];
	char           outfilename[PETSC_MAX_PATH_LEN];
	PetscBool      flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;

	flg = PETSC_FALSE;
	ierr = PetscOptionsGetString(PETSC_NULL,"-output_path",outputpath,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (flg == PETSC_FALSE) { 
		sprintf(outputpath,".");
	}

	flg = PETSC_FALSE;
	ierr = PetscOptionsGetString(PETSC_NULL,"-file_prefix",prefix,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (flg == PETSC_FALSE) { 
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Required to know the prefix of the checkpointed files");
	}
	
	sprintf(vname,"%s/%s.dmda-velocity",outputpath,prefix);
	sprintf(pname,"%s/%s.dmda-pressure",outputpath,prefix);
	sprintf(xname,"%s/%s.dmda-X",outputpath,prefix);

	flg = PETSC_FALSE;
	ierr = PetscOptionsGetString(PETSC_NULL,"-file_suffix",suffix,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (flg) { 
		sprintf(vname,"%s%s",vname,suffix);
		sprintf(pname,"%s%s",pname,suffix);
		sprintf(xname,"%s%s",xname,suffix);
	}
	PetscPrintf(PETSC_COMM_WORLD,"Reading files: (dmda-vel) %s\n\t\t(dmda-p) %s\n\t\t(vec-X) %s\n", vname,pname,xname);
	
	ierr = PhysCompStokesLoad_DM(vname,pname,&ctx);CHKERRQ(ierr);
	ierr = PhysCompStokesLoad_X(ctx,xname,PETSC_TRUE,&X);CHKERRQ(ierr);
	
	sprintf(outfilename,"%s",prefix);
	if (flg) { 
		sprintf(outfilename,"%s%s",outfilename,suffix);
	}
	sprintf(outfilename,"%s-vp",outfilename);
	PetscPrintf(PETSC_COMM_WORLD,"Writing file: %s/%s \n", outputpath,outfilename);
	ierr = PhysCompStokesWrite_DM_X(ctx,X,outputpath,outfilename);CHKERRQ(ierr);
	
	ierr = PhysCompDestroy_Stokes(&ctx);CHKERRQ(ierr);
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "test_LoadStokesFromCheckpoint_WriteToVTS"
PetscErrorCode test_LoadStokesFromCheckpoint_WriteToVTS(void)
{
	PhysCompStokes ctx;
	Vec            X;
	char           vname[PETSC_MAX_PATH_LEN],pname[PETSC_MAX_PATH_LEN],xuname[PETSC_MAX_PATH_LEN],xpname[PETSC_MAX_PATH_LEN];
	char           outputpath[PETSC_MAX_PATH_LEN];
	char           prefix[PETSC_MAX_PATH_LEN];
	char           suffix[PETSC_MAX_PATH_LEN];
	char           outfilename[PETSC_MAX_PATH_LEN];
	PetscBool      flg;
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	
	flg = PETSC_FALSE;
	ierr = PetscOptionsGetString(PETSC_NULL,"-output_path",outputpath,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (flg == PETSC_FALSE) { 
		sprintf(outputpath,".");
	}
	
	flg = PETSC_FALSE;
	ierr = PetscOptionsGetString(PETSC_NULL,"-file_prefix",prefix,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (flg == PETSC_FALSE) { 
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Required to know the prefix of the checkpointed files");
	}
	
	sprintf(vname,"%s/%s.dmda-velocity",outputpath,prefix);
	sprintf(pname,"%s/%s.dmda-pressure",outputpath,prefix);
	sprintf(xuname,"%s/%s.dmda-Xu",outputpath,prefix);
	sprintf(xpname,"%s/%s.dmda-Xp",outputpath,prefix);
	
	flg = PETSC_FALSE;
	ierr = PetscOptionsGetString(PETSC_NULL,"-file_suffix",suffix,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (flg) { 
		//sprintf(vname,"%s%s",vname,suffix);
		//sprintf(pname,"%s%s",pname,suffix);
		//sprintf(xuname,"%s%s",xuname,suffix);
		//sprintf(xpname,"%s%s",xpname,suffix);
		_strlcat(vname,suffix,PETSC_MAX_PATH_LEN-1);
		_strlcat(pname,suffix,PETSC_MAX_PATH_LEN-1);
		_strlcat(xuname,suffix,PETSC_MAX_PATH_LEN-1);
		_strlcat(xpname,suffix,PETSC_MAX_PATH_LEN-1);
	}
	PetscPrintf(PETSC_COMM_WORLD,"Reading files: (dmda-vel) %s\n\t\t(dmda-p) %s\n\t\t(vec-u : p) %s : %s\n", vname,pname,xuname,xpname);
	
	ierr = PhysCompStokesLoad_DM(vname,pname,&ctx);CHKERRQ(ierr);
	ierr = PhysCompStokesLoad_VP(ctx,xuname,xpname,PETSC_TRUE,&X);CHKERRQ(ierr);
	
	sprintf(outfilename,"%s",prefix);
	if (flg) { 
		//sprintf(outfilename,"%s%s",outfilename,suffix);
		_strlcat(outfilename,suffix,PETSC_MAX_PATH_LEN-1);
	}
	//sprintf(outfilename,"%s-vp",outfilename);
	_strlcat(outfilename,"-vp",PETSC_MAX_PATH_LEN-1);
	PetscPrintf(PETSC_COMM_WORLD,"Writing file: %s/%s \n", outputpath,outfilename);
	ierr = PhysCompStokesWrite_DM_X(ctx,X,outputpath,outfilename);CHKERRQ(ierr);
	
	//ierr = PhysCompOutput_StokesRawVelocityPressure(ctx,X,PETSC_TRUE,outputpath,suffix);CHKERRQ(ierr);
	
	ierr = PhysCompDestroy_Stokes(&ctx);CHKERRQ(ierr);
	ierr = VecDestroy(&X);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}



int main( int argc,char **argv )
{
	PetscErrorCode ierr;
	PetscLogDouble t0,t1,gt,lt;
	
	PetscInitialize(&argc,&argv,(char *)0,0);
	
//	ierr = test_DMDACheckPoint();CHKERRQ(ierr);

	PetscGetTime(&t0);
	ierr = test_LoadStokesFromCheckpoint_WriteToVTS();CHKERRQ(ierr);
	PetscGetTime(&t1);
	lt = t1-t0;
	ierr = MPI_Allreduce(&lt,&gt,1,MPIU_REAL,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"test_LoadStokesFromCheckpoint_WriteToVTS: CPU time %1.4e (sec)\n", gt);
	
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}
