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
**    filename:   model_ops_template.c
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
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin3d_stokes.h"
#include "ptatin3d_energy.h"
#include "ptatin_models.h"
#include "model_template_ctx.h"


PetscErrorCode ModelInitialize_Template(pTatinCtx c,void *ctx)
{
  ModelTemplateCtx *data;
  PetscBool        flg;
  PetscErrorCode   ierr;

  PetscFunctionBegin;

  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);
  data = (ModelTemplateCtx*)ctx;
  ierr = PetscOptionsGetReal(NULL,NULL,"-model_template_param1",&data->param1,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL, "-model_template_param2",&data->param2,&flg);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMeshGeometry_Template(pTatinCtx c,void *ctx)
{
  /* ModelTemplateCtx *data; */
  PhysCompStokes   stokes;
  DM               stokes_pack,dav,dap;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);


  /* data = (ModelTemplateCtx*)ctx; */
  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
  stokes_pack = stokes->stokes_pack;
  ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(dav,0.0,1.0,0.0,1.0,0.0,0.1);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialMaterialGeometry_Template(pTatinCtx c,void *ctx)
{
  /* ModelTemplateCtx *data; */

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);
  /* data = (ModelTemplateCtx*)ctx; */

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyInitialSolution_Template(pTatinCtx c,Vec X,void *ctx)
{
  /* ModelTemplateCtx *data; */

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);
  /* data = (ModelTemplateCtx*)ctx; */

  PetscFunctionReturn(0);
}

PetscErrorCode Template_VelocityBC(BCList bclist,DM dav,pTatinCtx c,ModelTemplateCtx *data)
{
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);
  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyBoundaryCondition_Template(pTatinCtx c,void *ctx)
{
  ModelTemplateCtx *data;
  PhysCompStokes   stokes;
  DM               stokes_pack,dav,dap;
  PetscErrorCode   ierr;

  PetscFunctionBegin;

  data = (ModelTemplateCtx*)ctx;

  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);
  PetscPrintf(PETSC_COMM_WORLD,"param1 = %lf \n", data->param1 );

  /* Define velocity boundary conditions */
  ierr = pTatinGetStokesContext(c,&stokes);CHKERRQ(ierr);
  stokes_pack = stokes->stokes_pack;
  ierr = DMCompositeGetEntries(stokes_pack,&dav,&dap);CHKERRQ(ierr);

  ierr = Template_VelocityBC(stokes->u_bclist,dav,c,data);CHKERRQ(ierr);

  /* Define boundary conditions for any other physics */

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyBoundaryConditionMG_Template(PetscInt nl,BCList bclist[],DM dav[],pTatinCtx c,void *ctx)
{
  ModelTemplateCtx *data;
  PetscInt         n;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);

  data = (ModelTemplateCtx*)ctx;
  /* Define velocity boundary conditions on each level within the MG hierarchy */
  for (n=0; n<nl; n++) {
    ierr = Template_VelocityBC(bclist[n],dav[n],c,data);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyMaterialBoundaryCondition_Template(pTatinCtx c,void *ctx)
{
  /* ModelTemplateCtx *data; */

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);
  /* data = (ModelTemplateCtx*)ctx; */

  PetscFunctionReturn(0);
}

PetscErrorCode ModelApplyUpdateMeshGeometry_Template(pTatinCtx c,Vec X,void *ctx)
{
  /* ModelTemplateCtx *data; */

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);
  /* data = (ModelTemplateCtx*)ctx; */

  PetscFunctionReturn(0);
}

PetscErrorCode ModelOutput_Template(pTatinCtx c,Vec X,const char prefix[],void *ctx)
{
  /* ModelTemplateCtx *data; */

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);
  /* data = (ModelTemplateCtx*)ctx; */

  /* ---- Velocity-Pressure Mesh Output ---- */
  /* [1] Standard viewer: v,p written out as binary in double */
  /*
  ierr = pTatin3d_ModelOutput_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
  */
  /* [2] Light weight viewer: Only v is written out. v and coords are expressed as floats */
  /*
  ierr = pTatin3d_ModelOutputLite_Velocity_Stokes(c,X,prefix);CHKERRQ(ierr);
  */
  /* [3] Write out v,p into PETSc Vec. These can be used to restart pTatin */
  /*
  ierr = pTatin3d_ModelOutputPetscVec_VelocityPressure_Stokes(c,X,prefix);CHKERRQ(ierr);
  */


  /* ---- Material Point Output ---- */
  /* [1] Basic viewer: Only reports coords, regionid and other internal data */
  /*
  ierr = pTatin3d_ModelOutput_MPntStd(c,prefix);CHKERRQ(ierr);
  */

  /* [2] Customized viewer: User defines specific fields they want to view - NOTE not .pvd file will be created */
  /*
  {
  DataBucket                materialpoint_db;
  const int                 nf = 4;
  const MaterialPointField  mp_prop_list[] = { MPField_Std, MPField_Stokes, MPField_StokesPl, MPField_Energy };
  char                      mp_file_prefix[256];

  ierr = pTatinGetMaterialPoints(c,&materialpoint_db,NULL);CHKERRQ(ierr);
  sprintf(mp_file_prefix,"%s_mpoints",prefix);
  ierr = SwarmViewGeneric_ParaView(materialpoint_db,nf,mp_prop_list,c->outputpath,mp_file_prefix);CHKERRQ(ierr);
  }
  */
  /* [3] Customized marker->cell viewer: Marker data is projected onto the velocity mesh. User defines specific fields */
  /*
  {
  const int                    nf = 3;
  const MaterialPointVariable  mp_prop_list[] = { MPV_viscosity, MPV_density, MPV_plastic_strain };

  ierr = pTatin3d_ModelOutput_MarkerCellFields(c,nf,mp_prop_list,prefix);CHKERRQ(ierr);
  }
  */

  PetscFunctionReturn(0);
}

PetscErrorCode ModelDestroy_Template(pTatinCtx c,void *ctx)
{
  ModelTemplateCtx *data;
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"[[%s]]\n", PETSC_FUNCTION_NAME);
  data = (ModelTemplateCtx*)ctx;

  /* Free contents of structure */

  /* Free structure */
  ierr = PetscFree(data);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode pTatinModelRegister_Template(void)
{
  ModelTemplateCtx *data;
  pTatinModel      m;
  PetscErrorCode   ierr;

  PetscFunctionBegin;

  /* Allocate memory for the data structure for this model */
  ierr = PetscMalloc(sizeof(ModelTemplateCtx),&data);CHKERRQ(ierr);
  ierr = PetscMemzero(data,sizeof(ModelTemplateCtx));CHKERRQ(ierr);

  /* set initial values for model parameters */
  data->param1 = 0.0;
  data->param2 = 0;

  /* register user model */
  ierr = pTatinModelCreate(&m);CHKERRQ(ierr);

  /* Set name, model select via -ptatin_model NAME */
  ierr = pTatinModelSetName(m,"template");CHKERRQ(ierr);

  /* Set model data */
  ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);

  /* Set function pointers */
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_INIT,                  (void (*)(void))ModelInitialize_Template);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MESH_GEOM,  (void (*)(void))ModelApplyInitialMeshGeometry_Template);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,   (void (*)(void))ModelApplyInitialMaterialGeometry_Template);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_SOLUTION,   (void (*)(void))ModelApplyInitialSolution_Template);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,              (void (*)(void))ModelApplyBoundaryCondition_Template);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BCMG,            (void (*)(void))ModelApplyBoundaryConditionMG_Template);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_MAT_BC,          (void (*)(void))ModelApplyMaterialBoundaryCondition_Template);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_UPDATE_MESH_GEOM,(void (*)(void))ModelApplyUpdateMeshGeometry_Template);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_OUTPUT,                (void (*)(void))ModelOutput_Template);CHKERRQ(ierr);
  ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_DESTROY,               (void (*)(void))ModelDestroy_Template);CHKERRQ(ierr);

  /* Insert model into list */
  ierr = pTatinModelRegister(m);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
