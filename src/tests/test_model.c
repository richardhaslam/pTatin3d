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
 **    filename:   test_model.c
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
/*

Run with:

>./test_model.app -model_template_param1 11.22 -model_template_param2 6699

*/


#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"
#include "ptatin_init.h"

extern PetscErrorCode pTatinModelRegister_Template(void);

/* define user model */
typedef struct {
	double alpha;
	double beta;
} UserModelCtx;


PetscErrorCode ModelTest_ApplyBoundaryCondition(pTatinCtx c,void *ctx)
{
	UserModelCtx *data = (UserModelCtx*)ctx;
	printf("alpha = %lf \n", data->alpha );
	PetscFunctionReturn(0);
}
PetscErrorCode ModelTest_ApplyInitialMaterialGeometry(pTatinCtx c,void *ctx)
{
	UserModelCtx *data = (UserModelCtx*)ctx;
	printf("beta = %lf \n", data->beta );
	PetscFunctionReturn(0);
}


PetscErrorCode test_model(void)
{
	UserModelCtx *data;
	pTatinModel m,model;
	PetscErrorCode ierr;

	data = malloc(sizeof(UserModelCtx));
	data->alpha = 11.23;
	data->beta = 66.99;

	/* register user model */
	ierr = pTatinModelCreate(&m);CHKERRQ(ierr);
	ierr = pTatinModelSetName(m,"test_model");CHKERRQ(ierr);

	ierr = pTatinModelSetUserData(m,data);CHKERRQ(ierr);

	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_BC,(void (*)(void))ModelTest_ApplyBoundaryCondition);CHKERRQ(ierr);
	ierr = pTatinModelSetFunctionPointer(m,PTATIN_MODEL_APPLY_INIT_MAT_GEOM,(void (*)(void))ModelTest_ApplyInitialMaterialGeometry);CHKERRQ(ierr);

	ierr = pTatinModelRegister(m);CHKERRQ(ierr);
//	ierr = pTatinModelRegister(m);CHKERRQ(ierr); /* test for duplicate model name insertition */

	ierr = pTatinModelGetByName("test_model",&model);CHKERRQ(ierr);
	printf("m     = %p \n", m );
	printf("model = %p \n", model );
	printf("data    = %p \n", data );
	printf("m->data = %p \n", m->model_data );


	{
		void *userdata;
		ierr = pTatinModelRegister_Template();CHKERRQ(ierr);
		ierr = pTatinModelGetByName("template",&model);CHKERRQ(ierr);

		ierr = pTatinModelGetByName("template",&model);CHKERRQ(ierr);
		ierr = pTatinModelGetUserData(model,&userdata);CHKERRQ(ierr);

		ierr = model->FP_pTatinModel_Initialize(0,userdata);CHKERRQ(ierr);
		ierr = model->FP_pTatinModel_ApplyBoundaryCondition(0,userdata);CHKERRQ(ierr);
		ierr = model->FP_pTatinModel_ApplyInitialMeshGeometry(0,userdata);CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}

int main(int argc,char **argv)
{
	PetscErrorCode ierr;

	ierr = pTatinInitialize(&argc,&argv,(char *)0,NULL);CHKERRQ(ierr);

	test_model();

	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}
