/*
 
Run with:

>./test_model.app -model_template_param1 11.22 -model_template_param2 6699
 
*/


#define _GNU_SOURCE
#include "petsc.h"
#include "ptatin3d.h"
#include "private/ptatin_impl.h"

extern PetscErrorCode pTatinModelRegister_Template(void);

/* define user model */
typedef struct {
	double alpha;
	double beta;
} UserModelCtx;


#undef __FUNCT__
#define __FUNCT__ "ModelTest_ApplyBoundaryCondition"
PetscErrorCode ModelTest_ApplyBoundaryCondition(pTatinCtx c,void *ctx)
{
	UserModelCtx *data = (UserModelCtx*)ctx;
	printf("alpha = %lf \n", data->alpha );
	PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "ModelTest_ApplyInitialMaterialGeometry"
PetscErrorCode ModelTest_ApplyInitialMaterialGeometry(pTatinCtx c,void *ctx)
{
	UserModelCtx *data = (UserModelCtx*)ctx;
	printf("beta = %lf \n", data->beta );
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "test_model"
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
	
	PetscInitialize(&argc,&argv,(char *)0,0);

	test_model();
	
	ierr = PetscFinalize();CHKERRQ(ierr);
	return 0;
}
