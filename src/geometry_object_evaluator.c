

#include "stdio.h"
#include "stdlib.h"
#include "ctype.h"
#include "math.h"
#include "petsc.h"
#include "geometry_object.h"
#include "geometry_object_evaluator.h"


#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvalCreate"
PetscErrorCode GeometryObjectEvalCreate(const char name[],GeometryObjectEval *G)
{
	GeometryObjectEval goe;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeometryObjectEval),&goe);CHKERRQ(ierr);
	ierr = PetscMemzero(goe,sizeof(struct _p_GeometryObjectEval));CHKERRQ(ierr);
	
	ierr = GeometryObjectCreate(name,&goe->go);CHKERRQ(ierr);
	
	*G = goe;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvalDestroy"
PetscErrorCode GeometryObjectEvalDestroy(GeometryObjectEval *G)
{
	GeometryObjectEval goe;
	PetscErrorCode ierr;
	
	
	if (!G) { PetscFunctionReturn(0); }
	goe = *G;
	
	if (goe->go) {
		ierr = GeometryObjectDestroy(&goe->go);CHKERRQ(ierr);
	}
	
	ierr = PetscFree(goe);CHKERRQ(ierr);	
	*G = PETSC_NULL;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvalSetGeometryObject"
PetscErrorCode GeometryObjectEvalSetGeometryObject(GeometryObjectEval goe,GeometryObject go)
{
	PetscErrorCode ierr;
	
	if (goe->go) {
		ierr = GeometryObjectDestroy(&goe->go);CHKERRQ(ierr);
	}
	/* increment reference counter as someone less might be using this GeometryObject */
	go->ref_cnt++;
	
	goe->go = go;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvalGetGeometryObject"
PetscErrorCode GeometryObjectEvalGetGeometryObject(GeometryObjectEval goe,GeometryObject *go)
{
	if (go) {
		*go = goe->go;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvalSetRegionIndex"
PetscErrorCode GeometryObjectEvalSetRegionIndex(GeometryObjectEval goe,int region_id)
{
	goe->region_index = region_id;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvalSetRegionValue"
PetscErrorCode GeometryObjectEvalSetRegionValue(GeometryObjectEval goe,double value)
{
	goe->region_value = value;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvalSetRegionFunction"
PetscErrorCode GeometryObjectEvalSetRegionFunction(GeometryObjectEval goe,double (*func)(double*,void*),void *data)
{
	goe->data = data;
	goe->evaluate_region_function = func;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvaluateRegionIndex"
PetscErrorCode GeometryObjectEvaluateRegionIndex(GeometryObjectEval goe,double pos[],int *region,PetscBool *assigned)
{
	int inside;
	PetscErrorCode ierr;
	
	inside = 0;
	ierr = GeometryObjectPointInside(goe->go,pos,&inside);CHKERRQ(ierr);
	
	if (assigned) {
			switch (inside) {
				case 0:
					*assigned = PETSC_FALSE;
					break;
				case 1:
					*assigned = PETSC_TRUE;
					break;
			}
	}
	
	*region = GEOM_OBJECT_EVAL_REGION_INDEX_NULL;
	if (inside == 1) {
		*region = goe->region_index;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvaluateRegionValue"
PetscErrorCode GeometryObjectEvaluateRegionValue(GeometryObjectEval goe,double pos[],double *value,PetscBool *assigned)
{
	int inside;
	PetscErrorCode ierr;
	
	
	inside = 0;
	ierr = GeometryObjectPointInside(goe->go,pos,&inside);CHKERRQ(ierr);
	
	if (assigned) {
		switch (inside) {
			case 0:
				*assigned = PETSC_FALSE;
				break;
			case 1:
				*assigned = PETSC_TRUE;
				break;
		}
	}
	
	*value = GEOM_OBJECT_EVAL_REGION_VALUE_NULL;
	if (inside == 1) {
		*value = goe->region_value;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvaluateRegionFunction"
PetscErrorCode GeometryObjectEvaluateRegionFunction(GeometryObjectEval goe,double pos[],double *value,PetscBool *assigned)
{
	int inside;
	PetscErrorCode ierr;
	
	
	inside = 0;
	ierr = GeometryObjectPointInside(goe->go,pos,&inside);CHKERRQ(ierr);
	
	if (assigned) {
		switch (inside) {
			case 0:
				*assigned = PETSC_FALSE;
				break;
			case 1:
				*assigned = PETSC_TRUE;
				break;
		}
	}
	
	*value = GEOM_OBJECT_EVAL_REGION_VALUE_NULL;
	if (inside == 1) {
		if (goe->evaluate_region_function) {
			*value = goe->evaluate_region_function(pos,goe->data);
		} else {
			SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"No method for go->evaluate_region_function provided");
		}
	}
	PetscFunctionReturn(0);
}

/* wrappers */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectFindByName"
PetscErrorCode GeometryObjectEvalFindByName(GeometryObjectEval G[],const char name[],GeometryObjectEval *g)
{
	GeometryObjectEval item;
	int i,v;
	
	*g = NULL;
	
	i = 0;
	item = G[i];
	while (item != NULL) {
		v = strcmp(name,G[i]->go->name);
		if (v == 0) {
			*g = item;
			break;
		}
		i++;
		item = G[i];
	}
	if (*g == NULL) {
		PetscPrintf(PETSC_COMM_SELF,"[Warning] GeomObjectEval with GeomObject name %s was not found in list\n",name);
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvalIdFindByName"
PetscErrorCode GeometryObjectEvalIdFindByName(GeometryObjectEval G[],const char name[],PetscInt *GoId)
{
	GeometryObjectEval item;
	int i,v;
	
	*GoId = -1;
	
	i = 0;
	item = G[i];
	while (item != NULL) {
		v = strcmp(name,G[i]->go->name);
		if (v == 0) {
			*GoId = i;
			break;
		}
		i++;
		item = G[i];
	}
	if (*GoId == -1) {
		PetscPrintf(PETSC_COMM_SELF,"[Warning] GeomObjectEval with GeomObject name %s was not found in list\n",name);
	}
	PetscFunctionReturn(0);
}