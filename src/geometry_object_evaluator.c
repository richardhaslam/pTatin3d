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
 **    filename:   geometry_object_evaluator.c
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


#include "stdio.h"
#include "stdlib.h"
#include "ctype.h"
#include "math.h"
#include "petsc.h"
#include "geometry_object.h"
#include "geometry_object_evaluator.h"


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
	*G = NULL;

	PetscFunctionReturn(0);
}

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

PetscErrorCode GeometryObjectEvalGetGeometryObject(GeometryObjectEval goe,GeometryObject *go)
{
	if (go) {
		*go = goe->go;
	}
	PetscFunctionReturn(0);
}

PetscErrorCode GeometryObjectEvalSetRegionIndex(GeometryObjectEval goe,int region_id)
{
	goe->region_index = region_id;
	PetscFunctionReturn(0);
}

PetscErrorCode GeometryObjectEvalSetRegionValue(GeometryObjectEval goe,double value)
{
	goe->region_value = value;
	PetscFunctionReturn(0);
}

PetscErrorCode GeometryObjectEvalSetRegionFunction(GeometryObjectEval goe,double (*func)(double*,void*),void *data)
{
	goe->data = data;
	goe->evaluate_region_function = func;
	PetscFunctionReturn(0);
}

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
