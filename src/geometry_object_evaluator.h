
#ifndef __geometry_object_eval_h__
#define __geometry_object_eval_h__


#define GEOM_OBJECT_EVAL_REGION_VALUE_NULL -1.0e32
#define GEOM_OBJECT_EVAL_REGION_INDEX_NULL -69


typedef struct _p_GeometryObjectEval *GeometryObjectEval;

struct _p_GeometryObjectEval {
  /* user methods */
	GeometryObject go;
	void    *data;
	int     region_index;
	double  region_value;
	double  (*evaluate_region_function)(double*,void*);
};


PetscErrorCode GeometryObjectEvalCreate(const char name[],GeometryObjectEval *G);
PetscErrorCode GeometryObjectEvalDestroy(GeometryObjectEval *G);

PetscErrorCode GeometryObjectEvalSetGeometryObject(GeometryObjectEval goe,GeometryObject go);
PetscErrorCode GeometryObjectEvalGetGeometryObject(GeometryObjectEval goe,GeometryObject *go);

PetscErrorCode GeometryObjectEvalSetRegionIndex(GeometryObjectEval goe,int region_id);
PetscErrorCode GeometryObjectEvalSetRegionValue(GeometryObjectEval goe,double value);
PetscErrorCode GeometryObjectEvalSetRegionFunction(GeometryObjectEval goe,double (*func)(double*,void*),void *data);

PetscErrorCode GeometryObjectEvaluateRegionIndex(GeometryObjectEval goe,double pos[],int *region,PetscBool *assigned);
PetscErrorCode GeometryObjectEvaluateRegionValue(GeometryObjectEval goe,double pos[],double *value,PetscBool *assigned);
PetscErrorCode GeometryObjectEvaluateRegionFunction(GeometryObjectEval goe,double pos[],double *value,PetscBool *assigned);

PetscErrorCode GeometryObjectEvalFindByName(GeometryObjectEval G[],const char name[],GeometryObjectEval *g);
PetscErrorCode GeometryObjectEvalIdFindByName(GeometryObjectEval G[],const char name[],PetscInt *GoId);

/* macro wrappers for geometry object */
#define GeometryObjectEvalRotate(goe,dir,angle) \
{ PetscErrorCode ierr = GeometryObjectRotate( (goe)->go,(dir),(axis) );CHKERRQ(ierr); }

#define GeometryObjectEvalPointInside(goe,pos,inside) \
{ PetscErrorCode ierr = GeometryObjectPointInside( (goe)->go,(pos),(inside) );CHKERRQ(ierr); }


#endif