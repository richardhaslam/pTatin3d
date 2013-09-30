
#include "stdio.h"
#include "stdlib.h"
#include "ctype.h"
#include "math.h"
#include "petsc.h"
#include "geometry_object.h"


const char *GeomTypeNames[] = {
	"GeomType_Box", 
	"GeomType_Cylinder", 
	"GeomType_Sphere",
	"GeomType_SetOperation",
	0
};

const char *GeomTypeSetOperatorNames[] = {
	"GeomType_SetUnion", 
	"GeomType_SetIntersection", 
	"GeomType_SetComplement",
	0
};

const char *GeomRotateAxisNames[] = {
	"ROTATE_AXIS_X", 
	"ROTATE_AXIS_Y", 
	"ROTATE_AXIS_Z",
	0
};

void PointRotate(double xin[],GeomRotateAxis axis,double angle,double xout[]);
void PointBackRotate(double xin[],GeomRotateAxis axis,double angle,double xout[]);

/* ------------------------------------------------------------------------------------------- */
/* API */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectCreate"
PetscErrorCode GeometryObjectCreate(const char name[],int region_index,double value,double (*fp)(double*),GeometryObject *G)
{
	GeometryObject go;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeometryObject),&go);CHKERRQ(ierr);
	ierr = PetscMemzero(go,sizeof(struct _p_GeometryObject));CHKERRQ(ierr);
	
	asprintf(&go->name,"%s",name);
	go->type = GeomType_NULL;
	go->ctx = PETSC_NULL;
	go->n_rotations = 0;
	go->ref_cnt = 0;
	
	go->region_index = region_index;
	go->value = value;
	go->evaluate_region_function = fp;
	
	go->geom_point_inside = PETSC_NULL;
	go->geom_transform_translate = PETSC_NULL;
	
	*G = go;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectDestroy"
PetscErrorCode GeometryObjectDestroy(GeometryObject *G)
{
	GeometryObject go;
	PetscErrorCode ierr;
	
	
	if (!G) { PetscFunctionReturn(0); }
	go = *G;
	if (go->ref_cnt > 0) {
		go->ref_cnt--;
		PetscFunctionReturn(0);
	}
	
	if (go->geom_destroy) {
		ierr = go->geom_destroy(go);CHKERRQ(ierr);
	}
	free(go->name);
	ierr = PetscFree(go);CHKERRQ(ierr);
	
	*G = PETSC_NULL;
	
	PetscFunctionReturn(0);
}

/* rotate */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectRotate"
PetscErrorCode GeometryObjectRotate(GeometryObject go,GeomRotateAxis dir,double angle)
{
	double angle2;
	
	
	if (go->n_rotations == GEOM_SHAPE_MAX_ROTATIONS-1) {
		printf("No more rotations permitted\n");
		exit(0);
	}
	
	angle2 = angle * M_PI / 180.0;
	
	go->rotation_axis[  go->n_rotations ] = dir;
	go->rotation_angle[ go->n_rotations ] = angle2;
	go->n_rotations++;
	
	PetscFunctionReturn(0);
}

/* translate */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectTransformTranslate"
PetscErrorCode GeometryObjectTransformTranslate(GeometryObject go,double shift[])
{
	/* call method */
	if (go->geom_transform_translate) {
		go->geom_transform_translate(go,shift);
	} else {
		printf("No method for go->geom_transform_translate provided\n");
		exit(0);
	}
	PetscFunctionReturn(0);
}

/* point inside */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPointInside"
PetscErrorCode GeometryObjectPointInside(GeometryObject go,double pos[],int *inside)
{
	int    nr;
	double posin[3],posout[3];
	
	
	posin[0] = pos[0];
	posin[1] = pos[1];
	posin[2] = pos[2];
	
	for (nr=go->n_rotations-1; nr>=0; nr--) {
		PointBackRotate(posin,go->rotation_axis[nr],go->rotation_angle[nr],posout);
		
		posin[0] = posout[0];
		posin[1] = posout[1];
		posin[2] = posout[2];
	}
	
	/* call method */
	if (go->geom_point_inside) {
		*inside = 0;
		go->geom_point_inside(go,posin,inside);
	} else {
		printf("No method for go->geom_point_inside provided\n");
		exit(0);
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvaluateRegionIndex"
PetscErrorCode GeometryObjectEvaluateRegionIndex(GeometryObject go,double pos[],int *region)
{
	int inside;
	PetscErrorCode ierr;
	
	
	inside = 0;
	ierr = GeometryObjectPointInside(go,pos,&inside);CHKERRQ(ierr);
	
	*region = GEOM_SHAPE_REGION_INDEX_NULL;
	if (inside == 1) {
		*region = go->region_index;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvaluateRegionValue"
PetscErrorCode GeometryObjectEvaluateRegionValue(GeometryObject go,double pos[],double *value)
{
	int inside;
	PetscErrorCode ierr;
	
	
	inside = 0;
	ierr = GeometryObjectPointInside(go,pos,&inside);CHKERRQ(ierr);
	
	*value = GEOM_SHAPE_REGION_VALUE_NULL;
	if (inside == 1) {
		*value = go->value;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectEvaluateRegionFunction"
PetscErrorCode GeometryObjectEvaluateRegionFunction(GeometryObject go,double pos[],double *value)
{
	int inside;
	PetscErrorCode ierr;
	
	
	inside = 0;
	ierr = GeometryObjectPointInside(go,pos,&inside);CHKERRQ(ierr);
	
	*value = GEOM_SHAPE_REGION_VALUE_NULL;
	if (inside == 1) {
		if (go->evaluate_region_function) {
			*value = go->evaluate_region_function(pos);
		} else {
			printf("No method for go->evaluate_region_function provided\n");
			exit(0);
		}
	}
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------------------------- */
/* helpers */
#undef __FUNCT__
#define __FUNCT__ "GeomTypeNameGetId"
PetscErrorCode GeomTypeNameGetId(const char name[],int *id)
{
	const char *item;
	int i,v;
	
	i = 0;
	item = GeomTypeNames[i];
	
	*id = -1;
	while ( item != 0 ) {
		v = strcmp(item,name);
		if (v == 0) {
			*id = i;
			PetscFunctionReturn(0);
		}
		i++;
		item = GeomTypeNames[i];
	}
	SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"GeomType %s not located in list of available GeomTypes\n",name);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectFindByName"
PetscErrorCode GeometryObjectFindByName(GeometryObject G[],const char name[],GeometryObject *g)
{
	GeometryObject item;
	int i,v;
	
	*g = NULL;
	
	i = 0;
	item = G[i];
	while (item != NULL) {
		v = strcmp(name,G[i]->name);
		if (v == 0) {
			*g = item;
			break;
		}
		i++;
		item = G[i];
	}
	if (*g == NULL) {
		printf("[Warning] GeomObject with name %s was not found in list\n",name);
	}
	PetscFunctionReturn(0);
}

void PointRotate(double xin[],GeomRotateAxis axis,double angle,double xout[])
{
	double RR[3][3],xin2[3];
	double c_ang,s_ang;
	
	c_ang = cos(angle);
	s_ang = sin(angle);
	
	switch (axis) {
			
		case ROTATE_AXIS_X:
		{
			RR[0][0] = 1.0;		RR[0][1] = 0.0;			RR[0][2] = 0.0;
			RR[1][0] = 0.0;		RR[1][1] = c_ang;		RR[1][2] = -s_ang;
			RR[2][0] = 0.0;		RR[2][1] = s_ang;		RR[2][2] =  c_ang;
		}
			break;
			
		case ROTATE_AXIS_Y:
		{
			RR[0][0] =  c_ang;		RR[0][1] = 0.0;		RR[0][2] = s_ang;
			RR[1][0] = 0.0;				RR[1][1] = 1.0;		RR[1][2] = 0.0;
			RR[2][0] = -s_ang;		RR[2][1] = 0.0;		RR[2][2] = c_ang;
		}
			break;
			
		case ROTATE_AXIS_Z:
		{
			RR[0][0] = c_ang;		RR[0][1] = -s_ang;	RR[0][2] = 0.0;
			RR[1][0] = s_ang;		RR[1][1] =  c_ang;	RR[1][2] = 0.0;
			RR[2][0] = 0.0;			RR[2][1] = 0.0;			RR[2][2] = 1.0;
		}
			break;
	}
	
	xin2[0] = xin[0];
	xin2[1] = xin[1];
	xin2[2] = xin[2];
	
	xout[0] = RR[0][0]*xin2[0] + RR[0][1]*xin2[1] + RR[0][2]*xin2[2];
	xout[1] = RR[1][0]*xin2[0] + RR[1][1]*xin2[1] + RR[1][2]*xin2[2];
	xout[2] = RR[2][0]*xin2[0] + RR[2][1]*xin2[1] + RR[2][2]*xin2[2];
}

void PointBackRotate(double xin[],GeomRotateAxis axis,double angle,double xout[])
{
	double angle2;
	
	angle2 = -angle;
	PointRotate(xin,axis,angle2,xout);
}

/* ------------------------------------------------------------------------------------------- */
/* Implementations */

/* BOX: implementation */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectGetContext_Box"
PetscErrorCode GeometryObjectGetContext_Box(GeometryObject go,GeomTypeBox *b)
{
	*b = (GeomTypeBox)(go->ctx);	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectDestroy_Box"
PetscErrorCode GeometryObjectDestroy_Box(GeometryObject go)
{
	GeomTypeBox box;
	PetscErrorCode ierr;
	
	
	ierr = GeometryObjectGetContext_Box(go,&box);CHKERRQ(ierr);
	ierr = PetscFree(box);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetFromOptions_Box"
PetscErrorCode GeometryObjectSetFromOptions_Box(GeometryObject go)
{
	GeomTypeBox box;
	PetscReal val;
	PetscBool flg;
	char name[1024];
	PetscErrorCode ierr;
	
	
	ierr = GeometryObjectGetContext_Box(go,&box);CHKERRQ(ierr);
	
	sprintf(name,"%s_",go->name);
	ierr = PetscOptionsGetReal(name,"-x0",&val,&flg);CHKERRQ(ierr);
	if (flg) { box->x0[0] = val; }
	ierr = PetscOptionsGetReal(name,"-x1",&val,&flg);CHKERRQ(ierr);
	if (flg) { box->x0[1] = val; }
	ierr = PetscOptionsGetReal(name,"-x2",&val,&flg);CHKERRQ(ierr);
	if (flg) { box->x0[2] = val; }
	
	ierr = PetscOptionsGetReal(name,"-Lx0",&val,&flg);CHKERRQ(ierr);
	if (flg) { box->Lx[0] = val; }
	ierr = PetscOptionsGetReal(name,"-Lx1",&val,&flg);CHKERRQ(ierr);
	if (flg) { box->Lx[1] = val; }
	ierr = PetscOptionsGetReal(name,"-Lx2",&val,&flg);CHKERRQ(ierr);
	if (flg) { box->Lx[2] = val; }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectTransformTranslate_Box"
PetscErrorCode GeometryObjectTransformTranslate_Box(GeometryObject go,double shift[])
{
	GeomTypeBox box;
	PetscErrorCode ierr;
	
	
	ierr = GeometryObjectGetContext_Box(go,&box);CHKERRQ(ierr);
	
	box->x0[0] = box->x0[0] + shift[0];
	box->x0[1] = box->x0[1] + shift[1];
	box->x0[2] = box->x0[2] + shift[2];
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPointInside_Box"
PetscErrorCode GeometryObjectPointInside_Box(GeometryObject go,double pos[],int *inside)
{
	GeomTypeBox box;
	PetscErrorCode ierr;
	
	
	ierr = GeometryObjectGetContext_Box(go,&box);CHKERRQ(ierr);
	
	if (pos[0] < box->x0[0])              { PetscFunctionReturn(0); }
	if (pos[0] > box->x0[0] + box->Lx[0]) { PetscFunctionReturn(0); }
	
	if (pos[1] < box->x0[1])              { PetscFunctionReturn(0); }
	if (pos[1] > box->x0[1] + box->Lx[1]) { PetscFunctionReturn(0); }
	
	if (pos[2] < box->x0[2])              { PetscFunctionReturn(0); }
	if (pos[2] > box->x0[2] + box->Lx[2]) { PetscFunctionReturn(0); }
	
	*inside = 1;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetType_Box"
PetscErrorCode GeometryObjectSetType_Box(GeometryObject go,double x0[],double Lx[])
{
	GeomTypeBox box;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeomTypeBox),&box);CHKERRQ(ierr);
	ierr = PetscMemzero(box,sizeof(struct _p_GeomTypeBox));CHKERRQ(ierr);
	box->x0[0] = x0[0];
	box->x0[1] = x0[1];
	box->x0[2] = x0[2];
	box->Lx[0] = Lx[0];
	box->Lx[1] = Lx[1];
	box->Lx[2] = Lx[2];
	
	go->type                     = GeomType_Box;
	go->geom_point_inside        = GeometryObjectPointInside_Box;
	go->geom_transform_translate = GeometryObjectTransformTranslate_Box;
	go->geom_destroy             = GeometryObjectDestroy_Box;
	go->ctx = (void*)box;
	
	ierr = GeometryObjectSetFromOptions_Box(go);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------------------------- */
/* SPHERE: implemenetation */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectDestroy_Sphere"
PetscErrorCode GeometryObjectDestroy_Sphere(GeometryObject go)
{
	GeomTypeSphere ctx;
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeSphere)go->ctx;
	ierr = PetscFree(ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetFromOptions_Sphere"
PetscErrorCode GeometryObjectSetFromOptions_Sphere(GeometryObject go)
{
	GeomTypeSphere ctx;
	PetscReal val;
	PetscBool flg;
	char name[1024];
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeSphere)go->ctx;
	
	sprintf(name,"%s_",go->name);
	ierr = PetscOptionsGetReal(name,"-Ox",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->origin[0] = val; }
	ierr = PetscOptionsGetReal(name,"-Oy",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->origin[1] = val; }
	ierr = PetscOptionsGetReal(name,"-Oz",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->origin[2] = val; }
	
	ierr = PetscOptionsGetReal(name,"-rad",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->radius = val; }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectTransformTranslate_Sphere"
PetscErrorCode GeometryObjectTransformTranslate_Sphere(GeometryObject go,double shift[])
{
	GeomTypeSphere ctx;
	
	
	ctx = (GeomTypeSphere)go->ctx;
	
	ctx->origin[0] = ctx->origin[0] + shift[0];
	ctx->origin[1] = ctx->origin[1] + shift[1];
	ctx->origin[2] = ctx->origin[2] + shift[2];
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPointInside_Sphere"
PetscErrorCode GeometryObjectPointInside_Sphere(GeometryObject go,double pos[],int *inside)
{
	GeomTypeSphere ctx;
	double r2,sep2;
	
	
	ctx = (GeomTypeSphere)go->ctx;
	
	r2 = ctx->radius * ctx->radius;
	sep2 =  (pos[0]-ctx->origin[0])*(pos[0]-ctx->origin[0]);
	sep2 += (pos[1]-ctx->origin[1])*(pos[1]-ctx->origin[1]);
	sep2 += (pos[2]-ctx->origin[2])*(pos[2]-ctx->origin[2]);
	if (sep2 < r2) {
		*inside = 1;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetType_Sphere"
PetscErrorCode GeometryObjectSetType_Sphere(GeometryObject go,double origin[],double radius)
{
	GeomTypeSphere ctx;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeomTypeSphere),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_GeomTypeSphere));CHKERRQ(ierr);
	ctx->origin[0] = origin[0];
	ctx->origin[1] = origin[1];
	ctx->origin[2] = origin[2];
	ctx->radius = radius;
	
	go->type                     = GeomType_Sphere;
	go->geom_point_inside        = GeometryObjectPointInside_Sphere;
	go->geom_transform_translate = GeometryObjectTransformTranslate_Sphere;
	go->geom_destroy             = GeometryObjectDestroy_Sphere;
	go->ctx = (void*)ctx;
	
	ierr = GeometryObjectSetFromOptions_Sphere(go);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------------------------- */
/* CYLINDER: implementation */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectDestroy_Cylinder"
PetscErrorCode GeometryObjectDestroy_Cylinder(GeometryObject go)
{
	GeomTypeCylinder ctx;
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeCylinder)go->ctx;
	ierr = PetscFree(ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetFromOptions_Cylinder"
PetscErrorCode GeometryObjectSetFromOptions_Cylinder(GeometryObject go)
{
	GeomTypeCylinder ctx;
	PetscReal val;
	PetscInt ival;
	PetscBool flg;
	char name[1024];
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeCylinder)go->ctx;
	
	sprintf(name,"%s_",go->name);
	ierr = PetscOptionsGetReal(name,"-Ox",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->x0[0] = val; }
	ierr = PetscOptionsGetReal(name,"-Oy",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->x0[1] = val; }
	ierr = PetscOptionsGetReal(name,"-Oz",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->x0[2] = val; }
	
	ierr = PetscOptionsGetReal(name,"-rad",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->radius = val; }
	
	ierr = PetscOptionsGetReal(name,"-len",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->length = val; }
	
	ierr = PetscOptionsGetInt(name,"-axis",&ival,&flg);CHKERRQ(ierr);
	if (flg) { ctx->axis = (GeomRotateAxis)ival; }
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectTransformTranslate_Cylinder"
PetscErrorCode GeometryObjectTransformTranslate_Cylinder(GeometryObject go,double shift[])
{
	GeomTypeCylinder ctx;
	
	
	ctx = (GeomTypeCylinder)go->ctx;
	
	ctx->x0[0] = ctx->x0[0] + shift[0];
	ctx->x0[1] = ctx->x0[1] + shift[1];
	ctx->x0[2] = ctx->x0[2] + shift[2];
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPointInside_Cylinder"
PetscErrorCode GeometryObjectPointInside_Cylinder(GeometryObject go,double pos[],int *inside)
{
	GeomTypeCylinder ctx;
	double r2,sep2;
	
	
	ctx = (GeomTypeCylinder)go->ctx;
	
	r2 = ctx->radius * ctx->radius;
	
	
	*inside = 0;
	switch (ctx->axis) {
			
		case ROTATE_AXIS_Z: // check in x-y
			
			if (pos[2] < ctx->x0[2])               { PetscFunctionReturn(0); }
			if (pos[2] > ctx->x0[2] + ctx->length) { PetscFunctionReturn(0); }
			
			sep2 =  (pos[0]-ctx->x0[0])*(pos[0]-ctx->x0[0]);
			sep2 += (pos[1]-ctx->x0[1])*(pos[1]-ctx->x0[1]);
			if (sep2 < r2) {
				*inside = 1;
			}
			
			break;
			
		case ROTATE_AXIS_Y: // check in x-z
			
			if (pos[1] < ctx->x0[1])               { PetscFunctionReturn(0); }
			if (pos[1] > ctx->x0[1] + ctx->length) { PetscFunctionReturn(0); }
			
			sep2 =  (pos[0]-ctx->x0[0])*(pos[0]-ctx->x0[0]);
			sep2 += (pos[2]-ctx->x0[2])*(pos[2]-ctx->x0[2]);
			if (sep2 < r2) {
				*inside = 1;
			}
			
			break;
			
		case ROTATE_AXIS_X: // check in y-z
			
			if (pos[0] < ctx->x0[0])               { PetscFunctionReturn(0); }
			if (pos[0] > ctx->x0[0] + ctx->length) { PetscFunctionReturn(0); }
			
			sep2  = (pos[1]-ctx->x0[1])*(pos[1]-ctx->x0[1]);
			sep2 += (pos[2]-ctx->x0[2])*(pos[2]-ctx->x0[2]);
			if (sep2 < r2) {
				*inside = 1;
			}
			
			break;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetType_Cylinder"
PetscErrorCode GeometryObjectSetType_Cylinder(GeometryObject go,double origin[],double radius,double length,GeomRotateAxis axis)
{
	GeomTypeCylinder ctx;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeomTypeCylinder),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_GeomTypeCylinder));CHKERRQ(ierr);
	ctx->x0[0]  = origin[0];
	ctx->x0[1]  = origin[1];
	ctx->x0[2]  = origin[2];
	ctx->radius = radius;
	ctx->length = length;
	ctx->axis   = axis;
	
	go->type                     = GeomType_Cylinder;
	go->geom_point_inside        = GeometryObjectPointInside_Cylinder;
	go->geom_transform_translate = GeometryObjectTransformTranslate_Cylinder;
	go->geom_destroy             = GeometryObjectDestroy_Cylinder;
	go->ctx = (void*)ctx;
	
	ierr = GeometryObjectSetFromOptions_Cylinder(go);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------------------------- */
/* SetOperation: Implementation */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectDestroy_SetOperation"
PetscErrorCode GeometryObjectDestroy_SetOperation(GeometryObject go)
{
	GeomTypeSetOperation ctx;
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeSetOperation)go->ctx;
	
	if (ctx->A) { ierr = GeometryObjectDestroy(&ctx->A);CHKERRQ(ierr); }
	if (ctx->B) { ierr = GeometryObjectDestroy(&ctx->B);CHKERRQ(ierr); }
	ierr = PetscFree(ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectTransformTranslate_SetOperation"
PetscErrorCode GeometryObjectTransformTranslate_SetOperation(GeometryObject go,double shift[])
{
	GeomTypeSetOperation ctx;
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeSetOperation)go->ctx;
	
	if (ctx->A) {
		ierr = GeometryObjectTransformTranslate(ctx->A,shift);CHKERRQ(ierr);
	}
	if (ctx->B) {
		ierr = GeometryObjectTransformTranslate(ctx->B,shift);CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPointInside_SetOperation"
PetscErrorCode GeometryObjectPointInside_SetOperation(GeometryObject go,double xp[],int *inside)
{
	GeomTypeSetOperation ctx;
	int is1,is2;
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeSetOperation)go->ctx;
	
	ierr = GeometryObjectPointInside(ctx->A,xp,&is1);CHKERRQ(ierr);
	ierr = GeometryObjectPointInside(ctx->B,xp,&is2);CHKERRQ(ierr);
	
	switch (ctx->set_type) {
			
		case GeomType_SetUnion:
			*inside = 0;
			if ( (is1 == 1) || (is2 == 1) ) {
				*inside = 1;
			}
			break;
			
		case GeomType_SetIntersection:
			*inside = 0;
			if ( (is1 == 1) && (is2 == 1) ) {
				*inside = 1;
			}
			break;
			
		case GeomType_SetComplement:
			*inside = 0;
			if (is1 == 1) {
				*inside = 1;
			}
			if ( (is1 == 1) && (is2 == 1) ) {
				*inside = 0;
			}
			if (is2 == 1) {
				*inside = 0;
			}
			break;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetType_SetOperation"
PetscErrorCode GeometryObjectSetType_SetOperation(GeometryObject go,GeomTypeSetOperator op_type,GeometryObject A,GeometryObject B)
{
	GeomTypeSetOperation ctx;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeomTypeSetOperation),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_GeomTypeSetOperation));CHKERRQ(ierr);
	ctx->set_type = op_type;
	ctx->A = A;
	ctx->B = B;
	if (A) { A->ref_cnt++; }
	if (B) { B->ref_cnt++; }
	
	go->type                     = GeomType_SetOperation;
	go->geom_point_inside        = GeometryObjectPointInside_SetOperation;
	go->geom_transform_translate = GeometryObjectTransformTranslate_SetOperation;
	go->geom_destroy             = GeometryObjectDestroy_SetOperation;
	go->ctx = (void*)ctx;
	
	PetscFunctionReturn(0);
}


