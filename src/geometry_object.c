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
 **    Filename:      geometry_object.c
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
	"GeomType_EllipticCylinder",
	"GeomType_Ellipsoid",
	"GeomType_InfLayer",
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
void PointTranslate(double xin[],double shift[],double xout[]);
void PointBackTranslate(double xin[],double shift[],double xout[]);
/* ------------------------------------------------------------------------------------------- */
/* API */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectCreate"
PetscErrorCode GeometryObjectCreate(const char name[],GeometryObject *G)
{
	GeometryObject go;
	int k;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeometryObject),&go);CHKERRQ(ierr);
	ierr = PetscMemzero(go,sizeof(struct _p_GeometryObject));CHKERRQ(ierr);
	
	asprintf(&go->name,"%s",name);
	go->type = GeomType_NULL;
	go->ctx = PETSC_NULL;
	go->n_rotations = 0;
	go->ref_cnt = 0;
	
	for (k=0; k<GEOM_SHAPE_MAX_ROTATIONS; k++) {
		go->rotation_angle[k] = 0.0;
		go->rotation_axis[k] = ROTATE_AXIS_UNDEFINED;
	}
	
	go->geom_point_inside = PETSC_NULL;
	go->geom_transform_translate = PETSC_NULL;
	
	go->centroid[0]=0.0;
	go->centroid[1]=0.0;
	go->centroid[2]=0.0;
	
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
		SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"No more rotations permitted");
	}
	
	angle2 = angle * M_PI / 180.0;
	
	go->rotation_axis[  go->n_rotations ] = dir;
	go->rotation_angle[ go->n_rotations ] = angle2;
	go->n_rotations++;
	
	PetscFunctionReturn(0);
}


/* point inside */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPointInside"
PetscErrorCode GeometryObjectPointInside(GeometryObject go,double pos[],int *inside)
{
	int    nr;
	double posin[3],posout[3];
	
	PointBackTranslate(pos,go->centroid,posout);
	
	for (nr=go->n_rotations-1; nr>=0; nr--) {
		posin[0] = posout[0];
		posin[1] = posout[1];
		posin[2] = posout[2];
		PointBackRotate(posin,go->rotation_axis[nr],go->rotation_angle[nr],posout);
	}
	
	/* call method */
	if (go->geom_point_inside) {
		*inside = 0;
		go->geom_point_inside(go,posout,inside);
	} else {
		SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"No method for go->geom_point_inside provided");
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetCentroid"
PetscErrorCode GeometryObjectSetCentroid(GeometryObject go,double cx[])
{
	go->centroid[0] = cx[0];
	go->centroid[1] = cx[1];
	go->centroid[2] = cx[2];
	
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
	SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"GeomType %s not located in list of available GeomTypes",name);
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
		PetscPrintf(PETSC_COMM_SELF,"[Warning] GeomObject with name %s was not found in list\n",name);
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectIdFindByName"
PetscErrorCode GeometryObjectIdFindByName(GeometryObject G[],const char name[],PetscInt *GoId)
{
	GeometryObject item;
	int i,v;
	
	*GoId = -1;
	
	i = 0;
	item = G[i];
	while (item != NULL) {
		v = strcmp(name,G[i]->name);
		if (v == 0) {
			*GoId = i;
			break;
		}
		i++;
		item = G[i];
	}
	if (*GoId == -1) {
		PetscPrintf(PETSC_COMM_SELF,"[Warning] GeomObject with name %s was not found in list\n",name);
	}
	PetscFunctionReturn(0);
}

void PointTranslate(double xin[],double shift[],double xout[])
{
	int i;
	for (i=0; i<3; i++){
		xout[i] = xin[i] + shift[i];
	}
}

void PointBackTranslate(double xin[],double shift[],double xout[])
{
	double shift2[3];
	int    i;
	for (i=0; i<3; i++){
		shift2[i] = -shift[i];
	}
	PointTranslate(xin,shift2,xout);
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
			
		case ROTATE_AXIS_UNDEFINED:
			PetscPrintf(PETSC_COMM_SELF,"ERROR: Cannot rotate as axis isn't defined");
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
	ierr = PetscOptionsGetReal(name,"-xc0",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[0] = val; }
	ierr = PetscOptionsGetReal(name,"-xc1",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[1] = val; }
	ierr = PetscOptionsGetReal(name,"-xc2",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[2] = val; }
	
	ierr = PetscOptionsGetReal(name,"-Lx0",&val,&flg);CHKERRQ(ierr);
	if (flg) { box->Lx[0] = val; }
	ierr = PetscOptionsGetReal(name,"-Lx1",&val,&flg);CHKERRQ(ierr);
	if (flg) { box->Lx[1] = val; }
	ierr = PetscOptionsGetReal(name,"-Lx2",&val,&flg);CHKERRQ(ierr);
	if (flg) { box->Lx[2] = val; }
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPointInside_Box"
PetscErrorCode GeometryObjectPointInside_Box(GeometryObject go,double pos[],int *inside)
{
	GeomTypeBox box;
	PetscErrorCode ierr;
	PetscInt i; 
	
	ierr = GeometryObjectGetContext_Box(go,&box);CHKERRQ(ierr);
	*inside = 0;
	
	for (i=0; i<3; i++){
		if ((pos[i]*pos[i]) > (box->Lx[i]*box->Lx[i]*0.25)) { PetscFunctionReturn(0); }
	}
		
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
	box->Lx[0] = Lx[0];
	box->Lx[1] = Lx[1];
	box->Lx[2] = Lx[2];
	go->centroid[0] = x0[0];
	go->centroid[1] = x0[1];
	go->centroid[2] = x0[2];
	
	go->type                     = GeomType_Box;
	go->geom_point_inside        = GeometryObjectPointInside_Box;
	go->geom_destroy             = GeometryObjectDestroy_Box;
	go->ctx = (void*)box;
	
	ierr = GeometryObjectSetFromOptions_Box(go);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetType_BoxCornerReference"
PetscErrorCode GeometryObjectSetType_BoxCornerReference(GeometryObject go,double x0[],double Lx[])
{
	double c0[3];
	PetscErrorCode ierr;
	
	c0[0] = x0[0] + 0.5 * Lx[0];
	c0[1] = x0[1] + 0.5 * Lx[1];
	c0[2] = x0[2] + 0.5 * Lx[2];
	ierr = GeometryObjectSetType_Box(go,c0,Lx);CHKERRQ(ierr);
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
	ierr = PetscOptionsGetReal(name,"-xc0",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[0] = val; }
	ierr = PetscOptionsGetReal(name,"-xc1",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[1] = val; }
	ierr = PetscOptionsGetReal(name,"-xc2",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[2] = val; }
	
	ierr = PetscOptionsGetReal(name,"-rad",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->radius = val; }
	
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
	sep2 =  pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
	
	if (sep2 < r2) {
		*inside = 1;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetType_Sphere"
PetscErrorCode GeometryObjectSetType_Sphere(GeometryObject go,double x0[],double radius)
{
	GeomTypeSphere ctx;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeomTypeSphere),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_GeomTypeSphere));CHKERRQ(ierr);
    
	go->centroid[0] = x0[0];
	go->centroid[1] = x0[1];
	go->centroid[2] = x0[2];

	ctx->radius = radius;
	
	go->type                     = GeomType_Sphere;
	go->geom_point_inside        = GeometryObjectPointInside_Sphere;
	go->geom_destroy             = GeometryObjectDestroy_Sphere;
	go->ctx = (void*)ctx;
	
	ierr = GeometryObjectSetFromOptions_Sphere(go);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------------- */
/* ELLIPSOID: implementation */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectDestroy_Ellipsoid"
PetscErrorCode GeometryObjectDestroy_Ellipsoid(GeometryObject go)
{
	GeomTypeEllipsoid ctx;
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeEllipsoid)go->ctx;
	ierr = PetscFree(ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetFromOptions_Ellipsoid"
PetscErrorCode GeometryObjectSetFromOptions_Ellipsoid(GeometryObject go)
{
	GeomTypeEllipsoid ctx;
	PetscReal val;
	PetscBool flg;
	char name[1024];
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeEllipsoid)go->ctx;
	
	sprintf(name,"%s_",go->name);
	ierr = PetscOptionsGetReal(name,"-xc0",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[0] = val; }
	ierr = PetscOptionsGetReal(name,"-xc1",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[1] = val; }
	ierr = PetscOptionsGetReal(name,"-xc2",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[2] = val; }
	
	ierr = PetscOptionsGetReal(name,"-ra",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->radia = val; }
	ierr = PetscOptionsGetReal(name,"-rb",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->radib = val; }
	ierr = PetscOptionsGetReal(name,"-rc",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->radic = val; }
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPointInside_Ellipsoid"
PetscErrorCode GeometryObjectPointInside_Ellipsoid(GeometryObject go,double pos[],int *inside)
{
	GeomTypeEllipsoid ctx;
	double ra2,rb2,rc2,sep2;
	
	*inside = 0;
	ctx = (GeomTypeEllipsoid)go->ctx;
	
	ra2 = ctx->radia * ctx->radia;
    rb2 = ctx->radib * ctx->radib;
    rc2 = ctx->radic * ctx->radic;
	
	sep2 =  pos[0]*pos[0]/ra2+pos[1]*pos[1]/rb2+pos[2]*pos[2]/rc2;
	
	if (sep2 < 1) {
		*inside = 1;
	}
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetType_Ellipsoid"
PetscErrorCode GeometryObjectSetType_Ellipsoid(GeometryObject go,double x0[],double radia,double radib,double radic)
{
	GeomTypeEllipsoid ctx;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeomTypeEllipsoid),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_GeomTypeEllipsoid));CHKERRQ(ierr);
    
    go->centroid[0] = x0[0];
	go->centroid[1] = x0[1];
	go->centroid[2] = x0[2];

	ctx->radia = radia	;
    ctx->radib = radib	;
    ctx->radic = radic	;
	
	go->type                     = GeomType_Ellipsoid;
	go->geom_point_inside        = GeometryObjectPointInside_Ellipsoid;
	go->geom_destroy             = GeometryObjectDestroy_Ellipsoid;
	go->ctx = (void*)ctx;
	
	ierr = GeometryObjectSetFromOptions_Ellipsoid(go);CHKERRQ(ierr);
	
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
	ierr = PetscOptionsGetReal(name,"-xc0",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[0] = val; }
	ierr = PetscOptionsGetReal(name,"-xc1",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[1] = val; }
	ierr = PetscOptionsGetReal(name,"-xc2",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[2] = val; }
	
	ierr = PetscOptionsGetReal(name,"-rad",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->radius = val; }
	
	ierr = PetscOptionsGetReal(name,"-len",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->length = val; }
	
	ierr = PetscOptionsGetInt(name,"-axis",&ival,&flg);CHKERRQ(ierr);
	if (flg) { ctx->axis = (GeomRotateAxis)ival; }
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPointInside_Cylinder"
PetscErrorCode GeometryObjectPointInside_Cylinder(GeometryObject go,double pos[],int *inside)
{
	GeomTypeCylinder ctx;
	double r2,l2;
	
	
	ctx = (GeomTypeCylinder)go->ctx;
	
	r2 = ctx->radius * ctx->radius;
	l2 = ctx->length * ctx->length * 0.25;
	
	*inside = 0;
	switch (ctx->axis) {
			
		case ROTATE_AXIS_Z: // check in x-y
			if ((pos[2]*pos[2]) > l2) { PetscFunctionReturn(0); }
			if (pos[0]*pos[0]+pos[1]*pos[1] > r2) { PetscFunctionReturn(0); } 
			break;
			
		case ROTATE_AXIS_Y: // check in x-z
			if ((pos[1]*pos[1]) > l2) { PetscFunctionReturn(0); }
			if (pos[0]*pos[0]+pos[2]*pos[2] > r2) { PetscFunctionReturn(0); } 
			break;			
			
		case ROTATE_AXIS_X: // check in y-z
			if ((pos[0]*pos[0]) > l2) { PetscFunctionReturn(0); }
			if (pos[1]*pos[1]+pos[2]*pos[2] > r2) { PetscFunctionReturn(0); } 
			break;			

		case ROTATE_AXIS_UNDEFINED:
			PetscPrintf(PETSC_COMM_SELF,"ERROR: Cannot rotate as axis isn't defined");
			break;
			
	}
	*inside = 1;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetType_Cylinder"
PetscErrorCode GeometryObjectSetType_Cylinder(GeometryObject go,double x0[],double radius,double length,GeomRotateAxis axis)
{
	GeomTypeCylinder ctx;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeomTypeCylinder),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_GeomTypeCylinder));CHKERRQ(ierr);
    go->centroid[0] = x0[0];
	go->centroid[1] = x0[1];
	go->centroid[2] = x0[2];

	ctx->radius = radius;
	ctx->length = length;
	ctx->axis   = axis;
	
	go->type                     = GeomType_Cylinder;
	go->geom_point_inside        = GeometryObjectPointInside_Cylinder;
	go->geom_destroy             = GeometryObjectDestroy_Cylinder;
	go->ctx = (void*)ctx;
	
	ierr = GeometryObjectSetFromOptions_Cylinder(go);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------------- */
/* Elliptic CYLINDER: implementation */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectDestroy_EllipticCylinder"
PetscErrorCode GeometryObjectDestroy_EllipticCylinder(GeometryObject go)
{
	GeomTypeEllipticCylinder ctx;
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeEllipticCylinder)go->ctx;
	ierr = PetscFree(ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetFromOptions_EllipticCylinder"
PetscErrorCode GeometryObjectSetFromOptions_EllipticCylinder(GeometryObject go)
{
	GeomTypeEllipticCylinder ctx;
	PetscReal val;
	PetscInt ival;
	PetscBool flg;
	char name[1024];
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeEllipticCylinder)go->ctx;
	
	sprintf(name,"%s_",go->name);
	ierr = PetscOptionsGetReal(name,"-xc0",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[0] = val; }
	ierr = PetscOptionsGetReal(name,"-xc1",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[1] = val; }
	ierr = PetscOptionsGetReal(name,"-xc2",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[2] = val; }
	
	ierr = PetscOptionsGetReal(name,"-ra",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->radia = val; }

	ierr = PetscOptionsGetReal(name,"-rb",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->radib = val; }
	
	ierr = PetscOptionsGetReal(name,"-len",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->length = val; }
	
	ierr = PetscOptionsGetInt(name,"-axis",&ival,&flg);CHKERRQ(ierr);
	if (flg) { ctx->axis = (GeomRotateAxis)ival; }
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPointInside_EllipticCylinder"
PetscErrorCode GeometryObjectPointInside_EllipticCylinder(GeometryObject go,double pos[],int *inside)
{
	GeomTypeEllipticCylinder ctx;
	double ra2,rb2;
	
	
	ctx = (GeomTypeEllipticCylinder)go->ctx;
	
	ra2 = ctx->radia * ctx->radia;
	rb2 = ctx->radib * ctx->radib;
	
	*inside = 0;
	switch (ctx->axis) {
			
		case ROTATE_AXIS_Z: // check in x-y
			if ((pos[2]*pos[2]) > (ctx->length*ctx->length*0.25)) { PetscFunctionReturn(0); }
			if ((pos[0]*pos[0])/ra2+pos[1]*pos[1]/rb2 > 1.0) { PetscFunctionReturn(0); } 
			break;
			
		case ROTATE_AXIS_Y: // check in x-z

			if ((pos[1]*pos[1]) > (ctx->length*ctx->length*0.25)) { PetscFunctionReturn(0); }
			if (pos[2]*pos[2]/ra2+pos[0]*pos[0]/rb2 > 1.0) { PetscFunctionReturn(0); } 
			break;			
			
		case ROTATE_AXIS_X: // check in y-z
			if ((pos[0]*pos[0]) > (ctx->length*ctx->length*0.25)) { PetscFunctionReturn(0); }
			if (pos[1]*pos[1]/ra2+pos[2]*pos[2]/rb2 > 1.0) { PetscFunctionReturn(0); } 
			break;	
			
		case ROTATE_AXIS_UNDEFINED:
			PetscPrintf(PETSC_COMM_SELF,"ERROR: Cannot rotate as axis isn't defined");
			break;
			
	}
	*inside = 1;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetType_EllipticCylinder"
PetscErrorCode GeometryObjectSetType_EllipticCylinder(GeometryObject go,double x0[],double radia,double radib,double length,GeomRotateAxis axis)
{
	GeomTypeEllipticCylinder ctx;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeomTypeEllipticCylinder),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_GeomTypeEllipticCylinder));CHKERRQ(ierr);
    
    go->centroid[0] = x0[0];
	go->centroid[1] = x0[1];
	go->centroid[2] = x0[2];

	ctx->radia  = radia;
	ctx->radib  = radib;
	ctx->length = length;
	ctx->axis   = axis;
	
	go->type                     = GeomType_EllipticCylinder;
	go->geom_point_inside        = GeometryObjectPointInside_EllipticCylinder;
	go->geom_destroy             = GeometryObjectDestroy_EllipticCylinder;
	go->ctx = (void*)ctx;
	
	ierr = GeometryObjectSetFromOptions_EllipticCylinder(go);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------------------------- */


/* ------------------------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------------------------- */
/* Elliptic CYLINDER: implementation */
#undef __FUNCT__
#define __FUNCT__ "GeometryObjectDestroy_InfLayer"
PetscErrorCode GeometryObjectDestroy_InfLayer(GeometryObject go)
{
	GeomTypeInfLayer ctx;
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeInfLayer)go->ctx;
	ierr = PetscFree(ctx);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetFromOptions_InfLayer"
PetscErrorCode GeometryObjectSetFromOptions_InfLayer(GeometryObject go)
{
	GeomTypeInfLayer ctx;
	PetscReal val;
	PetscInt ival;
	PetscBool flg;
	char name[1024];
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeInfLayer)go->ctx;
	
	sprintf(name,"%s_",go->name);
	ierr = PetscOptionsGetReal(name,"-xc0",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[0] = val; }
	ierr = PetscOptionsGetReal(name,"-xc1",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[1] = val; }
	ierr = PetscOptionsGetReal(name,"-xc2",&val,&flg);CHKERRQ(ierr);
	if (flg) { go->centroid[2] = val; }
	
	ierr = PetscOptionsGetReal(name,"-thick",&val,&flg);CHKERRQ(ierr);
	if (flg) { ctx->thickness = val; }
	
	ierr = PetscOptionsGetInt(name,"-axis",&ival,&flg);CHKERRQ(ierr);
	if (flg) { ctx->axis = (GeomRotateAxis)ival; }
	
	PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "GeometryObjectPointInside_InfLayer"
PetscErrorCode GeometryObjectPointInside_InfLayer(GeometryObject go,double pos[],int *inside)
{
	GeomTypeInfLayer ctx;
		
	
	ctx = (GeomTypeInfLayer)go->ctx;
	*inside = 0;
	switch (ctx->axis) {
			
		case ROTATE_AXIS_Z: // check in x-y
			if ((pos[2]*pos[2]) > (ctx->thickness*ctx->thickness*0.25)) { PetscFunctionReturn(0); }
			break;
			
		case ROTATE_AXIS_Y: // check in x-z
			if ((pos[1]*pos[1]) > (ctx->thickness*ctx->thickness*0.25)) { PetscFunctionReturn(0); }
			break;			
			
		case ROTATE_AXIS_X: // check in y-z
			if ((pos[0]*pos[0]) > (ctx->thickness*ctx->thickness*0.25)) { PetscFunctionReturn(0); }
			break;			
	
		case ROTATE_AXIS_UNDEFINED:
			PetscPrintf(PETSC_COMM_SELF,"ERROR: Cannot rotate as axis isn't defined");
			break;
			
	}
	*inside = 1;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetType_InfLayer"
PetscErrorCode GeometryObjectSetType_InfLayer(GeometryObject go,double x0[],double thickness,GeomRotateAxis axis)
{
	GeomTypeInfLayer ctx;
	PetscErrorCode ierr;
	
	
	ierr = PetscMalloc(sizeof(struct _p_GeomTypeInfLayer),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_GeomTypeInfLayer));CHKERRQ(ierr);
    go->centroid[0] = x0[0];
	go->centroid[1] = x0[1];
	go->centroid[2] = x0[2];

	ctx->thickness  = thickness;
	ctx->axis       = axis;
	
	go->type                     = GeomType_InfLayer;
	go->geom_point_inside        = GeometryObjectPointInside_InfLayer;
	go->geom_destroy             = GeometryObjectDestroy_InfLayer;
	go->ctx = (void*)ctx;
	
	ierr = GeometryObjectSetFromOptions_InfLayer(go);CHKERRQ(ierr);
	
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
#define __FUNCT__ "GeometryObjectPointInside_SetOperation"
PetscErrorCode GeometryObjectPointInside_SetOperation(GeometryObject go,double xp[],int *inside)
{
	GeomTypeSetOperation ctx;
	int is1,is2;
	double xpc[3];
	PetscErrorCode ierr;
	
	
	ctx = (GeomTypeSetOperation)go->ctx;
	
	/* Undo the coordinate shift performed by the first call to GeometryObjectPointInside() which was applied to the set */
	PointTranslate(xp,go->centroid,xpc);	

	ierr = GeometryObjectPointInside(ctx->A,xpc,&is1);CHKERRQ(ierr);
	ierr = GeometryObjectPointInside(ctx->B,xpc,&is2);CHKERRQ(ierr);
	
	switch (ctx->operator_type) {
			
		case GeomSet_Union:
			*inside = 0;
			if ( (is1 == 1) || (is2 == 1) ) {
				*inside = 1;
			}
			break;
			
		case GeomSet_Intersection:
			*inside = 0;
			if ( (is1 == 1) && (is2 == 1) ) {
				*inside = 1;
			}
			break;
			
		case GeomSet_Complement:
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
PetscErrorCode GeometryObjectSetType_SetOperation(GeometryObject go,GeomSetOperatorType type,double x0[],GeometryObject A,GeometryObject B)
{
	GeomTypeSetOperation ctx;
	PetscErrorCode ierr;
	
	go->centroid[0] = x0[0];
	go->centroid[1] = x0[1];
	go->centroid[2] = x0[2];
	
	ierr = PetscMalloc(sizeof(struct _p_GeomTypeSetOperation),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_GeomTypeSetOperation));CHKERRQ(ierr);
	ctx->operator_type = type;
	ctx->A = A;
	ctx->B = B;
	if (A) { A->ref_cnt++; }
	if (B) { B->ref_cnt++; }
	
	go->type                     = GeomType_SetOperation;
	go->geom_point_inside        = GeometryObjectPointInside_SetOperation;
	go->geom_destroy             = GeometryObjectDestroy_SetOperation;
	go->ctx = (void*)ctx;
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GeometryObjectSetType_SetOperationDefault"
PetscErrorCode GeometryObjectSetType_SetOperationDefault(GeometryObject go,GeomSetOperatorType type,GeometryObject A,GeometryObject B)
{
	GeomTypeSetOperation ctx;
	PetscErrorCode ierr;
	
	go->centroid[0] = 0.5 * ( A->centroid[0] + B->centroid[0] );
	go->centroid[1] = 0.5 * ( A->centroid[1] + B->centroid[1] );
	go->centroid[2] = 0.5 * ( A->centroid[2] + B->centroid[2] );
	
	ierr = PetscMalloc(sizeof(struct _p_GeomTypeSetOperation),&ctx);CHKERRQ(ierr);
	ierr = PetscMemzero(ctx,sizeof(struct _p_GeomTypeSetOperation));CHKERRQ(ierr);
	ctx->operator_type = type;
	ctx->A = A;
	ctx->B = B;
	if (A) { A->ref_cnt++; }
	if (B) { B->ref_cnt++; }
	
	go->type                     = GeomType_SetOperation;
	go->geom_point_inside        = GeometryObjectPointInside_SetOperation;
	go->geom_destroy             = GeometryObjectDestroy_SetOperation;
	go->ctx = (void*)ctx;
	
	PetscFunctionReturn(0);
}


