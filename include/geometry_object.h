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
 **    filename:   geometry_object.h
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


#ifndef __geometry_object_h__
#define __geometry_object_h__

typedef enum {
	GeomType_Box              = 0,
	GeomType_Cylinder         = 1,
	GeomType_Sphere           = 2,
	GeomType_EllipticCylinder = 3,
	GeomType_Ellipsoid        = 4,
	GeomType_InfLayer         = 5,
	GeomType_SetOperation     = 6,
	GeomType_HalfSpace        = 7,
	GeomType_NULL
} GeomType;
extern const char *GeomTypeNames[];

typedef enum {
	GeomSet_Union        = 0,
	GeomSet_Intersection = 1,
	GeomSet_Complement   = 2,
	GeomSet_Undefined,
} GeomSetOperatorType;
extern const char *GeomTypeSetOperatorNames[];

typedef enum {
	ROTATE_AXIS_X = 0,
	ROTATE_AXIS_Y = 1,
	ROTATE_AXIS_Z = 2,
	ROTATE_AXIS_UNDEFINED
} GeomRotateAxis;
extern const char *GeomRotateAxisNames[];

typedef enum {
	SIGN_POSITIVE = 0,
	SIGN_NEGATIVE = 1,
	SIGN_UNDEFINED
} GeomSign;
extern const char *GeomSignNames[];


#define GEOM_SHAPE_MAX_ROTATIONS 25

typedef struct _p_GeometryObject *GeometryObject;
struct _p_GeometryObject {
	char     *name;
	GeomType type;
	int      ref_cnt;
	void     *ctx;
	double   centroid[3];
	int            n_rotations;
	double         rotation_angle[GEOM_SHAPE_MAX_ROTATIONS];
	GeomRotateAxis rotation_axis[GEOM_SHAPE_MAX_ROTATIONS];
	/* operations */
	PetscErrorCode (*geom_point_inside)(GeometryObject,double*,int*);
	PetscErrorCode (*geom_transform_translate)(GeometryObject,double*);
	PetscErrorCode (*geom_destroy)(GeometryObject);
};


/*
 Geometry specific data type implementations can be defined here.
 These should be made private so no one tries to directly access members from these structs.
*/
typedef struct _p_GeomTypeBox *GeomTypeBox;
struct _p_GeomTypeBox {
	double Lx[3];
};

typedef struct _p_GeomTypeCylinder *GeomTypeCylinder;
struct _p_GeomTypeCylinder {
	GeomRotateAxis axis;
	double radius,length;
};

typedef struct _p_GeomTypeSphere *GeomTypeSphere;
struct _p_GeomTypeSphere {
	double radius;
};

typedef struct _p_GeomTypeEllipticCylinder *GeomTypeEllipticCylinder;
struct _p_GeomTypeEllipticCylinder {
    GeomRotateAxis axis;
	double radia,radib,length;
};

typedef struct _p_GeomTypeInfLayer *GeomTypeInfLayer;
struct _p_GeomTypeInfLayer {
    GeomRotateAxis axis;
	double thickness;
};

typedef struct _p_GeomTypeEllipsoid *GeomTypeEllipsoid;
struct _p_GeomTypeEllipsoid {
	double radia,radib,radic;
};

typedef struct _p_GeomTypeSetOperation *GeomTypeSetOperation;
struct _p_GeomTypeSetOperation {
	GeomSetOperatorType operator_type;
	GeometryObject A,B;
};

typedef struct _p_GeomTypeHalfSpace *GeomTypeHalfSpace;
struct _p_GeomTypeHalfSpace {
	GeomSign       sign;
	GeomRotateAxis axis;
};

/*
 API
*/
PetscErrorCode GeometryObjectCreate(const char name[],GeometryObject *G);
PetscErrorCode GeometryObjectDestroy(GeometryObject *G);
PetscErrorCode GeometryObjectView(GeometryObject G);

PetscErrorCode GeometryObjectRotate(GeometryObject go,GeomRotateAxis dir,double angle);
PetscErrorCode GeometryObjectPointInside(GeometryObject go,double pos[],int *inside);
PetscErrorCode GeometryObjectSetCentroid(GeometryObject go,double cx[]);

/*
 Specific constructors for each implementation
*/
PetscErrorCode GeometryObjectSetType_Box(GeometryObject go,double x0[],double Lx[]);
PetscErrorCode GeometryObjectSetType_BoxCornerReference(GeometryObject go,double x0[],double Lx[]);

PetscErrorCode GeometryObjectSetType_SetOperation(GeometryObject go,GeomSetOperatorType type,double x0[],GeometryObject A,GeometryObject B);
PetscErrorCode GeometryObjectSetType_SetOperationDefault(GeometryObject go,GeomSetOperatorType type,GeometryObject A,GeometryObject B);

PetscErrorCode GeometryObjectSetType_Sphere(GeometryObject go,double origin[],double radius);

PetscErrorCode GeometryObjectSetType_Cylinder(GeometryObject go,double x0[],double radius,double Lx,GeomRotateAxis axis);

PetscErrorCode GeometryObjectSetType_EllipticCylinder(GeometryObject go,double x0[],double radia,double radib,double Lx,GeomRotateAxis axis);

PetscErrorCode GeometryObjectSetType_Ellipsoid(GeometryObject go,double x0[],double radia,double radib,double radic);

PetscErrorCode GeometryObjectSetType_InfLayer(GeometryObject go,double x0[],double Lx,GeomRotateAxis axis);

PetscErrorCode GeometryObjectSetType_HalfSpace(GeometryObject go,double x0[],GeomSign sign,GeomRotateAxis axis);

PetscErrorCode GeomTypeNameGetId(const char name[],int *id);
PetscErrorCode GeometryObjectFindByName(GeometryObject G[],const char name[],GeometryObject *g);
PetscErrorCode GeometryObjectIdFindByName(GeometryObject G[],const char name[],PetscInt *GoId);

#endif
