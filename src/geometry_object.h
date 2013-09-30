#ifndef __geometry_object_h__
#define __geometry_object_h__

typedef enum { 
	GeomType_Box=0, 
	GeomType_Cylinder, 
	GeomType_Sphere,
	GeomType_SetOperation,
	GeomType_NULL
} GeomType;
extern const char *GeomTypeNames[];

typedef enum { 
	GeomType_SetUnion=0, 
	GeomType_SetIntersection, 
	GeomType_SetComplement,
} GeomTypeSetOperator;
extern const char *GeomTypeSetOperatorNames[];


typedef enum { 
	ROTATE_AXIS_X = 0,
	ROTATE_AXIS_Y = 1,
	ROTATE_AXIS_Z = 2
} GeomRotateAxis;
extern const char *GeomRotateAxisNames[];



#define GEOM_SHAPE_MAX_ROTATIONS 25

#define GEOM_SHAPE_REGION_VALUE_NULL -1.0e32
#define GEOM_SHAPE_REGION_INDEX_NULL -69


typedef struct _p_GeometryObject *GeometryObject;
struct _p_GeometryObject {
	char     *name;
	GeomType type;
	int      ref_cnt;
	void     *ctx;
	//double   centroid[3];
	int      region_index;
	double   value;
	int            n_rotations;
	double         rotation_angle[GEOM_SHAPE_MAX_ROTATIONS];
	GeomRotateAxis rotation_axis[GEOM_SHAPE_MAX_ROTATIONS];
	/* operations */
	PetscErrorCode (*geom_point_inside)(GeometryObject,double*,int*);
	PetscErrorCode (*geom_transform_translate)(GeometryObject,double*);
	PetscErrorCode (*geom_destroy)(GeometryObject);
  /* user methods */
	//int    (*evaluate_region_index)(double*);
	//double (*evaluate_region_value)(double*);
	double (*evaluate_region_function)(double*);
};


/* 
 Geometry specific data type implementations can be defined here.
 These should be made private so no one tries to directly access members from these structs.
*/
typedef struct _p_GeomTypeBox *GeomTypeBox;
struct _p_GeomTypeBox {
	double x0[3],Lx[3];
};

typedef struct _p_GeomTypeCylinder *GeomTypeCylinder;
struct _p_GeomTypeCylinder {
	GeomRotateAxis axis;
	double x0[3],radius,length;
};

typedef struct _p_GeomTypeSphere *GeomTypeSphere;
struct _p_GeomTypeSphere {
	double origin[3],radius;
};

typedef struct _p_GeomTypeSetOperation *GeomTypeSetOperation;
struct _p_GeomTypeSetOperation {
	int set_type;
	GeometryObject A,B;
};

/* 
 API
*/
PetscErrorCode GeometryObjectCreate(const char name[],int region_index,double value,double (*fp)(double*),GeometryObject *G);
PetscErrorCode GeometryObjectDestroy(GeometryObject *G);

PetscErrorCode GeometryObjectRotate(GeometryObject go,GeomRotateAxis dir,double angle);
PetscErrorCode GeometryObjectTransformTranslate(GeometryObject go,double shift[]);
PetscErrorCode GeometryObjectPointInside(GeometryObject go,double pos[],int *inside);
PetscErrorCode GeometryObjectEvaluateRegionIndex(GeometryObject go,double pos[],int *region);
PetscErrorCode GeometryObjectEvaluateRegionValue(GeometryObject go,double pos[],double *value);
PetscErrorCode GeometryObjectEvaluateRegionFunction(GeometryObject go,double pos[],double *value);

/* 
 Specific constructors for each implementation
*/
PetscErrorCode GeometryObjectSetType_Box(GeometryObject go,double x0[],double Lx[]);
PetscErrorCode GeometryObjectSetType_SetOperation(GeometryObject go,GeomTypeSetOperator type,GeometryObject A,GeometryObject B);
PetscErrorCode GeometryObjectSetType_Sphere(GeometryObject go,double origin[],double radius);
PetscErrorCode GeometryObjectSetType_Cylinder(GeometryObject go,double x0[],double radius,double Lx,GeomRotateAxis axis);

PetscErrorCode GeomTypeNameGetId(const char name[],int *id);
PetscErrorCode GeometryObjectFindByName(GeometryObject G[],const char name[],GeometryObject *g);

#endif
