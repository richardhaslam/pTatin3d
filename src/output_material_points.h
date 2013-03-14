
#ifndef __ptatin_output_material_points_h__
#define __ptatin_output_material_points_h__

typedef enum { 
	MPV_region=0, 
	MPV_viscosity, 
	MPV_density, 
	MPV_plastic_strain, 
	MPV_yield_indicator, 
	MPV_diffusivity, 
	MPV_heat_source 
} MaterialPointVariable;

PetscErrorCode pTatinOutputParaViewMarkerFields_VTS(DM dau,DataBucket material_points,const int nvars,const MaterialPointVariable vars[],const char name[]);
PetscErrorCode pTatinOutputParaViewMarkerFields_PVTS(DM dau,const int nvars,const MaterialPointVariable vars[],const char prefix[],const char name[]);
PetscErrorCode pTatinOutputParaViewMarkerFields(DM pack,DataBucket material_points,const int nvars,const MaterialPointVariable vars[],const char path[],const char prefix[]);
PetscErrorCode pTatin3d_ModelOutput_MarkerCellFields(pTatinCtx ctx,const int nvars,const MaterialPointVariable vars[],const char prefix[]);

#endif
