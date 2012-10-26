#ifndef __ptatin3d_mesh_quality_metrics_h__
#define __ptatin3d_mesh_quality_metrics_h__

typedef enum {
	MESH_QUALITY_ASPECT_RATIO = 1,
	MESH_QUALITY_DISTORTION,
	MESH_QUALITY_DIAGONAL_RATIO
} MeshQualityMeasure;

PetscErrorCode DMDAComputeMeshQualityMetrics(DM dm,MeshQualityMeasure measure,PetscReal *value);

#endif
