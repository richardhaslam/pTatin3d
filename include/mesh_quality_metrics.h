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
 **    filename:   mesh_quality_metrics.h
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
#ifndef __ptatin3d_mesh_quality_metrics_h__
#define __ptatin3d_mesh_quality_metrics_h__

typedef enum {
  MESH_QUALITY_ASPECT_RATIO = 1,
  MESH_QUALITY_DISTORTION,
  MESH_QUALITY_DIAGONAL_RATIO,
    MESH_QUALITY_VERTEX_ANGLE,
    MESH_QUALITY_FACE_AREA_RATIO
} MeshQualityMeasure;

#define Q2_CELL_NODE_CENTER 13

#define Q2_FACE_NODE_WEST   12
#define Q2_FACE_NODE_EAST   14
#define Q2_FACE_NODE_SOUTH  10
#define Q2_FACE_NODE_NORTH  16
#define Q2_FACE_NODE_BACK   4
#define Q2_FACE_NODE_FRONT  22

#define Q2_VERTEX_0   0
#define Q2_VERTEX_1   2
#define Q2_VERTEX_2   6
#define Q2_VERTEX_3   8
#define Q2_VERTEX_4   18
#define Q2_VERTEX_5   20
#define Q2_VERTEX_6   24
#define Q2_VERTEX_7   26

typedef enum {
  Q2_FACE_WEST = 1,
  Q2_FACE_EAST,
  Q2_FACE_SOUTH,
  Q2_FACE_NORTH,
  Q2_FACE_BACK,
  Q2_FACE_FRONT
} Q2Face;

PetscErrorCode DMDAComputeMeshQualityMetric(DM dm,MeshQualityMeasure measure,PetscReal *value);
PetscErrorCode DMDAComputeMeshQualityMetricList(DM dm,const PetscInt nmeasures,MeshQualityMeasure measure[],PetscReal value[]);
PetscErrorCode DMDMeshQualityMetricGetInfo(DM dm,MeshQualityMeasure measure,PetscReal range[],PetscReal valid_range[]);

#endif

