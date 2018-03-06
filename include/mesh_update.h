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
 **    filename:   mesh_update.h
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

#ifndef __ptatin3d_mesh_update_h__
#define __ptatin3d_mesh_update_h__

PetscErrorCode DMDABilinearizeQ2Elements(DM dau);
PetscErrorCode UpdateMeshGeometry_FullLagrangian(DM dav,Vec velocity,PetscReal step);
PetscErrorCode UpdateMeshGeometry_VerticalLagrangianSurfaceRemesh(DM dav,Vec velocity,PetscReal step);
PetscErrorCode UpdateMeshGeometry_FullLagrangianWithVerticalSurfaceRemesh(DM dav,Vec velocity,PetscReal step);
PetscErrorCode UpdateMeshGeometry_DecoupledHorizontalVerticalMeshMovement(DM dav,Vec velocity,PetscReal step);
PetscErrorCode UpdateMeshGeometry_FullLag_ResampleJMax_RemeshJMIN2JMAX(DM dav,Vec vel_fluid,Vec vel_mesh,PetscReal step);
PetscErrorCode UpdateMeshGeometry_ComputeSurfaceCourantTimestep(DM dav,Vec velocity,PetscReal vert_displacement_max,PetscReal *step);
PetscErrorCode UpdateMeshGeometry_ApplyDiffusionJMAX(DM dav,PetscReal diffusivity,PetscReal timespan,PetscBool dirichlet_east,PetscBool dirichlet_west,PetscBool dirichlet_front,PetscBool dirichlet_back,PetscBool only_update_surface);

#endif

