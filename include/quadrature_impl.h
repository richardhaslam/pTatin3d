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
 **    filename:   quadrature_impl.h
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

#ifndef __private_ptatin_quadrature_impl_h__
#define __private_ptatin_quadrature_impl_h__

#include "petsc.h"
#include "../element_type_Q2.h"

struct _p_Quadrature {
	PetscInt       dim;
	QuadratureType type; /* line (dim=2), surface(dim=3), vol(dim=2,3) */
	PetscInt    npoints;
	PetscReal  *q_xi_coor;
	PetscReal  *q_weight;
	PetscInt   n_elements;
	DataBucket properties_db;
};

struct _p_SurfaceQuadrature {
	ConformingElementFamily e;
	HexElementFace face_id;
	/* quadrature */
	PetscInt    ngp;
	QPoint2d    gp2[9]; /* s,t coordinates */
	QPoint3d    gp3[9]; /* xi,eta,zeta coordinates */
	PetscInt    nfaces;
	PetscInt    *element_list; /* list of cells connected to the face */
	DataBucket  properties_db;
};	


#endif
