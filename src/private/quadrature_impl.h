
#ifndef __private_ptatin_quadrature_impl_h__
#define __private_ptatin_quadrature_impl_h__

struct _p_Quadrature {
	PetscInt       dim;
	QuadratureType type; /* line (dim=2), surface(dim=3), vol(dim=2,3) */
	PetscInt    npoints;
	PetscReal  *q_xi_coor;
	PetscReal  *q_weight;
	PetscInt   n_elements;
	DataBucket properties_db;
};

#endif
