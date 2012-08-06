

#ifndef __ptatin3d_MATERIAL_POINT_POINT_LOCATION_H__
#define __ptatin3d_MATERIAL_POINT_POINT_LOCATION_H__

typedef enum { _FALSE=0, _TRUE=1 } Truth;

void InverseMappingDomain_3dQ2( 
															 double tolerance, int max_its,
															 Truth use_nonzero_guess, 
															 Truth monitor, Truth log,
															 const double coords[], const int mx, const int my, const int mz, const int element[],
															 int np, MPntStd marker[] );

#endif
