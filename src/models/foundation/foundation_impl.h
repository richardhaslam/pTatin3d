
#ifndef __foundation_impl_h__
#define __foundation_impl_h__

#include "foundation.h"

struct _p_Foundation {
	pTatinCtx              ptatin_ctx;
  cJSON                  *jfile,*root;
  FND_MeshGeomType       mesh_geom_type;
  FND_MeshRefinementType mesh_ref_type;
  FND_MaterialRegionType material_region_type;
  FoundationUserVars     user;
};

#endif
