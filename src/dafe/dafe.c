
#include <petsc.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petsc/private/petscimpl.h>
#include <petsc/private/dmimpl.h>
#include <private/dafeimpl.h>

const char *DAFETypeNames[] = { "DAFE_Q2","DAFE_Q1","DAFE_Q2sub","DAFE_Q1sub","DAFE_P1sub","DAFE_P0sub",0};


#undef __FUNCT__
#define __FUNCT__ "DAFECheckType"
PetscErrorCode DAFECheckType(DM dm)
{
  PetscBool isdafe = PETSC_FALSE;
  PetscFunctionBegin;
  PetscObjectTypeCompare((PetscObject)dm,"dafe",&isdafe);
  if (!isdafe) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Only valid if DM is of type DAFE");
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DMCreateGlobalVector_DAFE"
PetscErrorCode DMCreateGlobalVector_DAFE(DM dm,Vec *vec)
{
  DM_DAFE *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMCreateGlobalVector(fe->da,vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreateLocalVector_DAFE"
PetscErrorCode DMCreateLocalVector_DAFE(DM dm,Vec *vec)
{
  DM_DAFE *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMCreateLocalVector(fe->da,vec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDestroy_DAFE"
PetscErrorCode DMDestroy_DAFE(DM dm)
{
  DM_DAFE *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  
  if (fe->ops) PetscFree(fe->ops);
  if (fe->da) { ierr = DMDestroy(&fe->da);CHKERRQ(ierr); }
  if (fe->parent_dafe) { ierr = DMDestroy(&fe->parent_dafe);CHKERRQ(ierr); }
  if (fe->lmx_range) PetscFree(fe->lmx_range);
  if (fe->lmy_range) PetscFree(fe->lmy_range);
  if (fe->lmz_range) PetscFree(fe->lmz_range);

  if (fe->corner_imin_range) PetscFree(fe->corner_imin_range);
  if (fe->corner_jmin_range) PetscFree(fe->corner_jmin_range);
  if (fe->corner_kmin_range) PetscFree(fe->corner_kmin_range);

  if (fe->element_node_map) {
    if (fe->element_basis_map == fe->element_node_map) {
      PetscFree(fe->element_node_map);
      fe->element_node_map = NULL;
      fe->element_basis_map = NULL;
    }
  }
  
  if (fe->element_node_map) PetscFree(fe->element_node_map);
  if (fe->element_basis_map) PetscFree(fe->element_basis_map);
  
  ierr = PetscFree(fe);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMSetup_DAFE"
PetscErrorCode DMSetup_DAFE(DM dm)
{
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMView_DAFE"
PetscErrorCode DMView_DAFE(DM dm, PetscViewer viewer)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscBool      iascii,ibinary,ishdf5,isvtk;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm,DM_CLASSID,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,2);
  ierr = PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &iascii);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERBINARY,&ibinary);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERVTK,   &isvtk);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERHDF5,  &ishdf5);CHKERRQ(ierr);
  if (iascii) {
    PetscInt Mp,Np,Pp;
    
    ierr = DMDAGetInfo(fe->da,0,0,0,0,&Mp,&Np,&Pp,0,0, 0,0,0, 0);CHKERRQ(ierr);

    PetscViewerASCIIPushTab(viewer);
    PetscViewerASCIIPrintf(viewer,"Dimension: %D\n",dm->dim);
    PetscViewerASCIIPrintf(viewer,"DAFEType: %s\n",DAFETypeNames[(int)fe->dafetype]);
    PetscViewerASCIIPrintf(viewer,"Global elements: %D x %D x %D\n",fe->mx,fe->my,fe->mz);
    PetscViewerASCIIPrintf(viewer,"Basis per element: %D\n",fe->nbasis);
    PetscViewerASCIIPrintf(viewer,"Nodes per element: %D\n",fe->nodes_per_el);
    PetscViewerASCIIPrintf(viewer,"DOFs: %D\n",fe->ncomponents);
    PetscViewerASCIIPrintf(viewer,"Partition: %D x %D x %D\n",Mp,Np,Pp);
    if (fe->parent_dafe) {
      PetscViewerASCIIPrintf(viewer,"[ParentDMDAFE]\n");
      ierr = DMView(fe->parent_dafe,viewer);CHKERRQ(ierr);
    }
    PetscViewerASCIIPopTab(viewer);
  } else if (ibinary) {
    SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"NO binary support");
  } else if (ishdf5) {
#if defined(PETSC_HAVE_HDF5)
    SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"NO HDF5 support");
#else
    SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"HDF5 not supported. Please reconfigure using --download-hdf5");
#endif
  } else if (isvtk) {
    Vec dummy;
    
    ierr = DMCreateGlobalVector(fe->da,&dummy);CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)dummy,"empty_field");CHKERRQ(ierr);
    ierr = VecView(dummy,viewer);CHKERRQ(ierr);
    ierr = VecDestroy(&dummy);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMCreate_DAFE"
PetscErrorCode DMCreate_DAFE(DM dm)
{
  DM_DAFE        *fe;
  DAFEOps        ops;
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  ierr     = PetscNewLog(dm,&fe);CHKERRQ(ierr);
  dm->data = fe;

  ierr = PetscNewLog(dm,&ops);CHKERRQ(ierr);
  fe->ops = ops;
  
  fe->dafetype = DAFE_Q1;
  fe->ncomponents = 0;
  
  dm->dim  = 0;
  dm->ops->view                            = DMView_DAFE;
  dm->ops->load                            = NULL;
  dm->ops->setfromoptions                  = NULL;
  dm->ops->clone                           = NULL;
  dm->ops->setup                           = DMSetup_DAFE;
  dm->ops->createdefaultsection            = NULL;
  dm->ops->createdefaultconstraints        = NULL;
  dm->ops->createglobalvector              = DMCreateGlobalVector_DAFE;
  dm->ops->createlocalvector               = DMCreateLocalVector_DAFE;
  dm->ops->getlocaltoglobalmapping         = NULL;
  dm->ops->createfieldis                   = NULL;
  dm->ops->createcoordinatedm              = NULL;
  dm->ops->getcoloring                     = NULL;
  dm->ops->creatematrix                    = NULL;
  dm->ops->createinterpolation             = NULL;
  dm->ops->getaggregates                   = NULL;
  dm->ops->getinjection                    = NULL;
  dm->ops->refine                          = NULL;
  dm->ops->coarsen                         = NULL;
  dm->ops->refinehierarchy                 = NULL;
  dm->ops->coarsenhierarchy                = NULL;
  dm->ops->globaltolocalbegin              = NULL;
  dm->ops->globaltolocalend                = NULL;
  dm->ops->localtoglobalbegin              = NULL;
  dm->ops->localtoglobalend                = NULL;
  dm->ops->destroy                         = DMDestroy_DAFE;
  dm->ops->createsubdm                     = NULL;
  dm->ops->getdimpoints                    = NULL;
  dm->ops->locatepoints                    = NULL;
  
  PetscFunctionReturn(0);
}

/* setters */
#undef __FUNCT__
#define __FUNCT__ "DMDAFESetType"
PetscErrorCode DMDAFESetType(DM dm,DAFEType type)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);
  fe->dafetype = type;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAFESetNumberOfComponents"
PetscErrorCode DMDAFESetNumberOfComponents(DM dm,PetscInt nc)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);
  fe->ncomponents = nc;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAFESetUniformCoordinates"
PetscErrorCode DMDAFESetUniformCoordinates(DM dm,PetscReal xmin,PetscReal xmax,PetscReal ymin,PetscReal ymax,PetscReal zmin,PetscReal zmax)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(fe->da,xmin,xmax,ymin,ymax,zmin,zmax);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* getters */
#undef __FUNCT__
#define __FUNCT__ "DMDAFEGetNumberOfComponents"
PetscErrorCode DMDAFEGetNumberOfComponents(DM dm,PetscInt *d)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);
  if (d) *d = fe->ncomponents;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAFEGetElementSize"
PetscErrorCode DMDAFEGetElementSize(DM dm,PetscInt *ne,PetscInt r[])
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);
  if (ne) *ne = fe->nelements;
  if (r) {
    if (dm->dim == 2) {
      r[0] = fe->mx;
      r[1] = fe->my;
    }
    if (dm->dim == 3) {
      r[0] = fe->mx;
      r[1] = fe->my;
      r[2] = fe->mz;
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAFEGetLocalElementSize"
PetscErrorCode DMDAFEGetLocalElementSize(DM dm,PetscInt *ne,PetscInt r[])
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);
  if (ne) *ne = fe->lnelements;
  if (r) {
    if (dm->dim == 2) {
      r[0] = fe->lmx;
      r[1] = fe->lmy;
    }
    if (dm->dim == 3) {
      r[0] = fe->lmx;
      r[1] = fe->lmy;
      r[2] = fe->lmz;
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAFEGetElements"
PetscErrorCode DMDAFEGetElements(DM dm,PetscInt *ne,PetscInt *np,const PetscInt **map)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);
  if (ne) *ne = fe->lnelements;
  if (np) *np = fe->nodes_per_el;
  if (map) *map = fe->element_node_map;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAFEGetBasis"
PetscErrorCode DMDAFEGetBasis(DM dm,PetscInt *ne,PetscInt *np,const PetscInt **map)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);
  if (ne) *ne = fe->lnelements;
  if (np) *np = fe->nbasis;
  if (map) *map = fe->element_basis_map;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAFEGetDA"
PetscErrorCode DMDAFEGetDA(DM dm,DM *da)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);
  if (da) *da = fe->da;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAFEGetParentDM"
PetscErrorCode DMDAFEGetParentDM(DM dm,DM *parent)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);
  if (parent) *parent = fe->parent_dafe;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAFEGetOps"
PetscErrorCode DMDAFEGetOps(DM dm,DAFEOps *ops)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);
  if (ops) *ops = fe->ops;
  PetscFunctionReturn(0);
}

/* helper constructors */
#undef __FUNCT__
#define __FUNCT__ "DMDAFECreateOverlappingDMDAFE"
PetscErrorCode DMDAFECreateOverlappingDMDAFE(DM dm,PetscInt ncomponents,DAFEType type,DM *_dmsub)
{
  DM             dmsub;
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  DM_DAFE        *fesub;
  PetscErrorCode ierr;
  
  
  ierr = DAFECheckType(dm);CHKERRQ(ierr);

  /* Create structures */
  ierr = DMCreate(PetscObjectComm((PetscObject)dm),&dmsub);CHKERRQ(ierr);
  ierr = DMSetType(dmsub,"dafe");CHKERRQ(ierr);
  ierr = DMSetDimension(dmsub,dm->dim);CHKERRQ(ierr);
  ierr = DMDAFESetType(dmsub,type);CHKERRQ(ierr);
  ierr = DMDAFESetNumberOfComponents(dmsub,ncomponents);CHKERRQ(ierr);
  
  fesub = (DM_DAFE*)dmsub->data;
  fesub->parent_dafe = dm;
  ierr = PetscObjectReference((PetscObject)dm);CHKERRQ(ierr);
  
  switch (fe->dafetype) {

    case DAFE_Q2:
      if (type == DAFE_Q2sub) {
        SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Cannot create overlapped Q2 DAFE from parent of type Q2");
      } else if (type == DAFE_Q1sub) {
        //ierr = DAFECreateQ1FromQ2();CHKERRQ(ierr);
      } else if (type == DAFE_P1sub) {
        //ierr = DAFECreateP1FromQ2();CHKERRQ(ierr);
      } else if (type == DAFE_P0sub) {
        //ierr = DAFECreateP0FromQ2();CHKERRQ(ierr);
      }
      break;
      
    case DAFE_Q1:
      if (type == DAFE_Q2sub) {
        SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Cannot create overlapped Q2 DAFE from parent of type Q1");
      } else if (type == DAFE_Q1sub) {
        //ierr = DAFECreateQ1FromQ1();CHKERRQ(ierr);
      } else if (type == DAFE_P1sub) {
        //ierr = DAFECreateP1FromQ1();CHKERRQ(ierr);
      } else if (type == DAFE_P0sub) {
        //ierr = DAFECreateP0FromQ1();CHKERRQ(ierr);
      }
      break;
      
    case DAFE_Q2sub:
      SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Cannot create overlapping DAFE from parent of type Q2sub");
      break;
    case DAFE_Q1sub:
      SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Cannot create overlapping DAFE from parent of type Q2sub");
      break;
    case DAFE_P1sub:
      SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Cannot create overlapping DAFE from parent of type Q2sub");
      break;
    case DAFE_P0sub:
      SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Cannot create overlapping DAFE from parent of type Q2sub");
      break;
  }
  
  *_dmsub = dmsub;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAFECreateSubDMDAFE"
PetscErrorCode DMDAFECreateSubDMDAFE(DM dm,PetscInt ncomponents,PetscInt refine_factor,DAFEType type,DM *_dmsub)
{
  DM_DAFE        *fe = (DM_DAFE*)dm->data;
  DM             dmsub;
  DM_DAFE        *fesub;
  PetscErrorCode ierr;
  ierr = DAFECheckType(dm);CHKERRQ(ierr);

  ierr = DMCreate(PetscObjectComm((PetscObject)dm),&dmsub);CHKERRQ(ierr);
  ierr = DMSetType(dmsub,"dafe");CHKERRQ(ierr);
  ierr = DMSetDimension(dmsub,dm->dim);CHKERRQ(ierr);
  ierr = DMDAFESetType(dmsub,type);CHKERRQ(ierr);
  ierr = DMDAFESetNumberOfComponents(dmsub,ncomponents);CHKERRQ(ierr);
  
  fesub = (DM_DAFE*)dmsub->data;
  fesub->parent_dafe = dm;
  ierr = PetscObjectReference((PetscObject)dm);CHKERRQ(ierr);

  fesub->mx = fe->mx * refine_factor;
  fesub->my = fe->my * refine_factor;
  fesub->mz = fe->mz * refine_factor;

  fesub->nelements = fesub->mx * fesub->my * fesub->mz;
  
  fesub->lmx = fe->lmx * refine_factor;
  fesub->lmy = fe->lmy * refine_factor;
  fesub->lmz = fe->lmz * refine_factor;

  fesub->lnelements = fesub->lmx * fesub->lmy * fesub->lmz;
  
  
  switch (fe->dafetype) {

    case DAFE_Q2:
      if (type == DAFE_Q2sub) {
        SETERRQ(PetscObjectComm((PetscObject)dmsub),PETSC_ERR_SUP,"Cannot create refined Q2 DAFE from parent of type Q2");
      } else if (type == DAFE_Q1sub) {
        ierr = DAFE_CreateQ1FromQ2_3d(dmsub,refine_factor);CHKERRQ(ierr);
      } else if (type == DAFE_P1sub) {
        if (ncomponents > 1) SETERRQ(PetscObjectComm((PetscObject)dmsub),PETSC_ERR_SUP,"nDOF must be equal to 1");
        ierr = DAFE_CreateP1FromQ2_3d(dmsub,refine_factor);CHKERRQ(ierr);
      } else if (type == DAFE_P0sub) {
        if (ncomponents > 1) SETERRQ(PetscObjectComm((PetscObject)dmsub),PETSC_ERR_SUP,"nDOF must be equal to 1");
        ierr = DAFE_CreateP0FromQ2_3d(dmsub,refine_factor);CHKERRQ(ierr);
      }
      break;
      
    case DAFE_Q1:
      if (type == DAFE_Q2sub) {
        SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Cannot create refined Q2 DAFE from parent of type Q1");
      } else if (type == DAFE_Q1sub) {
        SETERRQ(PetscObjectComm((PetscObject)dmsub),PETSC_ERR_SUP,"TODO Q1 DAFE from parent of type Q1");
        
      } else if (type == DAFE_P1sub) {
        if (ncomponents > 1) SETERRQ(PetscObjectComm((PetscObject)dmsub),PETSC_ERR_SUP,"nDOF must be equal to 1");
        ierr = DAFE_CreateP1FromQ2_3d(dmsub,refine_factor);CHKERRQ(ierr);
      } else if (type == DAFE_P0sub) {
        if (ncomponents > 1) SETERRQ(PetscObjectComm((PetscObject)dmsub),PETSC_ERR_SUP,"nDOF must be equal to 1");
        ierr = DAFE_CreateP0FromQ2_3d(dmsub,refine_factor);CHKERRQ(ierr);
      }
      break;
    
    case DAFE_Q2sub:
      SETERRQ(PetscObjectComm((PetscObject)dmsub),PETSC_ERR_SUP,"Cannot create refined DAFE from parent of type Q2sub");
      break;
    case DAFE_Q1sub:
      SETERRQ(PetscObjectComm((PetscObject)dmsub),PETSC_ERR_SUP,"Cannot create refined DAFE from parent of type Q1sub");
      break;
    case DAFE_P1sub:
      SETERRQ(PetscObjectComm((PetscObject)dmsub),PETSC_ERR_SUP,"Cannot create refined DAFE from parent of type P1sub");
      break;
    case DAFE_P0sub:
      SETERRQ(PetscObjectComm((PetscObject)dmsub),PETSC_ERR_SUP,"Cannot create refined DAFE from parent of type P0sub");
      break;
  }
  *_dmsub = dmsub;
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DMDAFECreate3d"
PetscErrorCode DMDAFECreate3d(MPI_Comm comm,PetscInt ncomponents,DAFEType dafetype,PetscInt mi,PetscInt mj,PetscInt mk,DM *_dm)
{
  PetscErrorCode ierr;
  DM             dm;
  DM_DAFE        *fe;
  
  ierr = DMCreate(comm,&dm);CHKERRQ(ierr);
  ierr = DMSetType(dm,"dafe");CHKERRQ(ierr);
  ierr = DMSetDimension(dm,3);CHKERRQ(ierr);
  ierr = DMDAFESetType(dm,dafetype);CHKERRQ(ierr);
  fe = (DM_DAFE*)dm->data;
  ierr = DMDAFESetNumberOfComponents(dm,ncomponents);CHKERRQ(ierr);

  fe->mx = mi;
  fe->my = mj;
  fe->mz = mk;
  fe->nelements = mi * mj * mk;
  if (dafetype == DAFE_Q2) {
    ierr = DAFE_CreateQ2_3d(dm,mi,mj,mk);CHKERRQ(ierr);
  } else if (dafetype == DAFE_Q1) {
    SETERRQ(comm,PETSC_ERR_SUP,"TODO Constructor support for DAFE_Q1");
  } else {
    SETERRQ(comm,PETSC_ERR_SUP,"Constructor only supports DAFE_Q2, DAFE_Q1");
  }
  *_dm = dm;
  PetscFunctionReturn(0);
}

/*
 Interpolate coordinates FROM the parent to child DAFE
 */
#undef __FUNCT__
#define __FUNCT__ "DMDAFEProjectCoordinates"
PetscErrorCode DMDAFEProjectCoordinates(DM sub)
{
  DM             parent = NULL;
  DM_DAFE        *fe = (DM_DAFE*)sub->data;
  PetscErrorCode ierr;
  
  if (!fe->ops->ProjectCoordinates) SETERRQ1(PetscObjectComm((PetscObject)sub),PETSC_ERR_SUP,"No method to project coordinates for DAFE type %s",DAFETypeNames[(PetscInt)fe->dafetype]);
  parent = fe->parent_dafe;
  if (!parent) SETERRQ(PetscObjectComm((PetscObject)sub),PETSC_ERR_SUP,"DAFE does not provide a parent DM");
  
  ierr = fe->ops->ProjectCoordinates(parent,sub);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/*
 Inject coordinates FROM the child DAFE to parent
 */
#undef __FUNCT__
#define __FUNCT__ "DMDAFEInjectCoordinates"
PetscErrorCode DMDAFEInjectCoordinates(DM dmchild)
{
  DM_DAFE        *fe = (DM_DAFE*)dmchild->data;
  DM_DAFE        *parent_fe;
  DM             parent;
  
  if (!fe->parent_dafe) SETERRQ(PetscObjectComm((PetscObject)dmchild),PETSC_ERR_USER,"Missing parent DAFE. Requires DAFE was created from another DAFE");
  parent = fe->parent_dafe;
  parent_fe = (DM_DAFE*)parent->data;

  SETERRQ(PetscObjectComm((PetscObject)dmchild),PETSC_ERR_SUP,"TODO ");
  PetscFunctionReturn(0);
}

