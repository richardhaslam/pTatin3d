
#ifndef __xdmf_writer_h__
#define __xdmf_writer_h__

typedef enum { XDMFNode=0 , XDMFEdge , XDMFFace , XDMFCell , XDMFGrid } XDMFCenter;

typedef enum { XDMFScalar=0 , XDMFVector , XDMFTensor , XDMFTensor6, XDMFMatrix , XDMFMultiCompVector } XDMFAttribute;

typedef enum { XDMFXML=0 , XDMFHDF5 , XDMFBinary } XDMFDataItemFormat;

typedef enum { XDMFNative=0 , XDMFBigEndian , XDMFLittleEndian  } XDMFDataItemEndian;

extern const char *XDMFCenterNames[];
extern const char *XDMFAttributeNames[];
extern const char *XDMFDataItemFormatNames[];
extern const char *XDMFDataItemEndianNames[];


PetscErrorCode XDMFMetaWriteTime(PetscViewer v,PetscReal value);
PetscErrorCode XDMFMetaWriteInformationString(PetscViewer v,const char name[],const char value[]);
PetscErrorCode XDMFMetaWriteInformationInt(PetscViewer v,const char name[],PetscInt value);
PetscErrorCode XDMFMetaWriteInformationReal(PetscViewer v,const char name[],PetscReal value);
PetscErrorCode XDMFMetaWriteInformationRealList(PetscViewer v,const char name[],PetscInt n,PetscReal values[]);
PetscErrorCode XDMFMetaXDMFOpen(MPI_Comm comm,const char name[],PetscViewer *v);
PetscErrorCode XDMFMetaXDMFClose(PetscViewer *v);
PetscErrorCode XDMFMetaDomainOpen(PetscViewer v,const char name[]);
PetscErrorCode XDMFMetaDomainClose(PetscViewer v);
PetscErrorCode XDMFGridOpen_DMDA(PetscViewer v,DM dm,const char path[],const char suffix[],const char meshname[],XDMFDataItemFormat format);
PetscErrorCode XDMFGridClose_DMDA(PetscViewer v);
PetscErrorCode XDMFWriteAttribute_DMDA(PetscViewer v,DM da,Vec x,
                                       const char path[],const char suffix[],const char meshname[],const char fieldname[],
                                       XDMFCenter c_type,XDMFDataItemFormat format);
PetscErrorCode XDMFDataWriteField_Generic(Vec x,
                                          const char path[],const char filename[],
                                          XDMFDataItemFormat format);
PetscErrorCode XDMFDataWriteField_GenericDMDA(DM dm,Vec x,
                                              const char path[],const char filename[],
                                              XDMFDataItemFormat format);
PetscErrorCode ptatin3d_StokesOutput_VelocityXDMF(pTatinCtx ctx,Vec X,const char suffix[]);

#endif
