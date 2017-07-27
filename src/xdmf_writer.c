

#include "petscviewer.h"

#include "ptatin3d.h"
#include "ptatin3d_defs.h"
#include "private/ptatin_impl.h"
#include "ptatin_version_info.h"

#include "xdmf_writer.h"

const char *XDMFCenterNames[] = { "Node" , "Edge" , "Face" , "Cell" , "Grid" ,  0  };
const char *XDMFAttributeNames[] = { "Scalar" , "Vector" , "Tensor" ,  "Tensor6" , "Matrix" , "MultiComponent" , 0  };
const char *XDMFDataItemFormatNames[] = { "XML" , "HDF" , "Binary" , 0 };
const char *XDMFDataItemEndianNames[] = { "Native" , "Big" , "Little" , 0 };





#undef __FUNCT__
#define __FUNCT__ "_XDMFMeta_XDMFOpenClose"
PetscErrorCode _XDMFMeta_XDMFOpenClose(MPI_Comm comm,const char name[],PetscBool open,PetscViewer *v)
{
    PetscErrorCode ierr;

    if (open) {
        ierr = PetscViewerASCIIOpen(comm,name,v);CHKERRQ(ierr);

        PetscViewerASCIIPrintf(*v,"<?xml version=\"1.0\" ?>\n");
        PetscViewerASCIIPrintf(*v,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        PetscViewerASCIIPrintf(*v,"<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n");
    } else {
        PetscViewerASCIIPrintf(*v,"</Xdmf>\n");
        ierr = PetscViewerDestroy(v);CHKERRQ(ierr);
        *v = NULL;
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_XDMFMeta_DomainOpenClose"
PetscErrorCode _XDMFMeta_DomainOpenClose(PetscViewer v,const char name[],PetscBool open)
{
    if (open) {
        if (name) {
            PetscViewerASCIIPrintf(v,"<Domain Name=\"%s\">\n",name);
        } else {
            PetscViewerASCIIPrintf(v,"<Domain>\n");
        }
    }
    else { PetscViewerASCIIPrintf(v,"</Domain>\n"); }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_XDMFMeta_GridOpenClose_DMDA"
PetscErrorCode _XDMFMeta_GridOpenClose_DMDA(PetscViewer v,DM da,const char suffix[],const char meshname[],XDMFDataItemFormat format,PetscBool open)
{
    PetscErrorCode ierr;
    
    if (open) {
        PetscInt M,N,P;

        if (!meshname) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"XDMF_AddDMDA requires meshname"); }

        ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);

        PetscViewerASCIIPrintf(v,"<Grid Name=\"%s\" GridType=\"Uniform\">\n",meshname);
        ierr = PetscViewerASCIIPushTab(v);CHKERRQ(ierr);

        PetscViewerASCIIPrintf(v,"<Topology TopologyType=\"3DSMesh\" Dimensions=\"%D %D %D\"/>\n",P,N,M);

        PetscViewerASCIIPrintf(v,"<Geometry GeometryType=\"XYZ\">\n");
        ierr = PetscViewerASCIIPushTab(v);CHKERRQ(ierr);

        PetscViewerASCIIPrintf(v,"<DataItem Dimensions=\"%D %D %D 3\"\n",P,N,M);
        PetscViewerASCIIPrintf(v," NumberType=\"Float\" Precision=\"8\"\n");

        switch (format) {
            case XDMFXML:
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"XML (ascii) data items are not implemented");
                break;

            case XDMFHDF5:
            {
                const char *vecname;
                Vec        coords;

                ierr = DMGetCoordinates(da,&coords);CHKERRQ(ierr);
                ierr = PetscObjectGetName((PetscObject)coords,&vecname);CHKERRQ(ierr);
                PetscViewerASCIIPrintf(v," Format=\"%s\">\n",XDMFDataItemFormatNames[XDMFHDF5]);
                if (suffix) { PetscViewerASCIIPrintf(v,"   %s_%s_coor.h5:/%s\n",suffix,meshname,vecname); }
                else {        PetscViewerASCIIPrintf(v,"   %s_coor.h5:/%s\n",meshname,vecname); }
                break;
            }

            case XDMFBinary:
                PetscViewerASCIIPrintf(v," Format=\"%s\" Endian=\"%s\">\n",XDMFDataItemFormatNames[XDMFBinary],XDMFDataItemEndianNames[XDMFBigEndian]);
                if (suffix) { PetscViewerASCIIPrintf(v,"   %s_%s_coor.pbvec\n",suffix,meshname); }
                else {        PetscViewerASCIIPrintf(v,"   %s_coor.pbvec\n",meshname); }
                break;

            default:
                PetscViewerASCIIPrintf(v," Format=\"%s\" Endian=\"%s\">\n",XDMFDataItemFormatNames[XDMFBinary],XDMFDataItemEndianNames[XDMFBigEndian]);
                if (suffix) { PetscViewerASCIIPrintf(v,"   %s_%s_coor.pbvec\n",suffix,meshname); }
                else {        PetscViewerASCIIPrintf(v,"   %s_coor.pbvec\n",meshname); }
                break;
        }

        PetscViewerASCIIPrintf(v,"</DataItem>\n");
        ierr = PetscViewerASCIIPopTab(v);CHKERRQ(ierr);
        PetscViewerASCIIPrintf(v,"</Geometry>\n");

    } else {
        ierr = PetscViewerASCIIPopTab(v);CHKERRQ(ierr);
        PetscViewerASCIIPrintf(v,"</Grid>\n");
    }

    PetscFunctionReturn(0);
}

/*
 
 Notes:
 
 [1] A PETSc Vec allows only contains PetscScalar data.
 Thus, when writing a vector into an xdmf file, the data item is assumed to be
 NumberType=Float Precision=8 <default>
 NumberType=Float Precision=4 <if petsc compiled with single precision>
 
 [2] PetscBinaryWrite performs byte swapping and writes all data as big-endian.
 Hence, when writing binary we must specify that the data item uses
 format=Binary Endian=Big
 
*/
#undef __FUNCT__
#define __FUNCT__ "_XDMFMeta_AddAttributeField_DMDA"
PetscErrorCode _XDMFMeta_AddAttributeField_DMDA(PetscViewer v,DM da,Vec x,
                                               const char suffix[],const char meshname[],const char fieldname[],
                                               XDMFCenter c_type,XDMFDataItemFormat format)
{
    XDMFAttribute     attr_type;
    PetscInt          d,ndof,M,N,P;
    PetscErrorCode    ierr;

    if (!meshname) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"XDMF_AddField requires meshname"); }
    if (!fieldname) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"XDMF_AddField requires fieldname"); }

    ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,&ndof,0,0,0,0,0);CHKERRQ(ierr);
    switch (ndof) {
        case 1:
            attr_type = XDMFScalar;
            break;
        case 3:
            attr_type = XDMFVector;
            break;
        case 6:
            attr_type = XDMFTensor6;
            break;
        case 9:
            attr_type = XDMFTensor;
            break;
        default:
            attr_type = XDMFMultiCompVector;
            break;
    }

    switch (attr_type) {
        case XDMFScalar:
            PetscViewerASCIIPrintf(v,"<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"%s\">\n",fieldname,XDMFAttributeNames[attr_type],XDMFCenterNames[c_type]);
            ierr = PetscViewerASCIIPushTab(v);CHKERRQ(ierr);

            PetscViewerASCIIPrintf(v,"<DataItem Dimensions=\"%D %D %D 1\"\n",P,N,M);
            break;
        case XDMFVector:
            PetscViewerASCIIPrintf(v,"<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"%s\">\n",fieldname,XDMFAttributeNames[attr_type],XDMFCenterNames[c_type]);
            ierr = PetscViewerASCIIPushTab(v);CHKERRQ(ierr);

            PetscViewerASCIIPrintf(v,"<DataItem Dimensions=\"%D %D %D 3\"\n",P,N,M);
            break;
        case XDMFTensor6:
            PetscViewerASCIIPrintf(v,"<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"%s\">\n",fieldname,XDMFAttributeNames[attr_type],XDMFCenterNames[c_type]);
            ierr = PetscViewerASCIIPushTab(v);CHKERRQ(ierr);

            PetscViewerASCIIPrintf(v,"<DataItem Dimensions=\"%D %D %D 6\"\n",P,N,M);
            break;
        case XDMFTensor:
            PetscViewerASCIIPrintf(v,"<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"%s\">\n",fieldname,XDMFAttributeNames[attr_type],XDMFCenterNames[c_type]);
            ierr = PetscViewerASCIIPushTab(v);CHKERRQ(ierr);

            PetscViewerASCIIPrintf(v,"<DataItem Dimensions=\"%D %D %D 9\"\n",P,N,M);
            break;
        case XDMFMultiCompVector:
            break;
        default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Unknown XDMF AttributeType specified");
            break;
    }

    if (attr_type != XDMFMultiCompVector) {
        PetscViewerASCIIPrintf(v," NumberType=\"Float\" Precision=\"8\"\n");

        if (format == XDMFHDF5) {
            const char *vecname;

            ierr = PetscObjectGetName((PetscObject)x,&vecname);CHKERRQ(ierr);
            PetscViewerASCIIPrintf(v," Format=\"HDF\">\n");
            //if (suffix) { PetscViewerASCIIPrintf(v,"   %s_%s_fields.h5:/%s\n",suffix,meshname,fieldname); }
            //else {        PetscViewerASCIIPrintf(v,"   %s_fields.h5:/%s\n",meshname,fieldname); }
            if (suffix) { PetscViewerASCIIPrintf(v,"   %s_%s_%s.h5:/%s\n",suffix,meshname,fieldname,vecname); }
            else {        PetscViewerASCIIPrintf(v,"   %s_%s.h5:/%s\n",meshname,fieldname,vecname); }
        } else {
            PetscViewerASCIIPrintf(v," Format=\"Binary\" Endian=\"Big\">\n");
            if (suffix) { PetscViewerASCIIPrintf(v,"   %s_%s_%s.pbvec\n",suffix,meshname,fieldname); }
            else {        PetscViewerASCIIPrintf(v,"   %s_%s.pbvec\n",meshname,fieldname); }
        }
        PetscViewerASCIIPrintf(v,"</DataItem>\n");
        ierr = PetscViewerASCIIPopTab(v);CHKERRQ(ierr);

        PetscViewerASCIIPrintf(v,"</Attribute>\n");
    }

    if (attr_type == XDMFMultiCompVector) {
        for (d=0; d<ndof; d++) {
            char fieldname_comp[PETSC_MAX_PATH_LEN];

            PetscSNPrintf(fieldname_comp,PETSC_MAX_PATH_LEN-1,"%s_%D",fieldname,d);
            PetscViewerASCIIPrintf(v,"<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"%s\">\n",fieldname_comp,XDMFAttributeNames[XDMFScalar],XDMFCenterNames[c_type]);
            ierr = PetscViewerASCIIPushTab(v);CHKERRQ(ierr);

            PetscViewerASCIIPrintf(v,"<DataItem Dimensions=\"%D %D %D 1\"\n",P,N,M);
            PetscViewerASCIIPrintf(v," NumberType=\"Float\" Precision=\"8\"\n");

            if (format == XDMFHDF5) {
                const char *vecname;

                ierr = PetscObjectGetName((PetscObject)x,&vecname);CHKERRQ(ierr);
                PetscViewerASCIIPrintf(v," Format=\"HDF\">\n");
                if (suffix) { PetscViewerASCIIPrintf(v,"   %s_%s_%s.h5:/%s_%D\n",suffix,meshname,fieldname,vecname,d); }
                else {        PetscViewerASCIIPrintf(v,"   %s_%s.h5:/%s_%D\n",meshname,fieldname,vecname,d); }
            } else {
                PetscViewerASCIIPrintf(v," Format=\"Binary\" Endian=\"Big\">\n");
                if (suffix) { PetscViewerASCIIPrintf(v,"   %s_%s_%s_%D.pbvec\n",suffix,meshname,fieldname,d); }
                else {        PetscViewerASCIIPrintf(v,"   %s_%s_%D.pbvec\n",meshname,fieldname,d); }
            }
            PetscViewerASCIIPrintf(v,"</DataItem>\n");
            ierr = PetscViewerASCIIPopTab(v);CHKERRQ(ierr);

            PetscViewerASCIIPrintf(v,"</Attribute>\n");
        }
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFDataWriteField_Generic"
PetscErrorCode XDMFDataWriteField_Generic(Vec x,
                                           const char path[],const char filename[],
                                           XDMFDataItemFormat format)
{
    MPI_Comm       comm;
    PetscViewer    viewer;
	char           name[PETSC_MAX_PATH_LEN];
    PetscErrorCode ierr;

    ierr = PetscObjectGetComm((PetscObject)x,&comm);CHKERRQ(ierr);

    switch (format) {

        case XDMFBinary:
        {
            PetscBool ismpiio;

            ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
            ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
            ierr = PetscViewerBinarySetSkipHeader(viewer,PETSC_TRUE);CHKERRQ(ierr);
            ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
            if (path) { PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s/%s.pbvec",path,filename); }
            else {      PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s.pbvec",filename); }
            ierr = PetscViewerFileSetName(viewer,name);CHKERRQ(ierr);

            ierr = PetscViewerBinaryGetUseMPIIO(viewer,&ismpiio);CHKERRQ(ierr);
            if (ismpiio) {
                PetscPrintf(comm,"*** XDMFDataWriteField_Generic using MPI-IO ***\n");
            }

            ierr = VecView(x,viewer);CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
            break;
        }

        case XDMFHDF5:
        {
            ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
            ierr = PetscViewerSetType(viewer,PETSCVIEWERHDF5);CHKERRQ(ierr);
            ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
            if (path) { PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s/%s.h5",path,filename); }
            else {      PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s.h5",filename); }
            ierr = PetscViewerFileSetName(viewer,name);CHKERRQ(ierr);
            ierr = VecView(x,viewer);CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
            break;
        }

        case XDMFXML:
          SETERRQ(comm,PETSC_ERR_SUP,"Format not supported");
    }

    PetscFunctionReturn(0);
}

/*
 As of PETSc 3.5.4, the header-less PetscViewerBinary using MPI-IO
 appears to be have a bug which cause VecView() to fail on vectors
 created from 3D DMDA's.
 Writing standard vectors using MPI-IO appears to function correctly.

 The source of this bug is unknown to me at this time.
 As a work around, I will write a flat vector in the natural
 ordering using VecView() of a standard (non-DM created) vector.
*/
#undef __FUNCT__
#define __FUNCT__ "XDMFDataWriteField_GenericDMDA"
PetscErrorCode XDMFDataWriteField_GenericDMDA(DM dm,Vec x,
                                          const char path[],const char filename[],
                                          XDMFDataItemFormat format)
{
    MPI_Comm       comm;
    PetscViewer    viewer;
    char           name[PETSC_MAX_PATH_LEN];
    Vec            xn;
    PetscErrorCode ierr;

	ierr = DMDACreateNaturalVector(dm,&xn);CHKERRQ(ierr);
    ierr = PetscObjectGetComm((PetscObject)xn,&comm);CHKERRQ(ierr);
	ierr = DMDAGlobalToNaturalBegin(dm,x,INSERT_VALUES,xn);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalEnd(dm,x,INSERT_VALUES,xn);CHKERRQ(ierr);

    switch (format) {

        case XDMFBinary:
        {
            PetscBool ismpiio;

            ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
            ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
            ierr = PetscViewerBinarySetSkipHeader(viewer,PETSC_TRUE);CHKERRQ(ierr);
            ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
            if (path) { PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s/%s.pbvec",path,filename); }
            else {      PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s.pbvec",filename); }
            ierr = PetscViewerFileSetName(viewer,name);CHKERRQ(ierr);
            ierr = PetscViewerBinaryGetUseMPIIO(viewer,&ismpiio);CHKERRQ(ierr);
            if (ismpiio) {
                PetscPrintf(comm,"*** XDMFDataWriteField_GenericDMDA using MPI-IO ***\n");
            }
            ierr = VecView(xn,viewer);CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
            break;
        }

        case XDMFHDF5:
        {
            ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
            ierr = PetscViewerSetType(viewer,PETSCVIEWERHDF5);CHKERRQ(ierr);
            ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
            if (path) { PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s/%s.h5",path,filename); }
            else {      PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s.h5",filename); }
            ierr = PetscViewerFileSetName(viewer,name);CHKERRQ(ierr);
            ierr = VecView(xn,viewer);CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
            break;
        }
        case XDMFXML:
          SETERRQ(comm,PETSC_ERR_SUP,"Format not supported");
    }
    ierr = VecDestroy(&xn);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFWriteAttribute_DMDA"
PetscErrorCode XDMFWriteAttribute_DMDA(PetscViewer v,DM da,Vec x,
                               const char path[],const char suffix[],const char meshname[],const char fieldname[],
                               XDMFCenter c_type,XDMFDataItemFormat format)
{
	char           filename[PETSC_MAX_PATH_LEN];
    PetscErrorCode ierr;

    if (!v) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"XDMF file pointer is corrupt"); }

    ierr = _XDMFMeta_AddAttributeField_DMDA(v,da,x,suffix,meshname,fieldname,c_type,format);CHKERRQ(ierr);

    if (suffix) {
        PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s_%s_%s",suffix,meshname,fieldname);
    } else {
        PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s_%s",meshname,fieldname);
    }
    ierr = XDMFDataWriteField_GenericDMDA(da,x,path,filename,format);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFWriteDataItemByReference_DMDA"
PetscErrorCode XDMFWriteDataItemByReference_DMDA(PetscViewer v,DM dm,const char fieldname[],const char filename[],XDMFDataItemFormat format)
{
    PetscInt       ndof,M,N,P;
    PetscErrorCode ierr;

    ierr = DMDAGetInfo(dm,0,&M,&N,&P,0,0,0,&ndof,0,0,0,0,0);CHKERRQ(ierr);

    PetscViewerASCIIPrintf(v,"<DataItem Name=\"ref_%s\" Dimensions=\"%D %D %D %D\"\n",fieldname,P,N,M,ndof);
    PetscViewerASCIIPrintf(v,"  NumberType=\"Float\" Precision=\"8\"\n");

    if (format == XDMFBinary) {
        PetscViewerASCIIPrintf(v,"  Format=\"%s\" Endian=\"%s\">\n",XDMFDataItemFormatNames[format],XDMFDataItemEndianNames[XDMFBigEndian]);
    } else {
        PetscViewerASCIIPrintf(v,"  Format=\"%s\" Endian=\"%s\">\n",XDMFDataItemFormatNames[format]);
    }
    PetscViewerASCIIPrintf(v,"  %s\n",filename);
    PetscViewerASCIIPrintf(v,"</DataItem>\n");

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "_XDMFMeta_AddAttributeFunctionField_DMDA"
PetscErrorCode _XDMFMeta_AddAttributeFunctionField_DMDA(PetscViewer v,const char function[],
                                                DM da,Vec x,
                                                const char suffix[],const char meshname[],const char fieldname[],
                                                XDMFCenter c_type,XDMFDataItemFormat format)
{
    XDMFAttribute     attr_type;
    PetscInt          ndof,M,N,P;
    char              reference_key[PETSC_MAX_PATH_LEN];
    PetscErrorCode    ierr;

    if (!meshname) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"XDMF_AddField requires meshname"); }
    if (!fieldname) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"XDMF_AddField requires fieldname"); }

    ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,&ndof,0,0,0,0,0);CHKERRQ(ierr);
    switch (ndof) {
        case 1:
            attr_type = XDMFScalar;
            break;
        case 3:
            attr_type = XDMFVector;
            break;
        case 6:
            attr_type = XDMFTensor6;
            break;
        case 9:
            attr_type = XDMFTensor;
            break;
        default:
            attr_type = XDMFMultiCompVector;
            break;
    }

    PetscViewerASCIIPrintf(v,"<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"%s\">\n",fieldname,XDMFAttributeNames[attr_type],XDMFCenterNames[c_type]);
    ierr = PetscViewerASCIIPushTab(v);CHKERRQ(ierr);

    PetscViewerASCIIPrintf(v,"<DataItem ItemType=\"Function\" Function=\"%s\" Dimensions=\"%D %D %D %D\">\n",function,P,N,M,ndof);
    ierr = PetscViewerASCIIPushTab(v);CHKERRQ(ierr);

    PetscViewerASCIIPrintf(v,"<DataItem Reference=\"XML\">\n");

    PetscSNPrintf(reference_key,PETSC_MAX_PATH_LEN-1,"/Xdmf/Domain/Grid[@Name=\"%s\"]/DataItem[@Name=\"ref_%s\"]",
                  meshname,fieldname);

    PetscViewerASCIIPrintf(v,"  %s\n",reference_key);

    PetscViewerASCIIPrintf(v,"</DataItem>\n");
    ierr = PetscViewerASCIIPopTab(v);CHKERRQ(ierr);

    PetscViewerASCIIPrintf(v,"</DataItem>\n");
    ierr = PetscViewerASCIIPopTab(v);CHKERRQ(ierr);

    PetscViewerASCIIPrintf(v,"</Attribute>\n");

    if (format == XDMFHDF5) {
        const char *vecname;

        ierr = PetscObjectGetName((PetscObject)x,&vecname);CHKERRQ(ierr);
        if (suffix) { PetscSNPrintf(reference_key,PETSC_MAX_PATH_LEN-1,"%s_%s_%s.h5:/%s",suffix,meshname,fieldname,vecname); }
        else {        PetscSNPrintf(reference_key,PETSC_MAX_PATH_LEN-1,"%s_%s.h5:/%s",meshname,fieldname,vecname); }
    } else {
        if (suffix) { PetscSNPrintf(reference_key,PETSC_MAX_PATH_LEN-1,"%s_%s_%s.pbvec",suffix,meshname,fieldname); }
        else {        PetscSNPrintf(reference_key,PETSC_MAX_PATH_LEN-1,"%s_%s.pbvec",meshname,fieldname); }
    }
    ierr = XDMFWriteDataItemByReference_DMDA(v,da,fieldname,reference_key,format);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFWriteAttributeFunction_DMDA"
PetscErrorCode XDMFWriteAttributeFunction_DMDA(PetscViewer v,const char function[],
                                               DM da,Vec x,
                                               const char path[],const char suffix[],const char meshname[],const char fieldname[],
                                               XDMFCenter c_type,XDMFDataItemFormat format)
{
    PetscErrorCode ierr;

    if (!v) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"XDMF file pointer is corrupt"); }

    ierr = _XDMFMeta_AddAttributeFunctionField_DMDA(v,function,da,x,suffix,meshname,fieldname,c_type,format);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFMetaXDMFOpen"
PetscErrorCode XDMFMetaXDMFOpen(MPI_Comm comm,const char name[],PetscViewer *v)
{
    PetscErrorCode ierr;

    ierr = _XDMFMeta_XDMFOpenClose(comm,name,PETSC_TRUE,v);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(*v);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFMetaXDMFClose"
PetscErrorCode XDMFMetaXDMFClose(PetscViewer *v)
{
    PetscErrorCode ierr;

    ierr = PetscViewerASCIIPopTab(*v);CHKERRQ(ierr);
    ierr = _XDMFMeta_XDMFOpenClose(MPI_COMM_NULL,NULL,PETSC_FALSE,v);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFMetaDomainOpen"
PetscErrorCode XDMFMetaDomainOpen(PetscViewer v,const char name[])
{
    PetscErrorCode ierr;

    if (!v) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"XDMF file pointer is corrupt"); }
    ierr = _XDMFMeta_DomainOpenClose(v,name,PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(v);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFMetaDomainClose"
PetscErrorCode XDMFMetaDomainClose(PetscViewer v)
{
    PetscErrorCode ierr;

    if (!v) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"XDMF file pointer is corrupt"); }
    ierr = PetscViewerASCIIPopTab(v);CHKERRQ(ierr);
    ierr = _XDMFMeta_DomainOpenClose(v,NULL,PETSC_FALSE);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFMetaWriteTime"
PetscErrorCode XDMFMetaWriteTime(PetscViewer v,PetscReal value)
{
    if (!v) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"XDMF file pointer is corrupt"); }
    PetscViewerASCIIPrintf(v,"<Time Value=\"%1.6e\"/>\n",value);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFMetaWriteInformationString"
PetscErrorCode XDMFMetaWriteInformationString(PetscViewer v,const char name[],const char value[])
{
    if (!v) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"XDMF file pointer is corrupt"); }
    PetscViewerASCIIPrintf(v,"<Information Name=\"%s\" Value=\"%s\"/>\n",name,value);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFMetaWriteInformationInt"
PetscErrorCode XDMFMetaWriteInformationInt(PetscViewer v,const char name[],PetscInt value)
{
    if (!v) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"XDMF file pointer is corrupt"); }
    PetscViewerASCIIPrintf(v,"<Information Name=\"%s\" Value=\"%D\"/>\n",name,value);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFMetaWriteInformationReal"
PetscErrorCode XDMFMetaWriteInformationReal(PetscViewer v,const char name[],PetscReal value)
{
    if (!v) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"XDMF file pointer is corrupt"); }
    PetscViewerASCIIPrintf(v,"<Information Name=\"%s\" Value=\"%1.6e\"/>\n",name,value);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFMetaWriteInformationRealList"
PetscErrorCode XDMFMetaWriteInformationRealList(PetscViewer v,const char name[],PetscInt n,PetscReal values[])
{
    PetscInt       i;
    PetscErrorCode ierr;

    if (!v) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"XDMF file pointer is corrupt"); }
    PetscViewerASCIIPrintf(v,"<Information Name=\"%s\">\n",name);
    ierr = PetscViewerASCIIPushTab(v);CHKERRQ(ierr);
    for (i=0; i<n; i++) {
        PetscViewerASCIIPrintf(v,"%1.6e ",values[i]);
    }
    PetscViewerASCIIPrintf(v,"\n");
    ierr = PetscViewerASCIIPopTab(v);CHKERRQ(ierr);
    PetscViewerASCIIPrintf(v,"</Information>\n");

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFGridOpen_DMDA"
PetscErrorCode XDMFGridOpen_DMDA(PetscViewer v,DM dm,const char path[],const char suffix[],const char meshname[],XDMFDataItemFormat format)
{
    char           filename[PETSC_MAX_PATH_LEN];
    Vec            coords;
    DM             cdm;
    PetscErrorCode ierr;

    if (!v) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"XDMF file pointer is corrupt"); }
    ierr = _XDMFMeta_GridOpenClose_DMDA(v,dm,suffix,meshname,format,PETSC_TRUE);CHKERRQ(ierr);

    /* write out coordinate vector */
    ierr = DMGetCoordinateDM(dm,&cdm);CHKERRQ(ierr);
    ierr = DMGetCoordinates(dm,&coords);CHKERRQ(ierr);
    if (suffix) {
        PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s_%s_coor",suffix,meshname);
    } else {
        PetscSNPrintf(filename,PETSC_MAX_PATH_LEN-1,"%s_coor",meshname);
    }
    ierr = XDMFDataWriteField_GenericDMDA(cdm,coords,path,filename,format);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "XDMFGridClose_DMDA"
PetscErrorCode XDMFGridClose_DMDA(PetscViewer v)
{
    PetscErrorCode ierr;

    if (!v) { SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"XDMF file pointer is corrupt"); }
    ierr = _XDMFMeta_GridOpenClose_DMDA(v,NULL,NULL,NULL,XDMFBinary,PETSC_FALSE);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ptatin3d_StokesOutput_VelocityXDMF"
PetscErrorCode ptatin3d_StokesOutput_VelocityXDMF(pTatinCtx ctx,Vec X,const char suffix[])
{
    PhysCompStokes stokes;
    DM             dmstokes,dmv,dmp;
    Vec            velocity,pressure;
    MPI_Comm       comm;
    PetscViewer    viewer;
	char           name[PETSC_MAX_PATH_LEN];
    char           infostr[PETSC_MAX_PATH_LEN];
    char           *model_name;
    PetscBool      useH5 = PETSC_FALSE;
    PetscLogDouble t0,t1;
    pTatinModel    model;
    PetscErrorCode ierr;

	ierr = pTatinGetStokesContext(ctx,&stokes);CHKERRQ(ierr);
    ierr = PhysCompStokesGetDMComposite(stokes,&dmstokes);
    ierr = DMCompositeGetEntries(dmstokes,&dmv,&dmp);CHKERRQ(ierr);
    ierr = DMCompositeGetAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);

    PetscObjectGetComm((PetscObject)dmstokes,&comm);
    if (suffix) { PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s/%s_dmda_vel.xmf",ctx->outputpath,suffix); }
    else {        PetscSNPrintf(name,PETSC_MAX_PATH_LEN-1,"%s/dmda_vel.xmf",ctx->outputpath); }

    ierr = XDMFMetaXDMFOpen(comm,name,&viewer);CHKERRQ(ierr);

    /* log info */
    ierr = PetscGetVersion(infostr,PETSC_MAX_PATH_LEN-1);CHKERRQ(ierr);
    ierr = XDMFMetaWriteInformationString(viewer,"PetscVersionInfo",infostr);CHKERRQ(ierr);

    //ierr = XDMFMetaWriteInformationString(viewer,"pTatinVersionInfo",PTATIN_VERSION_CNTR_REVISION);CHKERRQ(ierr);

    ierr = PetscGetProgramName(infostr,PETSC_MAX_PATH_LEN-1);CHKERRQ(ierr);
    ierr = XDMFMetaWriteInformationString(viewer,"GeneratedBy",infostr);CHKERRQ(ierr);

    ierr = XDMFMetaWriteInformationString(viewer,"Timestamp",ctx->formatted_timestamp);CHKERRQ(ierr);

    ierr = pTatinGetModel(ctx,&model);CHKERRQ(ierr);
    ierr = pTatinModelGetName(model,&model_name);CHKERRQ(ierr);
    ierr = XDMFMetaWriteInformationString(viewer,"ModelName",model_name);CHKERRQ(ierr);

    /* open domain and add dmda */
    ierr = XDMFMetaDomainOpen(viewer,"pTatin");CHKERRQ(ierr);
    ierr = XDMFMetaWriteInformationReal(viewer,"timestep",ctx->dt);CHKERRQ(ierr);
    ierr = XDMFMetaWriteInformationInt(viewer,"step",ctx->step);CHKERRQ(ierr);

    PetscOptionsGetBool(NULL,NULL,"-xdmf_use_hdf",&useH5,NULL);
    PetscTime(&t0);
    if (!useH5) {
        ierr = XDMFGridOpen_DMDA(viewer,dmv,ctx->outputpath,suffix,"stokes",XDMFBinary);CHKERRQ(ierr);
        ierr = XDMFWriteAttribute_DMDA(viewer,dmv,velocity,ctx->outputpath,suffix,"stokes","velocity",XDMFNode,XDMFBinary);CHKERRQ(ierr);
    } else {
        ierr = XDMFGridOpen_DMDA(viewer,dmv,ctx->outputpath,suffix,"stokes",XDMFHDF5);CHKERRQ(ierr);
        ierr = XDMFWriteAttribute_DMDA(viewer,dmv,velocity,ctx->outputpath,suffix,"stokes","velocity",XDMFNode,XDMFHDF5);CHKERRQ(ierr);
    }
    PetscTime(&t1);
    PetscPrintf(PETSC_COMM_WORLD,"ModelOutput_ExecuteXDMFWriter: time %1.4e <sec>\n",t1-t0);

    //ierr = XDMFWriteDataItemByReference_DMDA(viewer,dmv,"velocity","xxxxx",XDMFBinary);CHKERRQ(ierr);

    /*
    ierr = XDMFWriteAttributeFunction_DMDA(viewer,"10.0 * $0",
                                           dmv,velocity,
                                           ctx->outputpath,suffix,"stokes","velocity2",
                                           XDMFNode,XDMFBinary);CHKERRQ(ierr);
    */

    /* log timestep */
    ierr = XDMFMetaWriteTime(viewer,ctx->time);CHKERRQ(ierr);

    ierr = XDMFGridClose_DMDA(viewer);CHKERRQ(ierr);

    ierr = XDMFMetaDomainClose(viewer);CHKERRQ(ierr);
    ierr = XDMFMetaXDMFClose(&viewer);CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccess(dmstokes,X,&velocity,&pressure);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
