
#include <private/pcimpl.h>     /*I "petscpc.h" I*/
#include <petscksp.h>           /*I "petscksp.h" I*/
#include <sub_comm.h>

typedef struct {
	PetscInt nsubcomm_factor;
	PetscMPIInt nsubcomm_size;
	MPI_Subcomm subcomm;
	Mat A,B,Ared,Bred;
	Vec xred,yred;
	KSP ksp;
	IS isin;
	VecScatter scatter;
	Vec xtmp;
} PC_SemiRedundant;


/* helpers for gather matrix and scattering vectors */

#undef __FUNCT__
#define __FUNCT__ "MatCreateSemiRedundant"
PetscErrorCode MatCreateSemiRedundant(Mat A,MPI_Subcomm subcomm,MatReuse reuse,Mat *_red)
{
	PetscErrorCode ierr;
	IS             isrow,iscol;
	PetscInt       start,end,nr,nc;
	PetscInt       i,j,*nnz,*onnz;
	MPI_Comm       comm;
	Mat            Alocal,*_Alocal,red;
	
	
	ierr = MatGetSize(A,&nr,&nc);CHKERRQ(ierr);
	ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
	
	ierr = ISCreateStride(comm,nc,0,1,&iscol);CHKERRQ(ierr);
	
	start = end = 0;
	if (subcomm->parent_rank_active_in_subcomm) {

		if (reuse == MAT_INITIAL_MATRIX) {
			ierr = MatCreate(subcomm->sub_comm,&red);CHKERRQ(ierr);
			ierr = MatSetSizes(red,PETSC_DECIDE,PETSC_DECIDE,nr,nc);CHKERRQ(ierr);
			ierr = MatSetFromOptions(red);CHKERRQ(ierr);
		} else {
			red = *_red;
		}
		
		ierr = MatGetOwnershipRange(red,&start,&end);CHKERRQ(ierr);
		ierr = ISCreateStride(comm,(end-start),start,1,&isrow);CHKERRQ(ierr);
	} else {
		/* if rank not in subcomm, just fetch a local chunk of A */
		ierr = MatGetOwnershipRange(A,&start,&end);CHKERRQ(ierr);
		ierr = ISCreateStride(comm,1,start,1,&isrow);CHKERRQ(ierr);
	}
	
	ierr = MatGetSubMatrices(A,1,&isrow,&iscol,MAT_INITIAL_MATRIX,&_Alocal);CHKERRQ(ierr);
	Alocal = *_Alocal;
	
	/* insert entries */
	if (subcomm->parent_rank_active_in_subcomm) {
		PetscInt ncols,startc,endc,rowidx;
		const PetscInt *cols;
		const PetscScalar *vals;
		
		/* preallocation */
		ierr = MatGetOwnershipRange(red,&start,&end);CHKERRQ(ierr);
		ierr = MatGetOwnershipRangeColumn(red,&startc,&endc);CHKERRQ(ierr);
		
		PetscMalloc(sizeof(PetscInt)*(end-start),&nnz);
		PetscMalloc(sizeof(PetscInt)*(end-start),&onnz);
		PetscMemzero(nnz,sizeof(PetscInt)*(end-start));
		PetscMemzero(onnz,sizeof(PetscInt)*(end-start));
		
		for (i=0; i<(end-start); i++) {
			ierr = MatGetRow(Alocal,i,&ncols,&cols,PETSC_NULL);CHKERRQ(ierr);
			for (j=0; j<ncols; j++) {
				if ( (cols[j] >= startc) && (cols[j] < endc) ) {
					nnz[i]++;
				} else {
					onnz[i]++;
				}
			}
			ierr = MatRestoreRow(Alocal,i,&ncols,&cols,PETSC_NULL);CHKERRQ(ierr);
		}
		ierr = MatSeqAIJSetPreallocation(red,PETSC_NULL,nnz);CHKERRQ(ierr);
		ierr = MatMPIAIJSetPreallocation(red,PETSC_NULL,nnz,PETSC_NULL,onnz);CHKERRQ(ierr);
		
		PetscFree(nnz);
		PetscFree(onnz);
		
		/* insert */
		ierr = MatGetOwnershipRange(red,&start,&end);CHKERRQ(ierr);
		for (i=0; i<(end-start); i++) {
			ierr = MatGetRow(Alocal,i,&ncols,&cols,&vals);CHKERRQ(ierr);
			
			rowidx = i + start;
			ierr = MatSetValues(red,1,&rowidx,ncols,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
			
			ierr = MatRestoreRow(Alocal,i,&ncols,&cols,PETSC_NULL);CHKERRQ(ierr);
		}
		
		ierr = MatAssemblyBegin(red,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(red,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	
	ierr = ISDestroy(&isrow);CHKERRQ(ierr);
	ierr = ISDestroy(&iscol);CHKERRQ(ierr);
	ierr = MatDestroy(&Alocal);CHKERRQ(ierr);
	
	*_red = PETSC_NULL;
	if (subcomm->parent_rank_active_in_subcomm) {
		*_red = red;
	}
	
	return(0);
}

/* implementations for SemiRedundant */

#undef __FUNCT__  
#define __FUNCT__ "PCSetUp_SemiRedundant"
static PetscErrorCode PCSetUp_SemiRedundant(PC pc)
{
  PetscErrorCode   ierr;
  PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
	MPI_Comm comm;
	MPI_Subcomm sc;
	PetscMPIInt size;
	MatStructure str;
  MatReuse reuse;

	
	/* construction phase */
  if (!pc->setupcalled) {
	
		/* set up sub communicator */
		PetscObjectGetComm((PetscObject)pc,&comm);
		MPI_Subcomm_create_MethodA(comm,red->nsubcomm_factor,&sc);
		red->subcomm = sc;
		ierr = MPI_Comm_size(sc->sub_comm,&red->nsubcomm_size);CHKERRQ(ierr);
		
		red->Ared = PETSC_NULL;
		red->Bred = PETSC_NULL;
		ierr = PCGetOperators(pc,&red->A,&red->B,&str);CHKERRQ(ierr);
		
		red->ksp = PETSC_NULL;
		if (red->subcomm->parent_rank_active_in_subcomm) {
			const char     *prefix;

			ierr = KSPCreate(red->subcomm->sub_comm,&red->ksp);CHKERRQ(ierr);
			ierr = PetscObjectIncrementTabLevel((PetscObject)red->ksp,(PetscObject)pc,1);CHKERRQ(ierr);
			ierr = PetscLogObjectParent(pc,red->ksp);CHKERRQ(ierr);

			ierr = PCGetOptionsPrefix(pc,&prefix);CHKERRQ(ierr);
			ierr = KSPSetOptionsPrefix(red->ksp,prefix);CHKERRQ(ierr); 
			ierr = KSPAppendOptionsPrefix(red->ksp,"semiredundant_");CHKERRQ(ierr); 
		}
	}
	
	
	
	
	/* fetch redundant matrix */
	if (!pc->setupcalled) {

		ierr = MatCreateSemiRedundant(red->A,red->subcomm,MAT_INITIAL_MATRIX,&red->Ared);CHKERRQ(ierr);
		if (red->A != red->B) {
			ierr = MatCreateSemiRedundant(red->B,red->subcomm,MAT_INITIAL_MATRIX,&red->Bred);CHKERRQ(ierr);
		}

		if ( (red->A == red->B) && (red->Ared) ) {			
			red->Bred = red->Ared;
			ierr = PetscObjectReference((PetscObject)red->Ared);CHKERRQ(ierr);
		}
		
		red->xred = PETSC_NULL;
		red->yred = PETSC_NULL;
		if (red->Ared) {
			PetscInt m,n;
			
			ierr = MatGetLocalSize(red->Ared,&m,&n);CHKERRQ(ierr);
			/* create xred with empty local arrays, because xdup's arrays will be placed into it */
			ierr = VecCreateMPIWithArray(red->subcomm->sub_comm,m,PETSC_DECIDE,PETSC_NULL,&red->xred);CHKERRQ(ierr);

			ierr = MatGetVecs(red->Ared,PETSC_NULL,&red->yred);CHKERRQ(ierr);
		}
		
	} else {
		reuse = MAT_REUSE_MATRIX;

		ierr = MatCreateSemiRedundant(red->A,red->subcomm,reuse,&red->Ared);CHKERRQ(ierr);
		if (red->A != red->B) {
			ierr = MatCreateSemiRedundant(red->B,red->subcomm,reuse,&red->Bred);CHKERRQ(ierr);
		} else {
			red->Bred = red->Ared;
		}
	}
	
	/* setup scatters */
	if (!pc->setupcalled) {
		PetscInt st,ed;
		PetscInt i,n,N;
		Vec x;
	
		PetscObjectGetComm((PetscObject)pc,&comm);
		ierr = MatGetVecs(red->A,&x,PETSC_NULL);CHKERRQ(ierr);
		
		if (red->xred) {
			ierr = VecGetOwnershipRange(red->xred,&st,&ed);CHKERRQ(ierr);
			ierr = ISCreateStride(comm,ed-st,st,1,&red->isin);CHKERRQ(ierr);
		} else {
			ierr = VecGetOwnershipRange(x,&st,&ed);CHKERRQ(ierr);
			ierr = ISCreateStride(comm,1,st,1,&red->isin);CHKERRQ(ierr);
		}
		
		ierr = ISGetLocalSize(red->isin,&n);CHKERRQ(ierr);
		ierr = ISGetSize(red->isin,&N);CHKERRQ(ierr);
		ierr = VecCreate(((PetscObject)red->isin)->comm,&red->xtmp);CHKERRQ(ierr);
		ierr = VecSetSizes(red->xtmp,n,N);CHKERRQ(ierr);
		ierr = VecSetType(red->xtmp,((PetscObject)x)->type_name);CHKERRQ(ierr);
		ierr = VecScatterCreate(x,red->isin,red->xtmp,PETSC_NULL,&red->scatter);CHKERRQ(ierr);

		ierr = VecDestroy(&x);CHKERRQ(ierr);
	}
	
	
	/* common - no construction */
	ierr = PCGetOperators(pc,&red->A,&red->B,&str);CHKERRQ(ierr);
	if (red->Ared) {
		//MatView(red->Ared,PETSC_VIEWER_STDOUT_(red->subcomm->sub_comm));
		ierr = KSPSetOperators(red->ksp,red->Ared,red->Bred,str);CHKERRQ(ierr);
		
		if (pc->setfromoptionscalled){
			ierr = KSPSetFromOptions(red->ksp);CHKERRQ(ierr); 
		}
		
		ierr = KSPSetUp(red->ksp);CHKERRQ(ierr);
	}
	
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "PCApply_SemiRedundant"
static PetscErrorCode PCApply_SemiRedundant(PC pc,Vec x,Vec y)
{
  PetscErrorCode   ierr;
  PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
	
	Vec xtmp;
	PetscInt i,st,ed;
	VecScatter scatter;
  PetscScalar    *array;
	
	
  PetscFunctionBegin;
	
	xtmp = red->xtmp;
	scatter = red->scatter;
	
	/* pull in vector */
	ierr = VecScatterBegin(scatter,x,xtmp,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
	ierr = VecScatterEnd(scatter,x,xtmp,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
	
	/* solve */
	if (red->ksp) {

		/* define xred */
		ierr = VecGetArray(xtmp,&array);CHKERRQ(ierr);
		/* we created xred with empty local arrays, now we fill it in */
		ierr = VecPlaceArray(red->xred,(const PetscScalar*)array);CHKERRQ(ierr);
		
		ierr = KSPSolve(red->ksp,red->xred,red->yred);CHKERRQ(ierr);
		
		ierr = VecResetArray(red->xred);CHKERRQ(ierr);
		ierr = VecRestoreArray(xtmp,&array);CHKERRQ(ierr);
	}

	/* return vector */
	ierr = VecGetArray(xtmp,&array);CHKERRQ(ierr);
	if (red->yred) {
		PetscScalar *LA_yred;
		
		VecGetOwnershipRange(red->yred,&st,&ed);

		ierr = VecGetArray(red->yred,&LA_yred);CHKERRQ(ierr);
		for (i=0; i<ed-st; i++) {
			array[i] = LA_yred[i];
		}
		ierr = VecRestoreArray(red->yred,&LA_yred);CHKERRQ(ierr);
	} else {
		array[0] = 0.0;
	}
	ierr = VecRestoreArray(xtmp,&array);CHKERRQ(ierr);
	
	ierr = VecScatterBegin(scatter,xtmp,y,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
	ierr = VecScatterEnd(scatter,xtmp,y,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PCReset_SemiRedundant"
static PetscErrorCode PCReset_SemiRedundant(PC pc)
{
  PetscErrorCode   ierr;
  PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
	
  PetscFunctionBegin;
	if (red->xred) { ierr = VecDestroy(&red->xred);CHKERRQ(ierr);}
	if (red->yred) { ierr = VecDestroy(&red->yred);CHKERRQ(ierr);}
	if (red->Ared) { ierr = MatDestroy(&red->Ared);CHKERRQ(ierr);}
	if (red->Bred) { ierr = MatDestroy(&red->Bred);CHKERRQ(ierr);}

	ierr = ISDestroy(&red->isin);CHKERRQ(ierr);
	ierr = VecDestroy(&red->xtmp);CHKERRQ(ierr);
	ierr = VecScatterDestroy(&red->scatter);CHKERRQ(ierr);
  if (red->ksp) {  ierr = KSPReset(red->ksp);CHKERRQ(ierr);}
  
	PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PCDestroy_SemiRedundant"
static PetscErrorCode PCDestroy_SemiRedundant(PC pc)
{
  PetscErrorCode   ierr;
  PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
	
  PetscFunctionBegin;
  ierr = PCReset_SemiRedundant(pc);CHKERRQ(ierr);
	if (red->ksp) {
		ierr = KSPDestroy(&red->ksp);CHKERRQ(ierr);
	}
	ierr = MPI_Subcomm_free(&red->subcomm);CHKERRQ(ierr);
  ierr = PetscFree(pc->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PCSetFromOptions_SemiRedundant"
static PetscErrorCode PCSetFromOptions_SemiRedundant(PC pc)
{
  PetscErrorCode   ierr;
  PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
	
  PetscFunctionBegin;
  ierr = PetscOptionsHead("SemiRedundant options");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-pc_semiredundant_factor","Factor to reduce parent communication size by","PCSemiRedundantSetFactor",red->nsubcomm_factor,&red->nsubcomm_factor,0);CHKERRQ(ierr);
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PCView_SemiRedundant"
static PetscErrorCode PCView_SemiRedundant(PC pc,PetscViewer viewer)
{
  PetscErrorCode   ierr;
  PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
  PetscBool        iascii,isstring;
  PetscViewer    subviewer;
	
  PetscFunctionBegin;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  ierr = PetscTypeCompare((PetscObject)viewer,PETSCVIEWERSTRING,&isstring);CHKERRQ(ierr);
  if (iascii) {
		
		if (!red->subcomm) {
      ierr = PetscViewerASCIIPrintf(viewer,"  SemiRedundant preconditioner: Not yet setup\n");CHKERRQ(ierr);
		} else {

			/* ierr = PetscViewerASCIIPrintf(viewer,"  SemiRedundant preconditioner:\n");CHKERRQ(ierr); */
      ierr = PetscViewerASCIIPrintf(viewer,"  SemiRedundant: parent comm size reduction factor = %D\n",red->nsubcomm_factor);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"  SemiRedundant: subcomm_size = %D\n",red->nsubcomm_size);CHKERRQ(ierr);
			
			ierr = PetscViewerGetSubcomm(viewer,red->subcomm->sub_comm,&subviewer);CHKERRQ(ierr);
			ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
			if (red->subcomm->parent_rank_active_in_subcomm) {
				ierr = KSPView(red->ksp,subviewer);CHKERRQ(ierr);
			}
			ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
			ierr = PetscViewerRestoreSubcomm(viewer,red->subcomm->sub_comm,&subviewer);CHKERRQ(ierr);
			
		}
  } else if (isstring) { 
		ierr = PetscViewerStringSPrintf(viewer," SemiRedundant preconditioner");CHKERRQ(ierr);
  } else {
    SETERRQ1(((PetscObject)pc)->comm,PETSC_ERR_SUP,"Viewer type %s not supported for PC SemiRedundant",((PetscObject)viewer)->type_name);
  }
  PetscFunctionReturn(0);
}



EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PCCreate_SemiRedundant"
PetscErrorCode PCCreate_SemiRedundant(PC pc)
{
  PetscErrorCode   ierr;
  PC_SemiRedundant *red;
  PetscMPIInt      size;
  
  PetscFunctionBegin;
  ierr = PetscNewLog(pc,PC_SemiRedundant,&red);CHKERRQ(ierr);
  pc->data            = (void*)red; 
	
  red->nsubcomm_factor = 1;
  ierr = MPI_Comm_size(((PetscObject)pc)->comm,&size);CHKERRQ(ierr);
  red->nsubcomm_size   = size;
	
  pc->ops->apply           = PCApply_SemiRedundant;
  pc->ops->applytranspose  = 0;
  pc->ops->setup           = PCSetUp_SemiRedundant;
  pc->ops->destroy         = PCDestroy_SemiRedundant;
  pc->ops->reset           = PCReset_SemiRedundant;
  pc->ops->setfromoptions  = PCSetFromOptions_SemiRedundant;
  pc->ops->view            = PCView_SemiRedundant;    
	
  PetscFunctionReturn(0);
}
EXTERN_C_END

