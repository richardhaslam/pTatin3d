
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "sub_comm.h"

int _MPI_Subcomm_create(MPI_Subcomm *scomm)
{
	MPI_Subcomm comm;
	
	comm = malloc(sizeof(struct _p_MPI_Subcomm));
	memset(comm,0,sizeof(struct _p_MPI_Subcomm));
	
	*scomm = comm;
	
	return(0);
}

int MPI_Subcomm_free(MPI_Subcomm *scomm)
{
	MPI_Subcomm comm;
	int ierr;
	
	if (scomm) { comm = *scomm; }
	
	if (comm->sub_comm) { ierr = MPI_Comm_free(&comm->sub_comm); }
	if (comm->ranks_from_parent) { free(comm->ranks_from_parent); }
	free(comm);
	
	*scomm = NULL;
	
	return(0);
}

int MPI_Subcomm_get_active(MPI_Subcomm sc,int *a)
{
	*a = sc->parent_rank_active_in_subcomm;
	return(0);
}

int MPI_Subcomm_get_parent_comm(MPI_Subcomm sc,MPI_Comm *a)
{
	*a = sc->parent_comm;
	return(0);
}

int MPI_Subcomm_get_comm(MPI_Subcomm sc,MPI_Comm *a)
{
	*a = sc->sub_comm;
	return(0);
}

int MPI_Subcomm_get_active_ranks(MPI_Subcomm sc,int **a)
{
	*a = sc->ranks_from_parent;
	return(0);
}


int MPI_Subcomm_create_MethodA(MPI_Comm parent_comm,int parent_reduction_factor,MPI_Subcomm *scomm)
{
	MPI_Subcomm comm;
	MPI_Comm sub_comm;
	MPI_Group parent_group,sub_group;
	int nproc,rank,*subranks,nsubranks,c,i,active;
	int ierr;

	if (parent_reduction_factor < 1) {
		parent_reduction_factor = 1;
		printf("Warning:MPI_Subcomm_create_MethodA: parent_reduction_factor >=1\n");
	}
	
	ierr = MPI_Comm_size(parent_comm,&nproc);
	ierr = MPI_Comm_rank(parent_comm,&rank);
	nsubranks = 0;
	for (i=0; i<nproc; i++) {
		if (i%parent_reduction_factor == 0) {
			nsubranks++;
		}
	}
	subranks = malloc(sizeof(int)*nsubranks);
	c = 0;
	for (i=0; i<nproc; i++) {
		if (i%parent_reduction_factor == 0) {
			subranks[c] = i;
			c++;
		}
	}
	active = 0;
	for (i=0; i<nsubranks; i++) {
		if (rank == subranks[i]) {
			active = 1;
		}
	}
	
	ierr = _MPI_Subcomm_create(&comm);
	ierr = MPI_Comm_group(parent_comm,&parent_group);
	sub_group = NULL;
	if (active == 1) {
		ierr = MPI_Group_incl(parent_group, nsubranks, subranks, &sub_group);
	} else {
		ierr = MPI_Group_excl(parent_group, nsubranks, subranks, &sub_group);
	}
	ierr = MPI_Comm_create(parent_comm, sub_group, &sub_comm);

	{
		int sr,snp;
		ierr = MPI_Comm_size(sub_comm,&snp);
		ierr = MPI_Comm_rank(sub_comm,&sr);
		printf("parent[%d of %d]: sub[%d of %d]: active = %d \n",rank,nproc,sr,snp,active);
	}
	
	comm->parent_comm        = parent_comm;
	comm->sub_comm           = sub_comm;
	comm->nranks_from_parent = nsubranks;
	comm->ranks_from_parent  = subranks;
	comm->parent_rank_active_in_subcomm = active;

	/* should be able to free group safely as its been embedded insode sub_comm */
	ierr = MPI_Group_free(&sub_group);
	*scomm = comm;
	
	return(0);
}


