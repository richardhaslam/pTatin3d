
#ifndef __ptat3d_sub_comm_h__
#define __ptat3d_sub_comm_h__

typedef struct _p_MPI_Subcomm *MPI_Subcomm;

struct _p_MPI_Subcomm { 
	MPI_Comm parent_comm;
	MPI_Comm sub_comm;
	int nranks_from_parent;
	int *ranks_from_parent;
	int parent_rank_active_in_subcomm; /* 1:true, 0:false */
};

int MPI_Subcomm_free(MPI_Subcomm *scomm);
int MPI_Subcomm_get_active(MPI_Subcomm sc,int *a);
int MPI_Subcomm_get_parent_comm(MPI_Subcomm sc,MPI_Comm *a);
int MPI_Subcomm_get_comm(MPI_Subcomm sc,MPI_Comm *a);
int MPI_Subcomm_get_active_ranks(MPI_Subcomm sc,int **a);

int MPI_Subcomm_create_MethodA(MPI_Comm parent_comm,int parent_reduction_factor,MPI_Subcomm *scomm);


#endif
