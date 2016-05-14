
static const char help[] = "Test checkpointing functionality of DataBucket.\n";

#include <ptatin3d.h>
#include <ptatin_init.h>
#include <data_bucket.h>
#include <data_bucket_view.h>

typedef struct {
	long int  pid;
  double    coor;
	double    weight;
  short int region;
} QPoint;

const char QPointClassName[] = "QPoint";
const char QPointFieldName[] = "QPointFuncEval";

PetscBool view_data = PETSC_FALSE;

PetscErrorCode QuadratureEvaluate(DataBucket db,PetscReal *qsum)
{
  PetscErrorCode ierr;
  PetscReal sum,sum_g;
  int k,L;
	DataField dbField_a,dbField_b;
  double *eval;
  QPoint *points;
  
  MPI_Barrier(PETSC_COMM_WORLD);
	DataBucketGetSizes(db,&L,NULL,NULL);
  
  DataBucketGetDataFieldByName(db,QPointClassName,&dbField_a);
  DataBucketGetDataFieldByName(db,QPointFieldName,&dbField_b);
  DataFieldGetEntries(dbField_a,(void**)&points);
  DataFieldGetEntries(dbField_b,(void**)&eval);
  
  sum = 0.0;
  for (k=0; k<L; k++) {
    sum += points[k].weight * eval[k];
    if (view_data) printf("** (%d,pid %ld))xi = %1.4e : w = %1.4e \n",k,points[k].pid,points[k].coor,points[k].weight);
  }
  
  DataFieldRestoreEntries(dbField_b,(void**)&eval);
  DataFieldRestoreEntries(dbField_a,(void**)&points);
  
  ierr = MPI_Allreduce(&sum,&sum_g,1,MPIU_REAL,MPIU_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
  
  *qsum = sum_g;
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "DataBucketCheckpoint"
PetscErrorCode DataBucketCheckpoint(void)
{
  DataBucket db;
  PetscInt nqp;
  int k,len_local,len_rem,L;
  PetscMPIInt commsize,rank;
  PetscErrorCode ierr;
  double *eval;
  QPoint *points;
  double dxi;
	DataField  dbField_a,dbField_b;
  long int Lg,Ag;
  PetscReal qval;
  PetscLogDouble t0,t1;
  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&commsize);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  
  DataBucketCreate(&db);
	DataBucketRegisterField(db,QPointClassName,sizeof(QPoint),NULL);
	DataBucketRegisterField(db,QPointFieldName,sizeof(double),NULL);
	DataBucketFinalize(db);

  nqp = 10;
  PetscOptionsGetInt(NULL,NULL,"-nqp",&nqp,NULL);

  len_local = nqp / commsize;
  len_rem = 0;
  if (rank == (commsize-1)) {
    len_rem = nqp - len_local * commsize;
  }
  
  DataBucketSetSizes(db,len_local + len_rem,25+rank);
  
	DataBucketGetSizes(db,&L,NULL,NULL);

  DataBucketGetDataFieldByName(db,QPointClassName,&dbField_a);
  DataBucketGetDataFieldByName(db,QPointFieldName,&dbField_b);
  DataFieldGetEntries(dbField_a,(void**)&points);
  DataFieldGetEntries(dbField_b,(void**)&eval);

  dxi = 2.0 / ((double)nqp);
  for (k=0; k<L; k++) {
    points[k].region = 20;
    points[k].pid = len_local * rank + k;
    
    points[k].weight = 2.0 / ((double)nqp);
    points[k].coor   = -1.0 + 0.5*dxi + dxi*points[k].pid;
    if (view_data) printf("[%d] (%d,pid %ld))xi = %1.4e : w = %1.4e \n",rank,k,points[k].pid,points[k].coor,points[k].weight);
    
    eval[k] = sin(0.5 * M_PI* points[k].coor + 0.12345);
  }
  
  DataFieldRestoreEntries(dbField_b,(void**)&eval);
  DataFieldRestoreEntries(dbField_a,(void**)&points);
  
  ierr = QuadratureEvaluate(db,&qval);CHKERRQ(ierr);
  DataBucketGetGlobalSizes(PETSC_COMM_WORLD,db,&Lg,NULL,&Ag);
  PetscPrintf(PETSC_COMM_WORLD,"sum = %1.12e : npoints %ld : allocated %ld\n",qval,Lg,Ag);
  
  PetscTime(&t0);
  ierr = DataBucketView_Native(db,"test_",PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscTime(&t1);
  PetscPrintf(PETSC_COMM_WORLD,"DataBucketView_Native: %1.4e (sec)\n",t1-t0);
  
  DataBucketDestroy(&db);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DataBucketRestart"
PetscErrorCode DataBucketRestart(void)
{
  DataBucket db;
  PetscErrorCode ierr;
  long int Lg,Ag;
  PetscReal qval;
  PetscLogDouble t0,t1;

  DataBucketCreate(&db);
  PetscTime(&t0);
  ierr = DataBucketLoad_Native(db,"db.meta",PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscTime(&t1);
  
  ierr = QuadratureEvaluate(db,&qval);CHKERRQ(ierr);
  DataBucketGetGlobalSizes(PETSC_COMM_WORLD,db,&Lg,NULL,&Ag);
  PetscPrintf(PETSC_COMM_WORLD,"sum = %1.12e : npoints %ld : allocated %ld\n",qval,Lg,Ag);
  PetscPrintf(PETSC_COMM_WORLD,"DataBucketLoad_Native: %1.4e (sec)\n",t1-t0);

  DataBucketDestroy(&db);
 
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
  PetscBool write = PETSC_TRUE,flg = PETSC_FALSE;
  
	ierr = pTatinInitialize(&argc,&argv,0,help);CHKERRQ(ierr);

	PetscOptionsGetBool(NULL,NULL,"-view_data",&view_data,NULL);
  PetscOptionsGetBool(NULL,NULL,"-load_db",&flg,NULL);
  if (flg) { write = PETSC_FALSE; }
  if (write) {
    PetscPrintf(PETSC_COMM_WORLD,"Writing DB\n");
    ierr = DataBucketCheckpoint();CHKERRQ(ierr);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"Loading DB\n");
    ierr = DataBucketRestart();CHKERRQ(ierr);
  }
  
	ierr = pTatinFinalize();CHKERRQ(ierr);
	return 0;
}
