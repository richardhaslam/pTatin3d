/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    pTatin3d
 **    filename:   data_bucket.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
/*

gcc -O3 -g -c data_bucket.c

 // GENERIC //
 // insert into an exisitng location
 SwarmFieldInsertPoint( field, index, ctx ) 
 // remove data at index - replace with last point
 SwarmFieldRemovePoint( field, index, ctx )


 // USER PROVIDED //
 ParticleDataLoad_XXX
 ParticleDataView_XXX
 
 Usage:
 
 
 DataBucketGetKeyByName(db,"particle_data_standard",&K);
 DataFieldGetAccess( K, &gfield );
 DataFieldVerifyAccess( gfield,sizeof(ParticleData_Standard));

 // LOOP
 DataFieldAccessPoint(gfield,index,(void**)&stdParticle);

 // REPLACE 
 ParticleData_Standard my_point;
 ParticleDataLoad_Standard( &my_point, LOAD_CRAP_HERE );
 SwarmFieldInsertPoint( field, index, &my_point );
 
 DataFieldRestoreAccess( K, &gfield );
 
 // EXPAND - SHRINK set via 
 DataBucketSetSizes(db,newvalue, buffer); // -ve buffer val retains old value //
 
 
 
 Proposed memory management for adding / deleting particles
 
 A: active
 B: buffer / deactive
 D: deactive - formerly active point which was deleted

 initial memory layout
 AAAAAAAAAAAAAAAAA BBBBB

 1) fill up all buffer spots, then we reallocate
 AAAAAAAAAAAAAAAAA AAAAA (reallocate larger)
 AAAAAAAAAAAAAAAAA AAAAA BBBBB [RESULT]

 
 2) deactivate  buffer spots, then we reallocate 
 AAAAAAAAAAAAAA DDD BBBBB
 AAAAAAAAAAAA DDDDD BBBBB (reallocate smaller)
 AAAAAAAAAAAA BBBBB [RESULT]
 
*/


#include <data_bucket.h>
#include <mpiio_blocking.h>
#include <cjson_utils.h>

/* string helpers */
static void DataBucketStringInList( const char name[], const int N, const DataField gfield[], BTruth *val )
{
	int i;
	
	*val = BFALSE;
	for( i=0; i<N; i++ ) {
		if( strcmp( name, gfield[i]->name ) == 0 ) {
			*val = BTRUE;
			return;
		}
	}
}

static void DataBucketStringFindInList( const char name[], const int N, const DataField gfield[], int *index )
{
	int i;
	
	*index = -1;
	for( i=0; i<N; i++ ) {
		if( strcmp( name, gfield[i]->name ) == 0 ) {
			*index = i;
			return;
		}
	}
}

void DataFieldCreate( const char registration_function[], const char name[], const size_t size, const int L, DataField *DF )
{
	DataField df;
	
	df = malloc( sizeof(struct _p_DataField) );
	memset( df, 0, sizeof(struct _p_DataField) ); 
	
	if (asprintf( &df->registration_function, "%s", registration_function ) < 0) {printf("asprintf() error. Exiting ungracefully.\n"); exit(1);}
	if (asprintf( &df->name, "%s", name ) < 0) {printf("asprintf() error. Exiting ungracefully.\n"); exit(1);}
	df->atomic_size = size;
	df->L = L;
	
	df->data = malloc( size * L ); /* allocate something so we don't have to reallocate */
	memset( df->data, 0, size * L );
	
	*DF = df;
}

void DataFieldDestroy( DataField *DF )
{
	DataField df = *DF;
	
	free( df->registration_function );
	free( df->name );
	free( df->data );
	free(df);
	
	*DF = NULL;
}

/* data bucket */
void DataBucketCreate( DataBucket *DB )
{
	DataBucket db;
	
	
	db = malloc( sizeof(struct _p_DataBucket) );
	memset( db, 0, sizeof(struct _p_DataBucket) );

	db->finalised = BFALSE;
	
	/* create empty spaces for fields */
	db->L         = 0;
	db->buffer    = 1;
	db->allocated = 1;

	db->nfields   = 0;
	db->field     = malloc(sizeof(DataField));
	
	*DB = db;
}

void DataBucketDestroy( DataBucket *DB )
{
	DataBucket db = *DB;
	int f;
	
	/* release fields */
	for( f=0; f<db->nfields; f++ ) {
		DataFieldDestroy(&db->field[f]);	
	}

	/* this will catch the initially allocated objects in the event that no fields are registered */
	if(db->field!=NULL) {
		free(db->field);
	}
	
	free(db);
	
	*DB = NULL;
}

void _DataBucketRegisterField(
						DataBucket db,
						const char registration_function[],
						const char field_name[],
						size_t atomic_size, DataField *_gfield )
{
	BTruth val;
	DataField *field,fp;
	
	/* check we haven't finalised the registration of fields */
	/*
	if(db->finalised==BTRUE) {
		printf("ERROR: DataBucketFinalize() has been called. Cannot register more fields\n");
		ERROR();
	}
	*/
	 
	/* check for repeated name */
	DataBucketStringInList( field_name, db->nfields, (const DataField*)db->field, &val );
	if(val == BTRUE ) {
		printf("ERROR: Cannot add same field twice\n");
		ERROR();
	}

	/* create new space for data */
	field = realloc( db->field,     sizeof(DataField)*(db->nfields+1));
	db->field     = field;
	
	/* add field */
	DataFieldCreate( registration_function, field_name, atomic_size, db->allocated, &fp );
	db->field[ db->nfields ] = fp;
	
	db->nfields++;
	
	if(_gfield!=NULL){
		*_gfield = fp;
	}
}

void DataBucketGetDataFieldByName(DataBucket db,const char name[],DataField *gfield)
{
	int idx;
	BTruth found;
	
	DataBucketStringInList(name,db->nfields,(const DataField*)db->field,&found);
	if(found==BFALSE) {
		printf("ERROR: Cannot find DataField with name %s \n", name );
		ERROR();
	}
	DataBucketStringFindInList(name,db->nfields,(const DataField*)db->field,&idx);
		
	*gfield = db->field[idx];
}

void DataBucketQueryDataFieldByName(DataBucket db,const char name[],BTruth *found)
{
	*found = BFALSE;
	DataBucketStringInList(name,db->nfields,(const DataField*)db->field,found);
}

void DataBucketFinalize(DataBucket db)
{
	db->finalised = BTRUE;
}

void DataFieldGetNumEntries(DataField df, int *sum)
{
	*sum = df->L;
}

void DataFieldSetSize( DataField df, const int new_L )
{
	void *tmp_data;
	
	if( new_L <= 0 ) {
		printf("ERROR: Cannot set size of DataField to be <= 0 \n");
		ERROR();
	}
	if( new_L == df->L ) return;
	
	if( new_L > df->L ) {
		
		tmp_data = realloc( df->data, df->atomic_size * (new_L) );
		df->data = tmp_data;
		
		/* init new contents */
		memset( ( ((char*)df->data)+df->L*df->atomic_size), 0, (new_L-df->L)*df->atomic_size );
		
	}
	else {
		/* reallocate pointer list, add +1 in case new_L = 0 */
		tmp_data = realloc( df->data, df->atomic_size * (new_L+1) );
		df->data = tmp_data;
	}
	
	df->L = new_L;
}

void DataFieldZeroBlock( DataField df, const int start, const int end )
{
	if( start > end ) {
		printf("ERROR: Cannot zero a block of entries if start(%d) > end(%d) \n",start,end);
		ERROR();
	}
	if( start < 0 ) {
		printf("ERROR: Cannot zero a block of entries if start(%d) < 0 \n",start);
		ERROR();
	}
	if( end > df->L ) {
		printf("ERROR: Cannot zero a block of entries if end(%d) >= array size(%d) \n",end,df->L);
		ERROR();
	}
	
	memset( ( ((char*)df->data)+start*df->atomic_size), 0, (end-start)*df->atomic_size );
}

/*
 A negative buffer value will simply be ignored and the old buffer value will be used.
 */
void DataBucketSetSizes( DataBucket db, const int L, const int buffer )
{
	int current_allocated,new_used,new_unused,new_buffer,new_allocated,f;
	
	
	if( db->finalised == BFALSE ) {
		printf("ERROR: You must call DataBucketFinalize() before DataBucketSetSizes() \n");
		ERROR();
	}
	 
	current_allocated = db->allocated;
	
	new_used   = L;
	new_unused = current_allocated - new_used;
	new_buffer = db->buffer;
	if( buffer >= 0 ) { /* update the buffer value */
		new_buffer = buffer;
	}
	new_allocated = new_used + new_buffer;
	
	/* action */
	if ( new_allocated > current_allocated ) {
		/* increase size to new_used + new_buffer */
		for( f=0; f<db->nfields; f++ ) {
			DataFieldSetSize( db->field[f], new_allocated );
		}
		
		db->L         = new_used;
		db->buffer    = new_buffer;
		db->allocated = new_used + new_buffer;
	}
	else {
		if( new_unused > 2 * new_buffer ) {
			
			/* shrink array to new_used + new_buffer */
			for( f=0; f<db->nfields; f++ ) {
				DataFieldSetSize( db->field[f], new_allocated );
			}
			
			db->L         = new_used;
			db->buffer    = new_buffer;
			db->allocated = new_used + new_buffer;
		}
		else {
			db->L      = new_used;
			db->buffer = new_buffer;
		}
	}
	
	/* zero all entries from db->L to db->allocated */
	for( f=0; f<db->nfields; f++ ) {
		DataField field = db->field[f];
		DataFieldZeroBlock(field, db->L,db->allocated);
	}
}

void DataBucketSetInitialSizes( DataBucket db, const int L, const int buffer )
{
	int f;
	DataBucketSetSizes(db,L,buffer);
	
	for( f=0; f<db->nfields; f++ ) {
		DataField field = db->field[f];
		DataFieldZeroBlock(field,0,db->allocated);
	}
}

void DataBucketGetSizes( DataBucket db, int *L, int *buffer, int *allocated )
{
	if (L) { *L = db->L; }
	if (buffer) { *buffer = db->buffer; }
	if (allocated) { *allocated = db->allocated; }
}

void DataBucketGetGlobalSizes(MPI_Comm comm, DataBucket db, long int *L, long int *buffer, long int *allocated )
{
	int      ierr;
	long int _L,_buffer,_allocated;
	
	_L = (long int)db->L;
	_buffer = (long int)db->buffer;
	_allocated = (long int)db->allocated;
	
	if (L) {         
    ierr = MPI_Allreduce(&_L,L,1,MPI_LONG,MPI_SUM,comm);MPI_ERROR_CHECK(comm,ierr);
  }
	if (buffer) {    
    ierr = MPI_Allreduce(&_buffer,buffer,1,MPI_LONG,MPI_SUM,comm);MPI_ERROR_CHECK(comm,ierr);
  }
	if (allocated) { 
    ierr =  MPI_Allreduce(&_allocated,allocated,1,MPI_LONG,MPI_SUM,comm);MPI_ERROR_CHECK(comm,ierr);
  }
}

void DataBucketGetDataFields( DataBucket db, int *L, DataField *fields[] )
{
	if(L){      *L      = db->nfields; }
	if(fields){ *fields = db->field; }
}

void DataFieldGetAccess( const DataField gfield )
{
	if(gfield->active==BTRUE) {
		printf("ERROR: Field \"%s\" is already active. You must call DataFieldRestoreAccess()\n", gfield->name );
		ERROR();
	}
	gfield->active = BTRUE;
}

void DataFieldAccessPoint( const DataField gfield, const int pid, void **ctx_p )
{
#ifdef PTATIN_DATAFIELD_POINT_ACCESS_GUARD
	/* debug mode */
	/* check point is valid */
	if( pid < 0 ){ printf("ERROR: index must be >= 0\n"); ERROR();  }
	if( pid >= gfield->L ){ printf("ERROR: index must be < %d\n",gfield->L); ERROR(); }

	if(gfield->active==BFALSE) {
		printf("ERROR: Field \"%s\" is not active. You must call DataFieldGetAccess() before point data can be retrivied\n",gfield->name);
		ERROR();
	}
#endif
	
	//*ctx_p  = (void*)( ((char*)gfield->data) + pid * gfield->atomic_size);
	*ctx_p = __DATATFIELD_point_access(gfield->data,pid,gfield->atomic_size);
}

void DataFieldAccessPointOffset( const DataField gfield, const size_t offset, const int pid, void **ctx_p )
{
#ifdef PTATIN_DATAFIELD_POINT_ACCESS_GUARD
	/* debug mode */
	
	/* check point is valid */
	/* if( offset < 0 ){ printf("ERROR: offset must be >= 0\n"); ERROR();  } *//* Note compiler realizes this can never happen with an unsigned int */
	if( offset >= gfield->atomic_size ){ printf("ERROR: offset must be < %zu\n",gfield->atomic_size); ERROR(); }
	
	/* check point is valid */
	if( pid < 0 ){ printf("ERROR: index must be >= 0\n"); ERROR();  }
	if( pid >= gfield->L ){ printf("ERROR: index must be < %d\n",gfield->L); ERROR(); }

	if(gfield->active==BFALSE) {
		printf("ERROR: Field \"%s\" is not active. You must call DataFieldGetAccess() before point data can be retrivied\n",gfield->name);
		ERROR();
	}
#endif
	
	*ctx_p = __DATATFIELD_point_access_offset(gfield->data,pid,gfield->atomic_size,offset);
}

void DataFieldRestoreAccess( DataField gfield )
{
	if(gfield->active==BFALSE) {
		printf("ERROR: Field \"%s\" is not active. You must call DataFieldGetAccess()\n", gfield->name );
		ERROR();
	}
	gfield->active = BFALSE;
}

void DataFieldVerifyAccess( const DataField gfield, const size_t size)
{
#ifdef PTATIN_DATAFIELD_POINT_ACCESS_GUARD
	if(gfield->atomic_size != size ) {
        printf("ERROR: Field \"%s\" must be mapped to %zu bytes, your intended structure is %zu bytes in length.\n",
               gfield->name, gfield->atomic_size, size );
		ERROR();
	}
#endif
}

void DataFieldGetAtomicSize(const DataField gfield,size_t *size)
{
    if (size) { *size = gfield->atomic_size; }
}

void DataFieldGetEntries(const DataField gfield,void **data)
{
    if (data) {
        *data = gfield->data;
    }
}

void DataFieldRestoreEntries(const DataField gfield,void **data)
{
    if (data) {
        *data = NULL;
    }
}

/* y = x */
void DataBucketCopyPoint( const DataBucket xb, const int pid_x,
												  const DataBucket yb, const int pid_y )
{
	int f;
	for( f=0; f<xb->nfields; f++ ) {
		void *dest;
		void *src;
		
		DataFieldGetAccess( xb->field[f] );
		if (xb!=yb) { DataFieldGetAccess( yb->field[f] ); }
		
		DataFieldAccessPoint( xb->field[f],pid_x, &src );
		DataFieldAccessPoint( yb->field[f],pid_y, &dest );
		
		memcpy( dest, src, xb->field[f]->atomic_size );
		
		DataFieldRestoreAccess( xb->field[f] );
		if (xb!=yb) { DataFieldRestoreAccess( yb->field[f] ); }
	}
	
}

void DataBucketCreateFromSubset( DataBucket DBIn, const int N, const int list[], DataBucket *DB )
{
	int nfields;
	DataField *fields;
	DataBucketCreate(DB);
	int f,L,buffer,allocated,p;
	
	/* copy contents of DBIn */
	DataBucketGetDataFields(DBIn,&nfields,&fields);
	DataBucketGetSizes(DBIn,&L,&buffer,&allocated);
	
	for(f=0;f<nfields;f++) {
		DataBucketRegisterField(*DB,fields[f]->name,fields[f]->atomic_size,NULL);
	}
	DataBucketFinalize(*DB);
	
	DataBucketSetSizes(*DB,L,buffer);
	
	/* now copy the desired guys from DBIn => DB */
	for( p=0; p<N; p++ ) {
		DataBucketCopyPoint(DBIn,list[p], *DB,p);
	}
	
}

// insert into an exisitng location
void DataFieldInsertPoint( const DataField field, const int index, const void *ctx ) 
{

#ifdef PTATIN_DATAFIELD_POINT_ACCESS_GUARD
	/* check point is valid */
	if( index < 0 ){ printf("ERROR: index must be >= 0\n"); ERROR();  }
	if( index >= field->L ){ printf("ERROR: index must be < %d\n",field->L); ERROR(); }
#endif
	
//	memcpy( (void*)((char*)field->data + index*field->atomic_size), ctx, field->atomic_size );
	memcpy( __DATATFIELD_point_access(field->data,index,field->atomic_size), ctx, field->atomic_size );
}

// remove data at index - replace with last point
void DataBucketRemovePointAtIndex( const DataBucket db, const int index )
{
	int f;
	
#ifdef PTATIN_DATAFIELD_POINT_ACCESS_GUARD
	/* check point is valid */
	if( index < 0 ){ printf("ERROR: index must be >= 0\n"); ERROR(); }
	if( index >= db->allocated ){ printf("ERROR: index must be < %d\n",db->L+db->buffer); ERROR(); }
#endif	
	
	if (index >= db->L) { /* this point is not in the list - no need to error, but I will anyway */
		printf("ERROR: You should not be trying to remove point at index=%d since it's < db->L = %d \n", index, db->L );
		ERROR();
	}
	
#if 0	
	if (index == db->L-1) { /* last point in list */
		for( f=0; f<db->nfields; f++ ) {
			DataField field = db->field[f];

			DataFieldZeroPoint(field,index);
		}
	}
	else {
		for( f=0; f<db->nfields; f++ ) {
			DataField field = db->field[f];

			/* copy then remove */
			DataFieldCopyPoint( db->L-1,field, index,field ); 
			
			DataFieldZeroPoint(field,index);
		}
	}
#endif

	if (index != db->L-1) { /* not last point in list */
		for( f=0; f<db->nfields; f++ ) {
			DataField field = db->field[f];
			
			/* copy then remove */
			DataFieldCopyPoint( db->L-1,field, index,field ); 
			
			//DataFieldZeroPoint(field,index);
		}
	}
	
	/* decrement size */
	/* this will zero out an crap at the end of the list */
	DataBucketRemovePoint(db);
	
}

/* copy x into y */
void DataFieldCopyPoint( const int pid_x, const DataField field_x,
												 const int pid_y, const DataField field_y ) 
{

#ifdef PTATIN_DATAFIELD_POINT_ACCESS_GUARD
	/* check point is valid */
	if( pid_x < 0 ){ printf("ERROR: (IN) index must be >= 0\n"); ERROR(); }
	if( pid_x >= field_x->L ){ printf("ERROR: (IN) index must be < %d\n",field_x->L); ERROR(); }

	if( pid_y < 0 ){ printf("ERROR: (OUT) index must be >= 0\n"); ERROR(); }
	if( pid_y >= field_y->L ){ printf("ERROR: (OUT) index must be < %d\n",field_y->L); ERROR(); }

	if( field_y->atomic_size != field_x->atomic_size ) {
		printf("ERROR: atomic size must match \n"); ERROR();
	}
#endif	
	/*
	memcpy( (void*)((char*)field_y->data + pid_y*field_y->atomic_size), 
					(void*)((char*)field_x->data + pid_x*field_x->atomic_size), 
				  field_x->atomic_size );
	*/
	memcpy(		__DATATFIELD_point_access(field_y->data,pid_y,field_y->atomic_size),
						__DATATFIELD_point_access(field_x->data,pid_x,field_x->atomic_size),
						field_y->atomic_size );
	
}


// zero only the datafield at this point
void DataFieldZeroPoint( const DataField field, const int index ) 
{
#ifdef PTATIN_DATAFIELD_POINT_ACCESS_GUARD
	/* check point is valid */
	if( index < 0 ){ printf("ERROR: index must be >= 0\n"); ERROR(); }
	if( index >= field->L ){ printf("ERROR: index must be < %d\n",field->L); ERROR(); }
#endif	
	
//	memset( (void*)((char*)field->data + index*field->atomic_size), 0, field->atomic_size );
	memset( __DATATFIELD_point_access(field->data,index,field->atomic_size), 0, field->atomic_size );
}

// zero ALL data for this point
void DataBucketZeroPoint( const DataBucket db, const int index ) 
{
	int f;
	
	/* check point is valid */
	if( index < 0 ){ printf("ERROR: index must be >= 0\n"); ERROR(); }
	if( index >= db->allocated ){ printf("ERROR: index must be < %d\n",db->allocated); ERROR(); }
	
	for(f=0;f<db->nfields;f++){
		DataField field = db->field[f];
		
		DataFieldZeroPoint(field,index);
	}	
}

/* increment */
void DataBucketAddPoint( DataBucket db )
{
	DataBucketSetSizes( db, db->L+1, -1 );
}
/* decrement */
void DataBucketRemovePoint( DataBucket db )
{
	DataBucketSetSizes( db, db->L-1, -1 );
}

void DataBucketLoadFromFile(MPI_Comm comm,const char filename[], DataBucketViewType type, DataBucket *db)
{
  switch (type) {
    case DATABUCKET_VIEW_STDOUT:
      printf("ERROR: Cannot load using viewer type = stdout\n");
      MPI_ERROR_CHECK(comm,1);
      break;
      
    case DATABUCKET_VIEW_BINARY:
      printf("ERROR: Cannot load using viewer type = binary\n");
      MPI_ERROR_CHECK(comm,1);
      break;
      
    case DATABUCKET_VIEW_NATIVE:
      DataBucketLoad_NATIVE(comm,filename,db);
      break;
      
    default:
      printf("ERROR: Unknown viewer type\n");
      MPI_ERROR_CHECK(comm,1);
      break;
  }
}

void DataBucketView_STDOUT(MPI_Comm comm,DataBucket db,const char prefix[])
{
  int f;
  long int L,buffer,allocated;
  double memory_usage_total,memory_usage_total_local = 0.0;
  int rank,commsize;
  int ierr;

  ierr = MPI_Comm_size(comm,&commsize);MPI_ERROR_CHECK(comm,ierr);
  ierr = MPI_Comm_rank(comm,&rank);MPI_ERROR_CHECK(comm,ierr);
  
  DataBucketGetGlobalSizes(comm,db,&L,&buffer,&allocated);
  
  for( f=0; f<db->nfields; f++ ) {
    double memory_usage_f = (double)(db->field[f]->atomic_size * db->allocated) * 1.0e-6;
    
    memory_usage_total_local += memory_usage_f;
  }
  ierr = MPI_Allreduce(&memory_usage_total_local,&memory_usage_total,1,MPI_DOUBLE,MPI_SUM,comm);MPI_ERROR_CHECK(comm,ierr);
  
  if (rank == 0) {
    if (prefix) printf("DataBucketView <%s>:\n",prefix);
    else printf("DataBucketView:\n");
    printf("  L                  = %ld \n", L );
    printf("  buffer (max)       = %ld \n", buffer );
    printf("  allocated          = %ld \n", allocated );
    
    printf("  nfields registered = %d \n", db->nfields );
    for( f=0; f<db->nfields; f++ ) {
      double memory_usage_f = (double)(db->field[f]->atomic_size * db->allocated) * 1.0e-6;
      
      printf("    [%3d]: field name  ==>> %30s : Mem. usage = %1.2e (MB) : rank0\n", f, db->field[f]->name, memory_usage_f  );
    }
    
    printf("  Total mem. usage                                                      = %1.2e (MB) : <collective over %d ranks>\n", memory_usage_total, commsize );
  }
}

/*
 cJSON does not support long ints
*/
void DataBucketView_NATIVE(MPI_Comm comm,DataBucket db,const char prefix[])
{
  int commsize,rank;
  int ierr;
  int *pcount = NULL,*bcount = NULL,*acount = NULL,L,buffer,allocated;
  int f;
  char jfilename[2048];
  char fieldfilename[2048];
  FILE *fpbin = NULL;
  
  ierr = MPI_Comm_size(comm,&commsize);MPI_ERROR_CHECK(comm,ierr);
  ierr = MPI_Comm_rank(comm,&rank);MPI_ERROR_CHECK(comm,ierr);

  if (rank == 0) {
    pcount = (int*)malloc(sizeof(int)*commsize);
    bcount = (int*)malloc(sizeof(int)*commsize);
    acount = (int*)malloc(sizeof(int)*commsize);
  }
  
  /* create size array */
  DataBucketGetSizes(db,&L,&buffer,&allocated);
  
  ierr = MPI_Gather(&L,1,MPI_INT,pcount,1,MPI_INT,0,comm);MPI_ERROR_CHECK(comm,ierr);
  ierr = MPI_Gather(&buffer,1,MPI_INT,bcount,1,MPI_INT,0,comm);MPI_ERROR_CHECK(comm,ierr);
  ierr = MPI_Gather(&allocated,1,MPI_INT,acount,1,MPI_INT,0,comm);MPI_ERROR_CHECK(comm,ierr);
  
  sprintf(jfilename,"%s_db.json",prefix);
  sprintf(fieldfilename,"%s_db_data.bin",prefix);
  if (rank == 0) {
    cJSON *jso_file,*jso_db,*jso_part,*fields,*field,*content;
    cJSON *jso;
    
    /* create json meta data file */
    
    jso_file = cJSON_CreateObject();
    
    jso_db = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_file,"DataBucket",jso_db);
    
    jso = cJSON_CreateInt(db->nfields);    cJSON_AddItemToObject(jso_db,"nfields",jso);
    
    fields = cJSON_CreateArray();
    for (f=0; f<db->nfields; f++) {
      
      field = cJSON_CreateObject();
      content = cJSON_CreateString(db->field[f]->name);         cJSON_AddItemToObject(field,"fieldName",content);
      content = cJSON_CreateInt(db->field[f]->atomic_size);     cJSON_AddItemToObject(field,"atomicSize",content);
      content = cJSON_CreateString(db->field[f]->registration_function);     cJSON_AddItemToObject(field,"registrationFunction",content);
      content = cJSON_CreateString("nativeBinary");            cJSON_AddItemToObject(field,"dataFormat",content);
      
      content = cJSON_CreateString(fieldfilename);              cJSON_AddItemToObject(field,"fileName",content);
      
      cJSON_AddItemToArray(fields,field);
    }
    
    // add all fields to data bucket
    cJSON_AddItemToObject(jso_db,"fields",fields);

    jso_part = cJSON_CreateObject();
    cJSON_AddItemToObject(jso_db,"partition",jso_part);

    content = cJSON_CreateInt(commsize);                 cJSON_AddItemToObject(jso_part,"commSize",content);
    content = cJSON_CreateIntArray(pcount,commsize);     cJSON_AddItemToObject(jso_part,"length",content);
    content = cJSON_CreateIntArray(bcount,commsize);     cJSON_AddItemToObject(jso_part,"buffer",content);
    content = cJSON_CreateIntArray(acount,commsize);     cJSON_AddItemToObject(jso_part,"allocated",content);
    
    /* write json meta data file */
    {
      FILE *fp;
      char *jbuff = cJSON_Print(jso_file);
      
      fp = fopen(jfilename,"w");
      fprintf(fp,"%s\n",jbuff);
      fclose(fp);
      /*printf("%s\n",jbuff);*/
      free(jbuff);
    }
    
    cJSON_Delete(jso_file);
  }
  
  /* write raw binary data with the header */
  if (rank == 0) {
    fpbin = fopen(fieldfilename,"w");
  }
  for (f=0; f<db->nfields; f++) {
    /* write only the data being used - we do this so that we can load all the data written in parallel on 1 rank if required */
    ierr = MPIWrite_Blocking(fpbin,db->field[f]->data,db->L,db->field[f]->atomic_size,0,PETSC_FALSE,comm);
  }
  
  if (fpbin)  { fclose(fpbin); }
  if (pcount) { free(pcount); }
  if (bcount) { free(bcount); }
  if (acount) { free(acount); }
}

int _DataBucketRegisterFieldsFromFile_NATIVE(MPI_Comm comm,DataBucket db,cJSON *jso_root)
{
  int ierr,rank;
  int k,nf;
  cJSON *flist,*f_k;
  
  ierr = MPI_Comm_rank(comm,&rank);MPI_ERROR_CHECK(comm,ierr);
  flist = NULL;
  nf = 0;
  if (jso_root) {
    flist = cJSON_GetObjectItem(jso_root,"fields");
    if (!flist) { printf("<error> failed to locate key \"Fields\"\n"); return(1); }
    nf = cJSON_GetArraySize(flist);
    
    f_k = cJSON_GetArrayItemRoot(flist);
    for (k=0; k<nf; k++) {
      int found;
      char *field_name;
      int _atomic_size;
      size_t atomic_size;
      char *registration_function;
      
      cJSON_GetObjectValue_char(f_k,"fieldName",&found,&field_name);
      if (found == cJSON_False) { printf("<error> failed to locate key \"fieldName\"\n");  return(1); }
      
      cJSON_GetObjectValue_int(f_k,"atomicSize",&found,&_atomic_size);
      if (found == cJSON_False) { printf("<error> failed to locate key \"atomicSize\"\n");  return(1); }
      atomic_size = (size_t)_atomic_size;
      
      cJSON_GetObjectValue_char(f_k,"registrationFunction",&found,&registration_function);
      if (found == cJSON_False) { printf("<error> failed to locate key \"registrationFunction\"\n");  return(1); }
      
      _DataBucketRegisterField(db,(const char*)registration_function,(const char*)field_name,atomic_size,NULL);

      f_k = cJSON_GetArrayItemNext(f_k);
    }
  }
  
  /* broadcast from root */
  ierr = MPI_Bcast(&nf,1,MPI_INT,0,comm);MPI_ERROR_CHECK(comm,ierr);
  for (k=0; k<nf; k++) {
    char string_f[2048];
    char string_r[2048];
    int i,asize = 0;
    size_t size;
    
    for (i=0; i<2048; i++) {
      string_f[i] = '\0';
      string_r[i] = '\0';
    }
    
    if (rank == 0) { sprintf(string_f,"%s",db->field[k]->name); }
    ierr = MPI_Bcast(string_f,2048,MPI_CHAR,0,comm);MPI_ERROR_CHECK(comm,ierr);
    
    if (rank == 0) { asize = (int)db->field[k]->atomic_size; }
    ierr = MPI_Bcast(&asize,1,MPI_INT,0,comm);MPI_ERROR_CHECK(comm,ierr);
    size = (size_t)asize;

    //ierr = MPI_Bcast(db->field[k]->registration_function,1,MPI_CHAR,0,comm);MPI_ERROR_CHECK(comm,ierr);
    if (rank == 0) { sprintf(string_r,"%s",db->field[k]->registration_function); }
    ierr = MPI_Bcast(string_r,2048,MPI_CHAR,0,comm);MPI_ERROR_CHECK(comm,ierr);
    
    if (rank != 0) {
      _DataBucketRegisterField(db,(const char*)string_r,(const char*)string_f,size,NULL);
    }

  }
  return(0);
}

int _DataBuckeLoadFieldsFromFile_NATIVE(MPI_Comm comm,DataBucket db,cJSON *jso_root)
{
  int ierr,rank,commsize;
  int k,nf;
  cJSON *flist,*f_k,*part;
  char *filename,*dataformat;
  FILE *fpdata;
  int LBA[3],one2one = 1;
  int L_total,B_max;
  MPI_Status status;
  
  ierr = MPI_Comm_size(comm,&commsize);MPI_ERROR_CHECK(comm,ierr);
  ierr = MPI_Comm_rank(comm,&rank);MPI_ERROR_CHECK(comm,ierr);
  flist = NULL;
  nf = 0;
  if (jso_root) {
    flist = cJSON_GetObjectItem(jso_root,"fields");
    if (!flist) { printf("<error> failed to locate key \"fields\"\n"); return(1); }
    nf = cJSON_GetArraySize(flist);
    
    f_k = cJSON_GetArrayItemRoot(flist);
    for (k=0; k<nf; k++) {
      int found;
      
      cJSON_GetObjectValue_char(f_k,"fileName",&found,&filename);
      if (found == cJSON_False) { printf("<error> failed to locate key \"fileName\"\n"); return(1); }
      
      cJSON_GetObjectValue_char(f_k,"dataFormat",&found,&dataformat);
      if (found == cJSON_False) { printf("<error> failed to locate key \"dataFormat\"\n"); return(1); }
      
      f_k = cJSON_GetArrayItemNext(f_k);
    }
  }

  /* Determine if this is a valid load */
  part = NULL;
  if (jso_root) {
    int found;
    int commsize_file;
    
    part = cJSON_GetObjectItem(jso_root,"partition");
    if (!part) { printf("<error> failed to locate key \"partition\"\n"); return(1); }
    
    cJSON_GetObjectValue_int(part,"commSize",&found,&commsize_file);
    if (found == cJSON_False) { printf("<error> failed to locate key \"commSize\"\n"); return(1); }
    
    if ((commsize == 1) && (commsize != commsize_file)) {
      one2one = 0;
    }
    
    if ((commsize != 1) && (commsize != commsize_file)) {
      one2one = -1;
      printf("[ERROR][_DataBuckeLoadFieldsFromFile_NATIVE] It is only valid to load the data file on the same comm size as that which generated it, or on comm size = 1. Current comm size = %d : Input data generated with comm size = %d.\n",(int)commsize,(int)commsize_file);
    }
  }
  ierr = MPI_Bcast(&one2one,1,MPI_INT,0,comm);MPI_ERROR_CHECK(comm,ierr);
  if (one2one < 0) return(2);
  
  /* post receives - rank 0 will post sends for length and buffer next */
  if ((one2one == 1) && (rank != 0)) {
    ierr = MPI_Recv(LBA,3,MPI_INT,0,rank,comm,&status);MPI_ERROR_CHECK(comm,ierr);
  }

  L_total = 0;
  B_max = 0;
  if (jso_root) {
    int r,found;
    int commsize_file;
    int *L_file,*B_file,nvals;
    
    cJSON_GetObjectValue_int(part,"commSize",&found,&commsize_file);
    L_file = (int*)malloc(sizeof(int)*commsize_file);
    B_file = (int*)malloc(sizeof(int)*commsize_file);
    
    cJSON_GetObjectValue_intarray(part,"length",&found,&nvals,L_file);
    if (found == cJSON_False) { printf("<error> failed to locate key \"length\"\n"); }

    cJSON_GetObjectValue_intarray(part,"buffer",&found,&nvals,B_file);
    if (found == cJSON_False) { printf("<error> failed to locate key \"buffer\"\n"); }

    /* determine sizes if we are loading a parallel data set on commsize = 1 */
    for (r=0; r<commsize_file; r++) {
      L_total += L_file[r];
      if (B_file[r] > B_max) {
        B_max = B_file[r];
      }
    }

    LBA[0] = L_file[0];
    LBA[1] = B_file[0];
    LBA[2] = 0;

    if (one2one == 1) {
      for (r=1; r<commsize_file; r++) {
        LBA[0] = L_file[r];
        LBA[1] = B_file[r];
        LBA[2] = 0;
        ierr = MPI_Send(LBA,3,MPI_INT,r,r,comm);MPI_ERROR_CHECK(comm,ierr);
      }
    }
    
    free(L_file);
    free(B_file);
  }

  /* Special case if we are loading a parallel data set on commsize = 1 */
  if (one2one == 0) {
    LBA[0] = L_total;
    LBA[1] = B_max;
    LBA[2] = 0;
  }
  
  /* allocate space */
  DataBucketSetSizes(db,LBA[0],LBA[1]);
  
  /* broadcast from root */
  fpdata = NULL;
  if (rank == 0) {
    fpdata = fopen(filename,"r");
    if (!fpdata) { printf("<error> failed to open file \"%s\"\n",filename); return(3); }
  }

  /* load data from file */
  for (k=0; k<db->nfields; k++) {
    ierr = MPIRead_Blocking(fpdata,(void**)&db->field[k]->data,db->L,db->field[k]->atomic_size,0,PETSC_FALSE,comm);MPI_ERROR_CHECK(comm,ierr);
  }
  
  if (fpdata) { fclose(fpdata); }
  return(0);
}

void DataBucketLoad_NATIVE(MPI_Comm comm,const char jfilename[],DataBucket *_db)
{
  int ierr,ierr_l,ierr_g,nproc,rank;
  DataBucket db;
  cJSON *jfile = NULL,*jdb = NULL;
  
  ierr = MPI_Comm_size(comm,&nproc);MPI_ERROR_CHECK(comm,ierr);
  ierr = MPI_Comm_rank(comm,&rank);MPI_ERROR_CHECK(comm,ierr);
  
  if (rank == 0) {
    cJSON_FileView(jfilename,&jfile);
    if (!jfile) {
      printf("<error> failed to open JSON file \"%s\"\n",jfilename);
      *_db = NULL;
      return;
    }
    jdb = cJSON_GetObjectItem(jfile,"DataBucket");
  }
  
  DataBucketCreate(&db);
  
  /* load meta data */
  ierr_l = _DataBucketRegisterFieldsFromFile_NATIVE(comm,db,jdb);
  ierr = MPI_Allreduce(&ierr_l,&ierr_g,1,MPI_INT,MPI_MAX,comm);MPI_ERROR_CHECK(comm,ierr);
  if (ierr_g != 0) { MPI_Abort(comm,ierr_g); }

  DataBucketFinalize(db);

  /* load binary data */
  ierr_l = _DataBuckeLoadFieldsFromFile_NATIVE(comm,db,jdb);
  ierr = MPI_Allreduce(&ierr_l,&ierr_g,1,MPI_INT,MPI_MAX,comm);MPI_ERROR_CHECK(comm,ierr);
  if (ierr_g != 0) { MPI_Abort(comm,ierr_g); }
  
  if (jfile) { cJSON_Delete(jfile); }
  
  *_db = db;
}

int _DataBuckeLoadFieldsRedundantFromFile_NATIVE(MPI_Comm comm,DataBucket db,cJSON *jso_root)
{
  int ierr,rank,commsize;
  int k,nf;
  cJSON *flist,*f_k,*part;
  char *filename,*dataformat;
  FILE *fpdata;
  int LBA[3];
  int L_total,B_max;
  
  ierr = MPI_Comm_size(comm,&commsize);MPI_ERROR_CHECK(comm,ierr);
  ierr = MPI_Comm_rank(comm,&rank);MPI_ERROR_CHECK(comm,ierr);
  flist = NULL;
  nf = 0;
  if (jso_root) {
    flist = cJSON_GetObjectItem(jso_root,"fields");
    if (!flist) { printf("<error> failed to locate key \"fields\"\n"); return(1); }
    nf = cJSON_GetArraySize(flist);
    
    f_k = cJSON_GetArrayItemRoot(flist);
    for (k=0; k<nf; k++) {
      int found;
      
      cJSON_GetObjectValue_char(f_k,"fileName",&found,&filename);
      if (found == cJSON_False) { printf("<error> failed to locate key \"fileName\"\n"); return(1); }
      
      cJSON_GetObjectValue_char(f_k,"dataFormat",&found,&dataformat);
      if (found == cJSON_False) { printf("<error> failed to locate key \"dataFormat\"\n"); return(1); }
      
      f_k = cJSON_GetArrayItemNext(f_k);
    }
  }
  
  part = NULL;
  L_total = 0;
  B_max = 0;
  if (jso_root) {
    int r,found;
    int commsize_file;
    int *L_file,*B_file,nvals;

    part = cJSON_GetObjectItem(jso_root,"partition");
    if (!part) { printf("<error> failed to locate key \"partition\"\n"); return(1); }
    
    cJSON_GetObjectValue_int(part,"commSize",&found,&commsize_file);
    L_file = (int*)malloc(sizeof(int)*commsize_file);
    B_file = (int*)malloc(sizeof(int)*commsize_file);
    
    cJSON_GetObjectValue_intarray(part,"length",&found,&nvals,L_file);
    if (found == cJSON_False) { printf("<error> failed to locate key \"length\"\n"); return(1); }
    
    cJSON_GetObjectValue_intarray(part,"buffer",&found,&nvals,B_file);
    if (found == cJSON_False) { printf("<error> failed to locate key \"buffer\"\n"); return(1); }
    
    /* Sum total sizes */
    for (r=0; r<commsize_file; r++) {
      L_total += L_file[r];
      if (B_file[r] > B_max) {
        B_max = B_file[r];
      }
    }
    
    free(L_file);
    free(B_file);
  }
  
  LBA[0] = L_total;
  LBA[1] = B_max;
  LBA[2] = 0;
  ierr = MPI_Bcast(LBA,3,MPI_INT,0,comm);MPI_ERROR_CHECK(comm,ierr);
  
  /* allocate space */
  DataBucketSetSizes(db,LBA[0],LBA[1]);
  
  /* broadcast from root */
  fpdata = NULL;
  if (rank == 0) {
    fpdata = fopen(filename,"r");
    if (!fpdata) { printf("<error> failed to open file \"%s\"\n",filename); return(3); }
  }
  
  /* load data from file - ensure rank 0 reads everything - post read we broadcast */
  LBA[0] = L_total;
  LBA[1] = B_max;
  LBA[2] = 0;
  for (k=0; k<db->nfields; k++) {
    ierr = MPIRead_Blocking(fpdata,(void**)&db->field[k]->data,LBA[0],db->field[k]->atomic_size,0,PETSC_FALSE,comm);MPI_ERROR_CHECK(comm,ierr);

    ierr = MPI_Bcast(db->field[k]->data,db->L*db->field[k]->atomic_size,MPI_BYTE,0,comm);MPI_ERROR_CHECK(comm,ierr);
  }
  
  if (fpdata) { fclose(fpdata); }
  return(0);
}

void DataBucketLoadRedundant_NATIVE(MPI_Comm comm,const char jfilename[],DataBucket *_db)
{
  int ierr,ierr_l,ierr_g,nproc,rank;
  DataBucket db;
  cJSON *jfile = NULL,*jdb = NULL;
  
  ierr = MPI_Comm_size(comm,&nproc);MPI_ERROR_CHECK(comm,ierr);
  ierr = MPI_Comm_rank(comm,&rank);MPI_ERROR_CHECK(comm,ierr);
  
  if (rank == 0) {
    cJSON_FileView(jfilename,&jfile);
    if (!jfile) {
      printf("<error> failed to open JSON file \"%s\"\n",jfilename);
      *_db = NULL;
      return;
    }
    jdb = cJSON_GetObjectItem(jfile,"DataBucket");
  }
  
  DataBucketCreate(&db);
  
  /* load meta data */
  ierr_l = _DataBucketRegisterFieldsFromFile_NATIVE(comm,db,jdb);
  ierr = MPI_Allreduce(&ierr_l,&ierr_g,1,MPI_INT,MPI_MAX,comm);MPI_ERROR_CHECK(comm,ierr);
  if (ierr_g != 0) { MPI_Abort(comm,ierr_g); }
  
  DataBucketFinalize(db);
  
  /* load binary data */
  ierr_l = _DataBuckeLoadFieldsRedundantFromFile_NATIVE(comm,db,jdb);
  ierr = MPI_Allreduce(&ierr_l,&ierr_g,1,MPI_INT,MPI_MAX,comm);MPI_ERROR_CHECK(comm,ierr);
  if (ierr_g != 0) { MPI_Abort(comm,ierr_g); }
  
  if (jfile) { cJSON_Delete(jfile); }
  
  *_db = db;
}

void DataBucketLoadRedundantFromFile(MPI_Comm comm,const char filename[], DataBucketViewType type, DataBucket *db)
{
  switch (type) {
    case DATABUCKET_VIEW_STDOUT:
      printf("ERROR: Cannot load (redundant) using viewer type = stdout\n");
      MPI_ERROR_CHECK(comm,1);
      break;
      
    case DATABUCKET_VIEW_BINARY:
      printf("ERROR: Cannot load (redundant) using viewer type = binary\n");
      MPI_ERROR_CHECK(comm,1);
      break;
      
    case DATABUCKET_VIEW_NATIVE:
      DataBucketLoadRedundant_NATIVE(comm,filename,db);
      break;
      
    default:
      printf("ERROR: Unknown viewer type\n");
      MPI_ERROR_CHECK(comm,1);
      break;
  }
}

void DataBucketView(MPI_Comm comm,DataBucket db,const char prefix[],DataBucketViewType type)
{
  switch (type) {
    case DATABUCKET_VIEW_STDOUT:
      DataBucketView_STDOUT(comm,db,prefix);
      break;
      
    case DATABUCKET_VIEW_BINARY:
      printf("ERROR: Binary viewer is not implemented\n");
      MPI_ERROR_CHECK(comm,1);
      break;

    case DATABUCKET_VIEW_NATIVE:
      DataBucketView_NATIVE(comm,db,prefix);
      break;
      
    default:
      printf("ERROR: Unknown viewer type\n");
      MPI_ERROR_CHECK(comm,1);
      break;
  }
}

void DataBucketDuplicateFields(DataBucket dbA,DataBucket *dbB)
{
	DataBucket db2;
	int f;
	
	DataBucketCreate(&db2);
	
	/* copy contents from dbA into db2 */
	for (f=0; f<dbA->nfields; f++) {
		DataField field;
		size_t    atomic_size;
		char      *name;
		
		field = dbA->field[f];
		
		atomic_size = field->atomic_size;
		name        = field->name;
		
		DataBucketRegisterField(db2,name,atomic_size,NULL);
	}
	DataBucketFinalize(db2);
	DataBucketSetInitialSizes(db2,0,1000);
	
	/* set pointer */
	*dbB = db2;
}

/*
 Insert points from db2 into db1
 db1 <<== db2
 */
void DataBucketInsertValues(DataBucket db1,DataBucket db2)
{
	int n_mp_points1,n_mp_points2;
	int n_mp_points1_new,p;
	
	DataBucketGetSizes(db1,&n_mp_points1,0,0);
	DataBucketGetSizes(db2,&n_mp_points2,0,0);
	
	n_mp_points1_new = n_mp_points1 + n_mp_points2;
	DataBucketSetSizes(db1,n_mp_points1_new,-1);
	
	for (p=0; p<n_mp_points2; p++) {
		// db1 <<== db2 //
		DataBucketCopyPoint( db2,p, db1,(n_mp_points1 + p) );
	}
}

/* helpers for parallel send/recv */
void DataBucketCreatePackedArray(DataBucket db,size_t *bytes,void **buf)
{
    int       f;
    size_t    sizeof_marker_contents;
    void      *buffer;
    
    sizeof_marker_contents = 0;
    for (f=0; f<db->nfields; f++) {
        DataField df = db->field[f];
        
        sizeof_marker_contents += df->atomic_size;
    }
    
    buffer = malloc(sizeof_marker_contents);
    memset(buffer,0,sizeof_marker_contents);
    
    if (bytes) { *bytes = sizeof_marker_contents; }
    if (buf)   { *buf   = buffer; }
}

void DataBucketDestroyPackedArray(DataBucket db,void **buf)
{
    if (buf) {
        free(*buf);
        *buf = NULL;
    }
}

void DataBucketFillPackedArray(DataBucket db,const int index,void *buf)
{
    int    f;
    void   *data,*data_p;
    size_t asize,offset;
    
    offset = 0;
    for (f=0; f<db->nfields; f++) {
        DataField df = db->field[f];
        
        asize = df->atomic_size;
        
        data = (void*)( df->data );
        data_p = (void*)( (char*)data + index*asize );
        
        memcpy( (void*)((char*)buf + offset),  data_p,  asize);
        offset = offset + asize;
    }
}

void DataBucketInsertPackedArray(DataBucket db,const int idx,void *data)
{
    int f;
    void *data_p;
    size_t offset;
    
    offset = 0;
    for (f=0; f<db->nfields; f++) {
        DataField df = db->field[f];
        
        data_p = (void*)( (char*)data + offset );
        
        DataFieldInsertPoint(df, idx, (void*)data_p );
        offset = offset + df->atomic_size;
    }
}
