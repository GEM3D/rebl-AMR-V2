#include <stdlib.h>
#include "typedefs.h"
#include "tree.h"
#include "partition.h"
#define TOL "1.1"


// refer to manual for ZOLTAN_ID_TYPE and ZOLTAN_ID_PTR
// note that ZOLTAN_ID_PTR=*ZOLTAN_ID_TYPE

static int get_number_of_objects( void *data, int *ierr );

static int get_num_geometry( void *data, int *ierr );

static void get_object_list( void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim,
                             float *obj_wgts, int *ierr );

static void get_geometry_list( void *data, int sizeGID, int sizeLID, int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                               int num_dim, double *geom_vec, int *ierr );


Partition::Partition(int argcs,char *pArgs[],int meth,int sze)
{

myMesh.c=nullptr;
myMesh.w=nullptr;
myMesh.myGlobalIDs=nullptr;
/*
int initFlag=0;

    if ( MPI_Initialized( &initFlag ) != MPI_SUCCESS  && Com.myrank==0 )
    {
        cout << " Exit Code : " << ReblAmrGetErrorEnum( MPI_INIT_CHECK_FAIL ) << endl;
        exit(1 );
    }

    if ( initFlag == 0 )
    {
        if ( MPI_Init( &argcs, &pArgs ) != MPI_SUCCESS  && Com.myrank==0  )
        {
            cout << " Exit Code : " << ReblAmrGetErrorEnum( MPI_INIT_FAIL ) << endl;
            exit( 1 );
        }
    }
  MPIStartUp();
*/
 float ver;
 int rc;
//if(ZOLTAN_ON==1)
{
  rc = Zoltan_Initialize( argcs, pArgs, &ver );

  zz = Zoltan_Create( MPI_COMM_SELF );
 
 if ( rc != ZOLTAN_OK )
        {
            printf( RED "Problem in Initialzing Zoltan\n" RESET );
            exit( 0 );
        }

 method=meth;
 size=sze;

 weight=new real[size];
}
}


void Partition::construct(int argcs,char *pArgs[],int meth,int sze)
{

myMesh.c=nullptr;
myMesh.w=nullptr;
myMesh.myGlobalIDs=nullptr;
/*
int initFlag=0;

    if ( MPI_Initialized( &initFlag ) != MPI_SUCCESS  && Com.myrank==0 )
    {
        cout << " Exit Code : " << ReblAmrGetErrorEnum( MPI_INIT_CHECK_FAIL ) << endl;
        exit(1 );
    }

    if ( initFlag == 0 )
    {
        if ( MPI_Init( &argcs, &pArgs ) != MPI_SUCCESS  && Com.myrank==0  )
        {
            cout << " Exit Code : " << ReblAmrGetErrorEnum( MPI_INIT_FAIL ) << endl;
            exit( 1 );
        }
    }

 MPIStartUp();
*/ 

 float ver;
 int rc;
//if(ZOLTAN_ON==1)
{
 rc = Zoltan_Initialize( argcs, pArgs, &ver );

  zz = Zoltan_Create( MPI_COMM_SELF );

 if ( rc != ZOLTAN_OK )
        {
            printf( RED "Problem in Initialzing Zoltan\n" RESET );
            exit( 0 );
        }

 method=meth;
 size=sze;

 //cout<<" mesh size for zoltan  "<<size<<endl;

 weight=new real[size];
}
}




Partition::~Partition()
{
/*
if(zoltan_out!=NULL)
{
    
        Zoltan_LB_Free_Part( &(zoltan_out->importGlobalGids),&(zoltan_out->importLocalGids),&(zoltan_out->importProcs), &(zoltan_out->importToPart) );

        Zoltan_LB_Free_Part( &(zoltan_out->exportGlobalGids), &(zoltan_out->exportLocalGids), &(zoltan_out->exportProcs), &(zoltan_out->exportToPart) );

        Zoltan_Destroy( &zz );
}
*/
if(ZOLTAN_ON==1)
{
if(myMesh.c!=nullptr)
{
    delete[] myMesh.c;
}
if(myMesh.w!=nullptr)
{
    delete[] myMesh.w;
}
if(myMesh.myGlobalIDs!=nullptr)
{
    delete[] myMesh.myGlobalIDs;
}


if(weight!=nullptr)
{
   delete[] weight;
}
}
/*
MPI_Comm_free(&Com.mpicom);
*/
}



#if(0)
bool Partition::zoltanGeometricPartitioner( const uint size, const uint ncube_total, const uint offset )
{
    float ver;
    int rc;

/* General parameters */

    zoltanSetParams();
#if ( 0 )

    // this is added for serial version, it is to
    Zoltan_Set_Param( zz, "NUM_LOCAL_PARTS", "3" );

    Zoltan_Set_Param( zz, "IMBALANCE_TOL", TOL );

    Zoltan_Set_Param( zz, "DEBUG_LEVEL", "4" );
    // Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");

    Zoltan_Set_Param( zz, "LB_APPROACH", "PARTITION" );

    if ( method == 1 )
    {
        Zoltan_Set_Param( zz, "LB_METHOD", "HSFC" );
    }
    else if ( method == 2 )
    {
        Zoltan_Set_Param( zz, "LB_METHOD", "RCB" );
    }
    else if ( method == 3 )
    {
        Zoltan_Set_Param( zz, "LB_METHOD", "RIB" );
    }
    else
    {
        printf( "Specified Method Needs a Graph\n" );
        printf( RED "ERROR:Zoltan error \n" RESET );
        MPI_Finalize();
        exit( 0 );
    }
    //  Zoltan_Set_Param(zz, "LB_METHOD", "BLOCK");
    Zoltan_Set_Param( zz, "NUM_GID_ENTRIES", "1" );
    Zoltan_Set_Param( zz, "NUM_LID_ENTRIES", "1" );
    // Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
    Zoltan_Set_Param( zz, "RETURN_LISTS", "ALL" );

    // need to tell zoltan if you have weights associated with the edge or vertex, setting this to "0" means no weighting
    Zoltan_Set_Param( zz, "OBJ_WEIGHT_DIM", "1" );
#endif
// printf( "size=%d\n", size );

    MeshData myMesh;
    //	myMesh=(MeshData*)malloc(1.0*sizeof(MESH_DATA));

    myMesh.numGlobalPoints = ncube_total;
    myMesh.numMyPoints = size;
    cout<<"size= "<<size<<endl;
    myMesh.myGlobalIDs = nullptr;
    myMesh.myGlobalIDs = (ZOLTAN_ID_PTR)malloc( size * sizeof( ZOLTAN_ID_TYPE ) );

    //	myMesh.myGlobalIDs=new ZOLTAN_ID_TYPE[size];

    for ( unsigned int i = 0; i < size; i++ )
    {
        myMesh.myGlobalIDs[i] = offset + i;
    }

    myMesh.c = nullptr;
    myMesh.c = new real[3 * myMesh.numMyPoints];
    //  myMesh.c=(real*)malloc(3*myMesh.numMyPoints*sizeof(real));

    myMesh.w = nullptr;
    myMesh.w = new real[myMesh.numMyPoints];
    //   myMesh.w=(real*)malloc(myMesh.numMyPoints*sizeof(real));

    // due to modification in structure

    for ( unsigned int i = 0; i < myMesh.numMyPoints; i++ )
    {
        /*
        myMesh.c[3*i]=0.5*(cube[i].xyz[0]+cube[i].xyz[1]);
        myMesh.c[3*i+1]=0.5*(cube[i].xyz[2]+cube[i].xyz[3]);
        myMesh.c[3*i+2]=0.5*(cube[i].xyz[4]+cube[i].xyz[5]);
    */
        myMesh.c[3 * i] = XYZ.at( i ).x;
        myMesh.c[3 * i + 1] = XYZ.at( i ).y;
        myMesh.c[3 * i + 2] = XYZ.at( i ).z;

        myMesh.w[i] = weight[i];
        	printf(" %lf\n",myMesh.w[i]);
    }

    // printf("my_rank =%d number of my points %d glob_id global_id %d %d
    // \n",my_rank,myMesh->numMyPoints,myMesh->myGlobalIDs[0],myMesh->myGlobalIDs[1]);

    /* Query functions, to provide geometry to Zoltan */

    Zoltan_Set_Num_Obj_Fn( zz, get_number_of_objects, &myMesh );

    Zoltan_Set_Obj_List_Fn( zz, get_object_list, &myMesh );

    Zoltan_Set_Num_Geom_Fn( zz, get_num_geometry, &myMesh );

    Zoltan_Set_Geom_Multi_Fn( zz, get_geometry_list, &myMesh );

#if ( 1 )
    rc = Zoltan_LB_Partition( zz,                                /* input (all remaining fields are output) */
                              &( zoltan_out->changes ),          /* 1 if partitioning was changed, 0 otherwise */
                              &( zoltan_out->numGidEntries ),    /* Number of integers used for a global ID */
                              &( zoltan_out->numLidEntries ),    /* Number of integers used for a local ID */
                              &( zoltan_out->numImport ),        /* Number of vertices to be sent to me */
                              &( zoltan_out->importGlobalGids ), /* Global IDs of vertices to be sent to me */
                              &( zoltan_out->importLocalGids ),  /* Local IDs of vertices to be sent to me */
                              &( zoltan_out->importProcs ),      /* Process rank for source of each incoming vertex */
                              &( zoltan_out->importToPart ),     /* New partition for each incoming vertex */
                              &( zoltan_out->numExport ),        /* Number of vertices I must send to other processes*/
                              &( zoltan_out->exportGlobalGids ), /* Global IDs of the vertices I must send */
                              &( zoltan_out->exportLocalGids ),  /* Local IDs of the vertices I must send */
                              &( zoltan_out->exportProcs ),      /* Process to which I send each of the vertices */
                              &( zoltan_out->exportToPart ) );   /* Partition to which each vertex will belong */

    bool success = true;
    if ( rc != ZOLTAN_OK )
    {
        std::cout << RED "Warning: Zoltan Partition Rejected" RESET << endl;
        success = false;
        //        MPI_Finalize();
        //        exit( 0 );
    }

#endif
    // printf( "my_rank %d changes %d numExport %d numImport %d\n",my_rank,zoltan_out->changes,zoltan_out->numExport,zoltan_out->numImport);
/*
  deleting upon destruction of the class
    delete[] myMesh.c;
    delete[] myMesh.w;
    delete[] myMesh.myGlobalIDs;
*/
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//                               Serial Partitioning
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
bool Partition::zoltanGeometricPartitionerSerial( const uint size, const uint ncube_total, const uint offset, int comsize, Zoltan_Out *zoltan_out )
{
    float ver;
    int rc;


#if ( 1 )

#if(1)
    char csize[16];

    sprintf( csize, "%d", comsize );

   // std::cout<<"????????????????????????????????   "<<csize<<endl;

    //   printf("%d",comsize);
    // this is added for serial version, it is to
    Zoltan_Set_Param( zz, "NUM_LOCAL_PARTS", csize );

    Zoltan_Set_Param( zz, "IMBALANCE_TOL", TOL );

    Zoltan_Set_Param( zz, "DEBUG_LEVEL", "0");
    // Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");

    Zoltan_Set_Param( zz, "LB_APPROACH", "PARTITION" );

    if ( method == 1 )
    {
        Zoltan_Set_Param( zz, "LB_METHOD", "HSFC" );
    }
    else if ( method == 2 )
    {
        Zoltan_Set_Param( zz, "LB_METHOD", "RCB" );
    }
    else if ( method == 3 )
    {
        Zoltan_Set_Param( zz, "LB_METHOD", "RIB" );
    }
    else
    {
        printf( "Specified Method Needs a Graph\n" );
        printf( RED "ERROR:Zoltan error \n" RESET );
        MPI_Finalize();
        exit( 0 );
    }
    //  Zoltan_Set_Param(zz, "LB_METHOD", "BLOCK");
    Zoltan_Set_Param( zz, "NUM_GID_ENTRIES", "1" );
    Zoltan_Set_Param( zz, "NUM_LID_ENTRIES", "1" );
    // Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
    Zoltan_Set_Param( zz, "RETURN_LISTS", "ALL" );
    // need to tell zoltan if you have weights associated with the edge or vertex, setting this to "0" means no weighting
    Zoltan_Set_Param( zz, "OBJ_WEIGHT_DIM", "1" );
#endif
    printf( "size=%d\n", size );

//    MeshData myMesh;
    //	myMesh=(MeshData*)malloc(1.0*sizeof(MESH_DATA));

    cout<<"ccube_total "<<ncube_total<<endl;
    myMesh.numGlobalPoints = ncube_total;
    myMesh.numMyPoints = size;
    myMesh.myGlobalIDs = nullptr;
    //myMesh.myGlobalIDs = (ZOLTAN_ID_PTR)malloc( size * sizeof( ZOLTAN_ID_TYPE ) );
    myMesh.myGlobalIDs = new ZOLTAN_ID_TYPE [ size] ;

    //	myMesh.myGlobalIDs=new ZOLTAN_ID_TYPE[size];

    for ( unsigned int i = 0; i < size; i++ )
    {
        //myMesh.myGlobalIDs[i] = offset + i;
        myMesh.myGlobalIDs[i] =  i;
      // cout<<" ids "<<myMesh.myGlobalIDs[i]<<endl;
    }

    myMesh.c = nullptr;
    myMesh.c = new real[3 * myMesh.numMyPoints];
    //  myMesh.c=(real*)malloc(3*myMesh.numMyPoints*sizeof(real));

    myMesh.w = nullptr;
    myMesh.w = new real[myMesh.numMyPoints];
    //   myMesh.w=(real*)malloc(myMesh.numMyPoints*sizeof(real));

    // due to modification in structure

    for ( unsigned int i = 0; i < myMesh.numMyPoints; i++ )
    {
        /*
        myMesh.c[3*i]=0.5*(cube[i].xyz[0]+cube[i].xyz[1]);
        myMesh.c[3*i+1]=0.5*(cube[i].xyz[2]+cube[i].xyz[3]);
        myMesh.c[3*i+2]=0.5*(cube[i].xyz[4]+cube[i].xyz[5]);
    */
        myMesh.c[3 * i] = XYZ.at( i ).x;
        myMesh.c[3 * i + 1] = XYZ.at( i ).y;
        myMesh.c[3 * i + 2] = XYZ.at( i ).z;

        myMesh.w[i] = weight[i];
//        	printf("%d %lf %lf %lf %lf\n", i  ,myMesh.w[i],myMesh.c[3*i],myMesh.c[3*i+1],myMesh.c[3*i+2]);
    }

//     printf("my_rank =%d number of my points %d glob_id global_id %d %d \n",Com.myrank,myMesh->numMyPoints,myMesh->myGlobalIDs[0],myMesh->myGlobalIDs[1]);

    /* Query functions, to provide geometry to Zoltan */

    Zoltan_Set_Num_Obj_Fn( zz, get_number_of_objects, &myMesh );

    Zoltan_Set_Obj_List_Fn( zz, get_object_list, &myMesh );

    Zoltan_Set_Num_Geom_Fn( zz, get_num_geometry, &myMesh );

    Zoltan_Set_Geom_Multi_Fn( zz, get_geometry_list, &myMesh );

//    Zoltan_Out zoltan_out;
    // double t1,t2;
    // t1=MPI_Wtime();
#if ( 1 )
    rc = Zoltan_LB_Partition( zz,                                /* input (all remaining fields are output) */
                              &( zoltan_out->changes ),          /* 1 if partitioning was changed, 0 otherwise */
                              &( zoltan_out->numGidEntries ),    /* Number of integers used for a global ID */
                              &( zoltan_out->numLidEntries ),    /* Number of integers used for a local ID */
                              &( zoltan_out->numImport ),        /* Number of vertices to be sent to me */
                              &( zoltan_out->importGlobalGids ), /* Global IDs of vertices to be sent to me */
                              &( zoltan_out->importLocalGids ),  /* Local IDs of vertices to be sent to me */
                              &( zoltan_out->importProcs ),      /* Process rank for source of each incoming vertex */
                              &( zoltan_out->importToPart ),     /* New partition for each incoming vertex */
                              &( zoltan_out->numExport ),        /* Number of vertices I must send to other processes*/
                              &( zoltan_out->exportGlobalGids ), /* Global IDs of the vertices I must send */
                              &( zoltan_out->exportLocalGids ),  /* Local IDs of the vertices I must send */
                              &( zoltan_out->exportProcs ),      /* Process to which I send each of the vertices */
                              &( zoltan_out->exportToPart ) );   /* Partition to which each vertex will belong */
// t2=MPI_Wtime();
// cout<<RED<<"partition time "<<t2-t1<<RESET<<endl;

    bool success = true;
    if ( rc != ZOLTAN_OK )
    {
        std::cout << RED "Warning: Zoltan Partition Rejected" RESET << endl;
        success = false;
        //        MPI_Finalize();
        //        exit( 0 );
    }


#endif
 //   cout<<"zoltn Part Done"<<endl;
#endif
    // printf( "my_rank %d changes %d numExport %d numImport %d\n",my_rank,zoltan_out->changes,zoltan_out->numExport,zoltan_out->numImport);
/*
    delete[] myMesh.c;
    delete[] myMesh.w;
    delete[] myMesh.myGlobalIDs;
*/
return true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//  Interface functions for HSFC
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// for passing void and type conversion inside a function and how it works see The C++ Programming Language by Bjarne S. section 7.2
static int get_number_of_objects( void *data, int *ierr )
{
    MeshData *mesh = (MeshData *)data;
    *ierr = ZOLTAN_OK;
    printf("number of my points %d\n",mesh->numMyPoints);
    return mesh->numMyPoints;
}

static void get_object_list( void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim,
                             float *obj_wgts, int *ierr )
{

    MeshData *mesh = (MeshData *)data;
    *ierr = ZOLTAN_OK;

    // printf("?????????????????\n \n %d \n",mesh->numMyPoints);

    for ( unsigned int i = 0; i < mesh->numMyPoints; i++ )
    {
        globalID[i] = mesh->myGlobalIDs[i];
        localID[i] = i;
        obj_wgts[i] = mesh->w[i];
        //  mesh->w[i];
        // printf("local global %d %d\n",localID[i] ,globalID[i]);
    }
}

//======================================================================

static int get_num_geometry( void *data, int *ierr )
{
    *ierr = ZOLTAN_OK;
    return 3;
}

//======================================================================

static void get_geometry_list( void *data, int sizeGID, int sizeLID, int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                               int num_dim, double *geom_vec, int *ierr )
{

    int i;

    MeshData *mesh = (MeshData *)data;

    /*
      if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 2)){
        *ierr = ZOLTAN_FATAL;
        return;
      }
    */

    *ierr = ZOLTAN_OK;

    for ( i = 0; i < num_obj; i++ )
    {
        geom_vec[3 * i] = (double)mesh->c[3 * i];
        geom_vec[3 * i + 1] = (double)mesh->c[3 * i + 1];
        geom_vec[3 * i + 2] = (double)mesh->c[3 * i + 2];
    }
}
#if(1)
void Partition::zoltanSetParams()
{

    char csize[16];
    sprintf( csize, "%d", Com.comsize );


    Zoltan_Set_Param( zz, "NUM_LOCAL_PARTS", csize );

    Zoltan_Set_Param( zz, "IMBALANCE_TOL", TOL );

    Zoltan_Set_Param( zz, "DEBUG_LEVEL", "0" );

    Zoltan_Set_Param( zz, "LB_APPROACH", "PARTITION" );

    if ( method == 1 )
    {
        Zoltan_Set_Param( zz, "LB_METHOD", "HSFC" );
    }
    else if ( method == 2 )
    {
        Zoltan_Set_Param( zz, "LB_METHOD", "RCB" );
    }
    else if ( method == 3 )
    {
        Zoltan_Set_Param( zz, "LB_METHOD", "RIB" );
    }
    else
    {
        printf( "Specified Method Needs a Graph\n" );
        printf( RED "ERROR:Zoltan error \n" RESET );
        MPI_Finalize();
        exit( 0 );
    }
    //  Zoltan_Set_Param(zz, "LB_METHOD", "BLOCK");
    Zoltan_Set_Param( zz, "NUM_GID_ENTRIES", "1" );
    Zoltan_Set_Param( zz, "NUM_LID_ENTRIES", "1" );
    // Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
    Zoltan_Set_Param( zz, "RETURN_LISTS", "ALL" );

    // need to tell zoltan if you have weights associated with the edge or vertex, setting this to "0" means no weighting
    Zoltan_Set_Param( zz, "OBJ_WEIGHT_DIM", "1" );

}
#endif
void Partition::MPIStartUp()
{
  if ( MPI_SUCCESS != MPI_Comm_dup( MPI_COMM_WORLD, &Com.mpicom ) )
    {
        cout << BLUE << " Rank(" << Com.myrank << ") > Exit Code : " << MPI_DUP_FAIL << RESET << endl;
        cout << BLUE << ReblAmrGetErrorEnum( MPI_DUP_FAIL ) << RESET << endl;
        exit( 1 );
    }

    if ( MPI_SUCCESS != MPI_Comm_rank( Com.mpicom, &Com.myrank ) )
    {
        cout << BLUE << " Rank(" << Com.myrank << ") > Exit Code : " << MPI_GET_RANK_FAIL << RESET << endl;
        cout << BLUE << ReblAmrGetErrorEnum( MPI_GET_RANK_FAIL ) << RESET << endl;
        exit( 1 );
    }

    if ( MPI_SUCCESS != MPI_Comm_size( Com.mpicom, &Com.comsize ) )
    {
        cout << BLUE << " Rank(" << Com.myrank << ") > Exit Code : " << MPI_COMSIZE_FAIL << RESET << endl;
        cout << BLUE << ReblAmrGetErrorEnum( MPI_COMSIZE_FAIL ) << RESET << endl;
        exit( 1 );
    }
}

#if(0)
void Partition::partition(double *IDi, int totalvalue)
{
      float ver;
        int rc;
//        totalvalue = proc.size();

        rc = Zoltan_Initialize( argcs, pArgs, &ver );

        if ( rc != ZOLTAN_OK )
        {
            printf( RED "Problem in Initialzing Zoltan\n" RESET );
            MPI_Finalize();
            exit( 0 );
        }

        struct Zoltan_Struct *zz = NULL;
        //        zz = Zoltan_Create( MPI_COMM_WORLD );
        //        MPI_COMM_SELF is used to be able to utilize parallel function in serial
        zz = Zoltan_Create( MPI_COMM_SELF );
        // I define this sturcuture to mask sending all these arrays separately
        // and send only one structure instead

        Zoltan_Out zoltan_out;

        if ( ZOLTAN_GEOMETRIC_PARTITION )
        {
            zoltanGeometricPartitionerSerial( myvalue, totalvalue, offset, 1, zz, XYZ, weight, &zoltan_out, comsize );
        }
        //    if ( my_rank == 0 )
        // non scalable part for zoltan partitioning

        uint *ID = new unsigned int[proc.size()];

        for ( uint i = 0; i < proc.size(); i++ )
        {
            ID[i] = 0;
        }
        // elements to be exported are given,

        for ( int i = 0; i < zoltan_out.numExport; i++ )
        {
            ID[zoltan_out.exportLocalGids[i]] = zoltan_out.exportToPart[i];
        }

       Zoltan_LB_Free_Part( &zoltan_out.importGlobalGids, &zoltan_out.importLocalGids, &zoltan_out.importProcs, &zoltan_out.importToPart );

        Zoltan_LB_Free_Part( &zoltan_out.exportGlobalGids, &zoltan_out.exportLocalGids, &zoltan_out.exportProcs, &zoltan_out.exportToPart );

        Zoltan_Destroy( &zz );
}

#endif
