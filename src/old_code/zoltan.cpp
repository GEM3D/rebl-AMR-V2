#include <stdlib.h>
#include "typedefs.h"
#include "tree.h"
#define TOL "1.1"

/*

   Two Different Data Structures are Used for Graph and Mesh Data
         Only hypergraph can handle nonsymmetric graphs

 ----------------------------------------------------------------------
                  GEOMETRIC (Coordinate-Based) of Zoltan
  ----------------------------------------------------------------------

   set method=1 for HSFC: Hilbert Space Filling Curve
   set method=2 for  RCB: Recursive Coordinate Bisection
   set method=3 for  RIB: Recursive Inertial Bisection

    need to use coordinate transformation to impose the principal peorcessor topology
    direction, otherwise it this routine starts from z=0 plane, this
    might not necessarily be the best option as it might require too many
    elements to be moved unnecessarily, just rotate the xyz directions if you
    have more processor in x and y directions. e.g. 3 1 1 is needs transformation
    but 1 1 3 is ok, default is that nz>=ny and nz >=nx
 -----------------------------------------------------------------------
                       HyperGraph of Zoltan
 -----------------------------------------------------------------------

   For HyperGraph in Zoltan
   set method=0
 -----------------------------------------------------------------------
                            ParMETIS
 -----------------------------------------------------------------------

   set method=1 for PARTKWAY
   set method=2 for PARTGEOMKWAY (Hybrid Method uses geom+graph)
   set method=3 for ADEAPTIVEREPART
   set method=4 for PARTGEOM
   the *c in GraphData is only assigned when needed
   Default method is ADAPTIVEREPART as set per Zoltan

   #define TOL "1.1" defines the imbalance tolerance for
   Graph and Hypergraph  methods

 -----------------------------------------------------------------------
*/
/*!
 * \struct MeshData
 * \brief struct to supply information required by zoltan for geometric partitioners
 *
 */

typedef struct
{
    ZOLTAN_ID_TYPE numGlobalPoints;
    ZOLTAN_ID_TYPE numMyPoints;
    ZOLTAN_ID_PTR myGlobalIDs;
    real *c;
    real *w;
} MeshData;

/*!
 * \struct GraphData
 * \brief struct to supply information required by zoltan for geometric partitioners
 *
 */

typedef struct
{
    ZOLTAN_ID_TYPE numMyVertices; /* total vertices in in my partition */
                                  // ZOLTAN_ID_TYPE numAllNbors;   /* total number of neighbors of my vertices */
    ZOLTAN_ID_PTR vertexGID;      /* global ID of each of my vertices */
    ZOLTAN_ID_PTR nbrIndex;       /* nborIndex[i] is location of start of neighbors for vertex i */
    ZOLTAN_ID_PTR nbrGID;         /* nborGIDs[nborIndex[i]] is first neighbor of vertex i */
    int *nbrProc;                 /* process owning each nbor in nborGID */
    float *c;                     /* this one is for the center of the elements*/
} GraphData;

// refer to manual for ZOLTAN_ID_TYPE and ZOLTAN_ID_PTR
// note that ZOLTAN_ID_PTR=*ZOLTAN_ID_TYPE

static int get_number_of_objects( void *data, int *ierr );

static int get_num_geometry( void *data, int *ierr );

static void get_object_list( void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim,
                             float *obj_wgts, int *ierr );

static void get_geometry_list( void *data, int sizeGID, int sizeLID, int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                               int num_dim, double *geom_vec, int *ierr );

bool zoltanGeometricPartitioner( const uint size, const uint ncube_total, const uint offset, const int method, struct Zoltan_Struct *zz,
                                 const Center_coords &XYZ, real *weight, Zoltan_Out *zoltan_out )
{
    float ver;
    int rc;

/* General parameters */

#if ( 1 )

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

// printf( "size=%d\n", size );

#if ( 1 )
    MeshData myMesh;
    //	myMesh=(MeshData*)malloc(1.0*sizeof(MESH_DATA));

    myMesh.numGlobalPoints = ncube_total;
    myMesh.numMyPoints = size;
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
        //	printf(" %lf\n",myMesh.w[i]);
    }

    // printf("my_rank =%d number of my points %d glob_id global_id %d %d
    // \n",my_rank,myMesh->numMyPoints,myMesh->myGlobalIDs[0],myMesh->myGlobalIDs[1]);

    /* Query functions, to provide geometry to Zoltan */

    Zoltan_Set_Num_Obj_Fn( zz, get_number_of_objects, &myMesh );

    Zoltan_Set_Obj_List_Fn( zz, get_object_list, &myMesh );

    Zoltan_Set_Num_Geom_Fn( zz, get_num_geometry, &myMesh );

    Zoltan_Set_Geom_Multi_Fn( zz, get_geometry_list, &myMesh );

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
#endif

#endif
    bool success = true;
    if ( rc != ZOLTAN_OK )
    {
        std::cout << RED "Warning: Zoltan Partition Rejected" RESET << endl;
        success = false;
        //        MPI_Finalize();
        //        exit( 0 );
    }

    // printf( "my_rank %d changes %d numExport %d numImport %d\n",my_rank,zoltan_out->changes,zoltan_out->numExport,zoltan_out->numImport);

    delete[] myMesh.c;
    delete[] myMesh.w;
    delete[] myMesh.myGlobalIDs;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//                               Serial Partitioning
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool zoltanGeometricPartitionerSerial( const uint size, const uint ncube_total, const uint offset, const int method,
                                       struct Zoltan_Struct *zz, const Center_coords &XYZ, real *weight, Zoltan_Out *zoltan_out,
                                       int comsize )
{
    float ver;
    int rc;

/* General parameters */

#if ( 1 )

    // char csize=static_cast<char>(5);

    /* this works only for single digits not safe
       char csize;

       csize='0'+comsize; // this will '0'+int ---> 'int'
    */
    char csize[16];

    sprintf( csize, "%d", comsize );

    //   std::cout<<"????????????????????????????????   "<<csize<<endl;

    //   printf("%d",comsize);
    // this is added for serial version, it is to
    Zoltan_Set_Param( zz, "NUM_LOCAL_PARTS", csize );

    Zoltan_Set_Param( zz, "IMBALANCE_TOL", TOL );

    Zoltan_Set_Param( zz, "DEBUG_LEVEL", "2" );
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

//    printf( "size=%d\n", size );

#if ( 1 )
    MeshData myMesh;
    //	myMesh=(MeshData*)malloc(1.0*sizeof(MESH_DATA));

    myMesh.numGlobalPoints = ncube_total;
    myMesh.numMyPoints = size;
    myMesh.myGlobalIDs = nullptr;
    //myMesh.myGlobalIDs = (ZOLTAN_ID_PTR)malloc( size * sizeof( ZOLTAN_ID_TYPE ) );
    myMesh.myGlobalIDs = new ZOLTAN_ID_TYPE [ size] ;

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
        //	printf(" %lf\n",myMesh.w[i]);
    }

    // printf("my_rank =%d number of my points %d glob_id global_id %d %d
    // \n",my_rank,myMesh->numMyPoints,myMesh->myGlobalIDs[0],myMesh->myGlobalIDs[1]);

    /* Query functions, to provide geometry to Zoltan */

    Zoltan_Set_Num_Obj_Fn( zz, get_number_of_objects, &myMesh );

    Zoltan_Set_Obj_List_Fn( zz, get_object_list, &myMesh );

    Zoltan_Set_Num_Geom_Fn( zz, get_num_geometry, &myMesh );

    Zoltan_Set_Geom_Multi_Fn( zz, get_geometry_list, &myMesh );

    // double t1,t2;
    // t1=MPI_Wtime();
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
#endif

#endif
    bool success = true;
    if ( rc != ZOLTAN_OK )
    {
        std::cout << RED "Warning: Zoltan Partition Rejected" RESET << endl;
        success = false;
        //        MPI_Finalize();
        //        exit( 0 );
    }

    // printf( "my_rank %d changes %d numExport %d numImport %d\n",my_rank,zoltan_out->changes,zoltan_out->numExport,zoltan_out->numImport);

    delete[] myMesh.c;
    delete[] myMesh.w;
    delete[] myMesh.myGlobalIDs;
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
    // printf("number of my points %d\n",mesh->numMyPoints);
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
