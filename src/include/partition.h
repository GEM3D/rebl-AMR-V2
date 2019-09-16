#include <stdlib.h>
#include "typedefs.h"
#include "tree.h"
#include "zoltan.h"
#include "communicate.h"
#include "params.h"
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




class Partition{

private:
int size;
struct Zoltan_Struct *zz=NULL;
//Zoltan_Out *zoltan_out=NULL;
MpiCom Com; 
MeshData myMesh;
int method;
Center_coords XYZ;
real *weight=nullptr;

public:
Partition(int argcs,char *pArgs[],int meth,int sze);
Partition(){};
void construct(int argcs,char *pArgs[],int meth,int sze  );
void MPIStartUp();
void partition(double *ID, int totalvalue);

bool zoltanGeometricPartitioner( const uint size, const uint ncube_total, const uint offset);

void zoltanSetParams();

//bool zoltanGeometricPartitionerSerial( const uint size, const uint ncube_total, const uint offset, int comsizei, struct Zoltan_Struct *zz, Zoltan_Out *zoltan_out );
bool zoltanGeometricPartitionerSerial( const uint size, const uint ncube_total, const uint offset, int comsize, Zoltan_Out *zoltan_out );

void reSize(uint sze);
/*
friend int get_number_of_objects( void *data, int *ierr );

friend void get_object_list( void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim,
                             float *obj_wgts, int *ierr );
friend int  get_num_geometry( void *data, int *ierr );

friend void get_geometry_list( void *data, int sizeGID, int sizeLID, int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                               int num_dim, double *geom_vec, int *ierr );
*/
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
friend class ReblAmr;

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
friend class ReblAmrFull;


~Partition();

};


















