#ifndef _TEMPLATEFOREST_H_
#define _TEMPLATEFOREST_H_
#include "communicate.h"
#include "definitions.h"
#include "tree.h"
#include "typedefs.h"
#include "solver.h"
#include "interpolate.h"

/*!
 * \class TemplateForest
 * \brief  Template Class designed to unify tree and full_tree topologies 
 */

template <class T>
using treeList = std::list<T>; /*!< list container to hold the trees  */

template <size_t M> /*!< This map is defined to find the reference of the tree in the list that corresponds to the seed */
using seed = std::unordered_map<morton<M>, uint>;


//template <size_t N, typename Nvalue, size_t M, typename Mvalue, template< size_t L, typename Lvalue> class T=Tree >
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T >
class TemplateForest
{
    protected:
    real ancestorlength[3]; /*!< original length of the first generation (root) element*/
    real ancestorcoords[3]; /*!< centeroid of the  of the first generation (root) element*/
    uint npx;               /*!< discritization in x direction, this value for proc tree is 2, therefore forest needs its own value of npx*/
    uint npy;               /*!< discretization in y direction */
    uint npz;               /*!< discretization in z direction */
   
    uint px; /*!< number of cell centered elements in x-direction */
    uint py; /*!< number of cell centered elements in y-direction */
    uint pz; /*!< number of cell centered elements in z-direction */
    
    uint pxg; /*!< number of cell centered elements in x-direction plus ghost cells */
    uint pyg; /*!< number of cell centered elements in y-direction plus ghost cells */
    uint pzg; /*!< number of cell centered elements in z-direction plus ghost cells */
   

    MPI_Datatype sendType[9];
    MPI_Datatype contiguous; /*!< Contiguous Data type for the {f,p} pairs */
    MPI_Datatype tmpStride;


    private:
    MpiCom Com; /*!< Object to hold communicator data such as rank, size, etc  */
    int duplicated=0;
    treeList<Tree<N,Nvalue>> trees;                                         /*!< list of trees that each processor includes*/
    uint maxProcLevel; /*! < Maximum level in processor topology. i.e. tree*/
    bitlist<M> seeds; /*!< seeds: morton code for boxes to grow tree*/ /*maybe in the future add this to class tree*/
    uint maxseedlevel; /*!< finds tha maximum level of the seeds */
    Tree<M, real> geom;  
    vector<uint> destination; /*!< neighbor processes to communicate messages  */
    vector<uint> sendtag; 
    vector<uint> recvtag;
    vector <uint> nbrsOfNbrs;
    vector<bitset<M + N>> *message=nullptr;
    struct Zoltan_Struct *zz=nullptr; /*<Zoltan structure to be used for partitioning */ 
    Zoltan_Out zoltan_out;
    MPI_Comm graphComm;
    MPI_Request *request=nullptr;
    MPI_Request *request1=nullptr;
    char   nameAppendix[80];
    unsigned long long meshSize;
    solver<Nvalue> S;
    interpolate<Nvalue>Intrp;
 
    ofstream fdeb;
    /*!< Each seed will contain its own geometry points to search, this assumption implicitly coincides processor topology with geometry
    voxelization
    this way tree's root will be morton  code used in geometry voxelization and no extra operation is necessary */
    // seed<M> seeds; /*!< look up the tree that need to be transfered to other procs */
    // may be use tuple in the future*/
    // std::tuple<morton<M>,Tree<N,Nvalue>> lookup();

    public:
//    TemplateForest<N, Nvalue, M, Mvalue, T>( T &proc, real *length, real *coords, uint nx, uint ny, uint nz ); /*!<constructor */
          
    TemplateForest<N, Nvalue, M, Mvalue, T>( ){}; /*!<constructor */
    TemplateForest<N, Nvalue, M, Mvalue, T>(int argcs,char *pArgs[], T &proc, real *length, real *coords, uint nx, uint ny, uint nz ); /*!<constructor */
    void construct(int argcs,char *pArgs[], T &proc, real *length, real *coords, uint nx, uint ny, uint nz );
    uint getTotalSize();
    void assignSeeds( real *length,T &proc);
//    void formProcTop(Tree<M, Mvalue> &proc ,uint proclevel,real *geom_xyz,int geom_nn);
    void assignGeom( T &proc,real *geom_xyz, uint geom_nn );
    void encodeGeometry(); /*!< encode the geometry once to the deepest level  */
    void refineEachTree( uint nlevel ); /*!<Refines every Tree in the list, nlevels, balance is satisfired for each tree but not in the global scope*/
    void moveGeom( T &proc, real *geom_xyz,uint n, real x[3]); /*!< moves the geomerty with displacements specified in x[3] in x,y and z directions */
    void getListEachTree();
    bool isInSeed( morton<M> &key, uint *counter ); /*!< check and see if the forest incldues the particluar seed*/
    void flipAll( morton<N> &key, uint *mylevel, uint *direction );
    void getDirections( morton<N + M> &key, uint combinedlevel, vector<uint> &directions );
    void recoverAllZeroSingularity(morton<N+M> &key, const uint &combinedlevel); /*!< note that this function operates on the element key, this is done to remove redundant calc*/
    void combinedLevel(const morton<N+M> &key, uint *level );
    void findSeedLevelForRcvdMessage(const morton<N+M> &key, uint *mylevel, Tree<M,Mvalue> & proc);
    void findSeedLevelForRcvdMessage(const morton<N+M> &key, uint *mylevel, FullTree<M,Mvalue> & proc);
    void constructSeedKeyForRcvdMessage(const morton<N+M> &key,const uint &seedlevel, morton<M> &seedkey);
    void constructElementKeyForRcvdMessage(const morton<N+M> &key,const uint &seedlevel, morton<N> &elementkey);
    void removeAllZeroSingularity(morton<N+M> &key,const uint &combinedlevel);   /*!<we use one bit to communicate the level of the element ending in zero at element level  */  
    void getMaxSeedsLevel(T &proc);
    void findFlipLevel( morton<N + M> key, uint *mylevel, uint *changedirectionlevel,
                        uint *direction ); /*!<same as the function defined in class three  except that it workd on (M+N) bits */
    void flipForNbr( morton<N + M> &key, uint *mylevel, uint *changedirectionlevel,
                     uint *direction ); /*!<same as the function defined in class three  except that it workd on (M+N) bits */
    
    void getTotalMeshSize(); /*!< calculates the final number of the mesh generated with all processes */
    void getNbrSeedLevel( morton<N + M> &combinedkey, uint topologylevel, uint *nbrseedleve, Tree<M,Mvalue> &proc );
    void getNbrSeedLevel( morton<N + M> &combinedkey, uint topologylevel, uint *nbrseedleve, FullTree<M,Mvalue> &proc );
    uint forestsize();
    void getElemNbrs(Tree<M, Mvalue> &proc,const morton<M> key, bitvector<M> &nbr ); /*!< collects the neghbors of a given element */
    void comPatternConstruct( Tree<M, Mvalue> &proc ); // might be a bug in here
    void comPatternConstruct( FullTree<M, Mvalue> &proc, vector<uint> &Nbrs); 
    void fourToOneBalance( T &proc ); /*!< 4:1 balance enforced at forest including the other processors */
    void refineForestBalanced( uint nlevel,T &proc );
    void nonCollectiveNbrComm(  );

    // composite morton code functions
    // void flipForNbr( morton<M> &key, uint *mylevel, uint *changedirectionlevel, uint *direction );
    // typename bitlist<M>::iterator find( morton<M> key );

    void pushToDerefineEachTree( uint nlevel, Tree<M, uint> &proc ); 
    void retainFourToOneBalance( Tree<M, uint> &proc ); /*!< Eliminates the elements that their removal would destroy the 2:1 balance*/
    void checkZoltanPartConsistency(Tree<M,Mvalue>&proc ); /*! Zoltan might give zero elements to a processor,
                                                           this function checks to see each element has at least one 
                                                           element assigned to it and that the number of elements before and after partitioning is the same */

    void createCommGraph(uint Nnbr); /*!< level of neighbors */ // works for both of the objects if we use the proper getElemNbrs
    void createNbrsOfNbrs( );  // this one uses neighbor connectivity   
    void checkGraphConsistency(); /*! Sanity check for the graph, it checks he fact that if element "A" is a neighbor of element "B", element "B" should be a neighbor of element "A" as well  */
    void filterRefineList();
    void setMaxProcLevel(const uint refinelevel );   
  
    //uint findIndexInSeed( Tree<M,uint> &proc,  morton<M> &seedkey );
    uint findIndexInSeed( T &proc,  morton<M> &seedkey );

//  void appendToMessage(Tree<M, uint> &proc,const uint seednbrlevel ,  morton<M+N> &seednbrkey ,morton<M+N> &combinedkey,const uint combinedlevel);

    void appendToMessage( T &proc, morton<M> &seednbrkey ,morton<M+N> &combinedkey,const uint combinedlevel );


#if(0)
// derefine routines
   void zoltanGeomrepart(  Tree<M,Mvalue>& proc,uint setmethod );
// debug routnes
    void debug( Tree<M, Mvalue> &proc );
    void debugDerefine( Tree<M, Mvalue> &proc );
    void checkNbrsOfNbrsConsistency();
    bool checkWithNbrs(bool *sendbuf,bool *recvbuf); /*!< this is to enforce broadcast only with first degree neighbors */
    void rcvrMessageSize(int *sendbuf, int *recvbuf); 
    void constructCommWeak(const vector<int>Nbr);
// graph communication trial
#endif
    void MPIStartUp(); /*!< Creates an MPI communicator for ReblAmr and assigns ranks */
    void checkInputParams(int argcs,char *pArgs[]); /*!< Checks input parameters for consistency */
    void runInfo();
    void currentDateTime();
    template <size_t N1, typename Nvalue1, size_t M1, typename Mvalue1,class T1 >
    friend class templatePhdf5;

// utilities added for poisson solver
   void allocatePointersforEachBlock(T &proc);
   void initializeTrigonometricForest(T &proc);
   void swapGhostCells(T &proc);


   void addToRefineLisGradBased();
   void commit();
   void swapGhostCellsXdirectionSiblings(T &proc );
   void swapGhostCellsYdirectionSiblings(T &proc );
   void swapGhostCellsZdirectionSiblings(T &proc );
   void swapGhostCellsXdirectionCousins(T &proc );
   void swapGhostCellsYdirectionCousins(T &proc );
   void swapGhostCellsZdirectionCousins(T &proc );
   void flipForSibling(morton<N+M> &globalKey,uint combinedLevel, uint direction );

void extractSameGlobalLevelSiblingInfo( T &proc, const morton<M> seedkey, const morton<N > key, const uint topologylevel, const uint mylevel, const uint direction, morton<N+M> &globalKey, morton<M> &seednbrkey, morton<N> &nbrkey, uint &nbrseedlevel );

int extractSameGlobalLevelCousinInfo( T &proc, const morton<M> seedkey, const morton<N > key, const uint topologylevel, uint mylevel, uint direction,uint changedirectionlevel  ,morton<N+M> &globalKey, morton<M> &seednbrkey, morton<N> &nbrkey, uint &nbrseedlevel );


  bool isBoundary(T &proc, const morton<M> seedkey,const morton<N>key,const uint direction,const uint level);
//void extractRealNbrInfo(  morton<M> seednbrkey,  morton<N > nbrkey,uint estimatedNbrLevel, uint &realLevel );
int extractRealNbrInfo(  morton<M> seednbrkey, morton<N> nbrkey,const int estimatedNbrLevel, uint &realLevel );

  void getPointers( morton<M> seedKey, vector<morton<N>> &nbrs , vector<Q*> &ptrs );

void getPointersForNegtiveEstimate( morton<M> seednbrkey,uint nbrseedlevel,uint direction ,vector<morton<M>> &nbrs ,vector<Q *> &ptrs );
void getPointersForNegtiveEstimateCousin( morton<M> seednbrkey,uint nbrseedlevel,uint direction ,vector<morton<M>> &nbrs ,vector<Q *> &ptrs );
 void constructHigherLevelNbrs( morton<N+M> globalKey, uint combinedlevel, uint direction,vector<morton<N+M>>&nbrs);


  void swapSameLevel(int direction, Q *point0, Q* point1,morton<N+M> globalkey ,  morton<M> seednbrkey, vector<morton<N>> &nbrs, uint combinedlevel,int
&count, MPI_Request *sendReq,  MPI_Request *rcvReq);

  void swapSameLevelCousins(int direction, Q *point0, Q* point1,morton<N+M> globalkey ,  morton<M> seednbrkey, vector<morton<N>> &nbrs, uint combinedlevel,int
&count, MPI_Request *sendReq,  MPI_Request *rcvReq);


  void swapWithHigher(int direction, Q *point0, Q* point1,morton<N+M> globalkey ,  morton<M> seednbrkey, vector<morton<N>> &nbrs, uint combinedlevel,int
&count, MPI_Request *sendReq,  MPI_Request *rcvReq);

void swapWithHigherForNegtiveEstimate( int direction, Q *point0, vector<Q *> &ptrs, morton<N + M> globalKey,
                                                              morton<M> seednbrkey, vector<morton<M>> &nbrs, uint combinedlevel, int &count,
                                                              MPI_Request *sendReq, MPI_Request *rcvReq );


void swapWithHigher(int direction, Q* point0, vector<Q*> &ptrs,morton<N+M> globalKey , morton<M> seednbrkey, int *index, uint combinedlevel,int &count, MPI_Request *sendReq,  MPI_Request *rcvReq);

void swapWithLowerLevelCousin( int direction, Q *point0, vector<Q *> &ptrs, morton<N + M> globalKey,
                                                              morton<M> seednbrkey, vector<morton<N>> &nbrs, uint combinedlevel, int &count,
                                                              MPI_Request *sendReq, MPI_Request *rcvReq );

void swapWithHigherLevelCousin( int direction, Q *point0, vector<Q *> &ptrs, morton<N + M> globalKey,
                                                              morton<M> seednbrkey, int *index, uint combinedlevel, int &count,
                                                              MPI_Request *sendReq, MPI_Request *rcvReq );

void constructTrueNbrKeysSiblings( morton<M>seednbrkey,uint topologylevel, morton<N>nbrkey,uint mylevel, uint direction, vector<morton<N>> &nbr, int flag);
void constructTrueNbrKeysCousins( morton<M>seednbrkey,uint topologylevel, morton<N>nbrkey,uint mylevel, uint direction, vector<morton<N>> &nbr, int flag);
void solvePoisson(T &proc,const real &Epsilon,const char *bc, const real *Dirichlet, const real *Neumann);

void imposeBoundaryConditionXdirection( T &proc );
void imposeBoundaryConditionYdirection( T &proc );
void imposeBoundaryConditionZdirection( T &proc );
void imposeBoundaryConditions( T &proc );

void tagBorderElementXdirectionSiblings( T &proc );
void initializeTrigonometricForMMS( T &proc );
void orderAnalysis( T &proc, const real &Epsilon, const char*bc, const real* Dirichlet, const real* Neumann );
void calcError(T &proc);

void quadraticInterpTransverseDirection( T &proc );

void quadInterpTransDirXdirectionSiblings(T &proc );
void quadInterpTransDirYdirectionSiblings(T &proc );
void quadInterpTransDirZdirectionSiblings(T &proc );

void quadInterpTransDirXdirectionCousins(T &proc );
void quadInterpTransDirYdirectionCousins(T &proc );
void quadInterpTransDirZdirectionCousins(T &proc );

void preSwaptestTransInterpolation( T &proc );
void testTransInterpolationPrint2Display(Nvalue *point, int direction, int row, int column );
void postSwaptestTransInterpolation(T &proc);
void postSwapPrint2Display(Nvalue *point);

~TemplateForest(); /*!< Destructor of the object*/

   // ~TemplateForest(){}; /*!< Destructor of the object*/

};


#endif
