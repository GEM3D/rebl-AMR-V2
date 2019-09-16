#ifndef _FOREST_H_
#define _FOREST_H_
#include "communicate.h"
#include "definitions.h"
#include "tree.h"
#include "typedefs.h"

/*!
 * \class Forest
 * \brief  template class that is a forest of octrees with semi-structured process topology
 * \details
 *  includes a list of tree's
 *  and functionality for manipulation as well as
 *  exchange of the trees with neighboring processes
 *  The algorithm starts with a very coarse 16 bit
 *  semi-structured processor topology, for now only 16 bits
 *  are used but it can be modified according to the need
 *  each process will have one forest and trees will be distributed
 *  dynamically as the solution porogresses
 */

template <size_t N, typename value>
using treelist = std::list<Tree<N, value>>;

template <size_t M> /*This map is defined to find the reference of the tree in the list that corresponds to the seed */
using seed = std::unordered_map<morton<M>, uint>;

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
class Forest
{
    protected:
    real ancestorlength[3]; /*!< original length of the first generation (root) element*/
    real ancestorcoords[3]; /*!< centeroid of the  of the first generation (root) element*/
    uint npx;               /*!< discritization in x direction, this value for proc tree is 2, therefore forest needs its own value of npx*/
    uint npy;               /*!< discretization in y direction */
    uint npz;               /*!< discretization in z direction */

    private:
    MpiCom Com;
    treelist<N, Nvalue> trees;                                         /*!< list of trees that each processor includes*/
    bitlist<M> seeds; /*!< seeds: morton code for boxes to grow tree*/ /*maybe in the future add this to class tree*/
    uint maxseedlevel; /*!< finds tha maximum level of the seeds */
    Tree<M, real> geom; 
    vector<uint> destination;
    vector <uint> nbrsOfNbrs;
    struct Zoltan_Struct *zz=nullptr; 
    Zoltan_Out zoltan_out;
    MPI_Comm graphComm;
    /*!< Each seed will contain its own geometry points to search, this assumption implicitly coincides processor topology with geometry
    voxelization
    this way tree's root will be morton  code used in geometry voxelization and no extra operation is necessary */
    // seed<M> seeds; /*!< look up the tree that need to be transfered to other procs */
    // may be use tuple in the future*/
    // std::tuple<morton<M>,Tree<N,Nvalue>> lookup();

    public:
    Forest<N, Nvalue, M, Mvalue>( real *length, real *coords, Tree<M, Mvalue> &proc,const int fixedlevel, uint nx, uint ny, uint nz ); /*!<constructor */
    void assignSeeds( real *length, Tree<M, Mvalue> &proc,const int ficedlevel);
    uint getTotalSize();
//    void formProcTop(Tree<M, Mvalue> &proc ,uint proclevel,real *geom_xyz,int geom_nn);
    void refineEachTreeVoxel( uint nlevel ); /*!<Refines every Tree in the list, nlevels, balance is satisfired for each tree but not in the
                                           boundaries with other trees */
    void assignGeom( Tree<M, Mvalue> &proc,const uint fixedlevel ,real *geom_xyz, uint geom_nn );
    void refineEachTree( uint nlevel ); /*!<Refines every Tree in the list, nlevels, balance is satisfired for each tree but not in the global scope*/
    void getListEachTree();
    void fourToOneBalance( Tree<M, uint> &proc ); /*!< 4:1 balance enforced at forest including the other processors */
    bool isInSeed( morton<M> &key, uint *counter ); /*!< check and see if the forest incldues the particluar seed*/
    void flipAll( morton<N> &key, uint *mylevel, uint *direction );
    // composite morton code functions
    void findFlipLevel( morton<N + M> key, uint *mylevel, uint *changedirectionlevel,
                        uint *direction ); /*!<same as the function defined in class three  except that it workd on (M+N) bits */
    void flipForNbr( morton<N + M> &key, uint *mylevel, uint *changedirectionlevel,
                     uint *direction ); /*!<same as the function defined in class three  except that it workd on (M+N) bits */
    // bool getNbrSeedLevel(morton<N+M> &combinedkey,uint topologylevel,uint *nbrseedleve,Tree<M,uint>& proc);
    // void flipForNbr( morton<M> &key, uint *mylevel, uint *changedirectionlevel, uint *direction );
//    void combine(const morton<N> &key,const morton<M>  &seed, morton<M+N> &combinedkey); /*!<combined seed key with element key*/
    void getNbrSeedLevel( morton<N + M> &combinedkey, uint topologylevel, uint *nbrseedleve, Tree<M, uint> &proc );
    void debug( Tree<M, Mvalue> &proc );
    // typename bitlist<M>::iterator find( morton<M> key );
    void getElemNbrs(Tree<M, Mvalue> &proc,const morton<M> key, bitvector<M> &nbr ); /*!< collects the neghbors of a given element */
    void comPatternConstruct( Tree<M, Mvalue> &proc );
    void getDirections( morton<N + M> &key, uint combinedlevel, vector<uint> &directions );
    void encodeGeometry();
    void removeAllZeroSingularity(morton<N+M> &key,const uint &combinedlevel);    
    void getMaxSeedsLevel(Tree<M,Mvalue>&  proc);
    void findSeedLevelForRcvdMessage(const morton<N+M> &key, uint *mylevel, Tree<M,Mvalue> & proc);
    void recoverAllZeroSingularity(morton<N+M> &key, const uint &combinedlevel); /*!< note that this function operates on the element key, this is done to remove redundant calc*/
    void constructSeedKeyForRcvdMessage(const morton<N+M> &key,const uint &seedlevel, morton<M> &seedkey);
    void constructElementKeyForRcvdMessage(const morton<N+M> &key,const uint &seedlevel, morton<N> &elementkey);
    // void refine();
    void refineForestBalanced( uint nlevel,Tree<M,Mvalue>& proc );
    void combinedLevel(const morton<N+M> &key, uint *level );
    void zoltanGeomrepart(  Tree<M,Mvalue>& proc,uint setmethod );
   // void EncodeGeomEachTree( Tree<M, Mvalue> &proc, real *geom_xyz, uint n );
    uint forestsize();
    void retainFourToOneBalance( Tree<M, uint> &proc );
    void moveGeom( Tree<M, Mvalue> &proc, const uint fixedlevel, real *geom_xyz,uint n, real x[3]); /*!< moves the geomerty with displacements specified in x[3] in x,y and z directions */
    void pushToDerefineEachTree( uint nlevel, Tree<M, uint> &proc ); 
    void convertBitsToDouble(  morton<N+M> key, double *val);
    void convertDoubleToBits(  morton<N+M> &key, const double val);
    void createCommGraph(uint Nnbr); /*!< level of neighbors */
    void createNbrsOfNbrs( );     
 // debug routines
    void debugDerefine( Tree<M, Mvalue> &proc );
    void checkGraphConsistency();
    void checkNbrsOfNbrsConsistency();
    bool checkWithNbrs(bool *sendbuf,bool *recvbuf); /*!< this is to enforce broadcast only with first degree neighbors */
    void rcvrMessageSize(int *sendbuf, int *recvbuf); 
    void getTotalMeshSize();
    void checkZoltanPartConsistency(Tree<M,Mvalue>&proc );
// friends and such
    void constructCommWeak(const vector<int>Nbr);

    template <size_t N1, typename Nvalue1, size_t M1, typename Mvalue1>
    friend class Phdf5;

/*
    template <size_t N2>
    friend class Scale;
*/

// destructor 

    ~Forest() {};
};

#endif
