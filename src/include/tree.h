#ifndef _TREE_H_
#define _TREE_H_
#include "definitions.h"

/*!    \class Tree
 * \brief  This Class Generates a 4:1 balancerd AMR mesh
 *
 *
 *
 */

template <size_t N, typename value>
class Hdf5Xmf;

template <size_t N, typename value>
class Tree
{
    template <size_t N1, typename value2>
    friend class Hdf5Xmf; /*!< this is a friend class to write out in hdf5 format*/

    protected:
    real ancestorlength[3];    /*!< original length of the first generation (root) element*/
    real ancestorcoords[3];    /*!< centeroid of the  of the first generation (root) element*/
    morton<N> ancestorkey = 0; /*!< root value is always set as 0000000000000   */
    uint npx;                  /*!< discritization in x direction*/
    uint npy;                  /*!< discretization in y direction */
    uint npz;                  /*!< discretization in y direction */

    private:
    bitmap<N, value> mesh;                         /*!< base main container */
    std::unordered_map<morton<N>, int> refinelist; /*!<  list of elements tagged to be refined  */
    bitvector<N> mortonSTL;                        /*!< vector to store geometric points in morton code */
    // bitvector<N> derefinelist; /*!< list of elements tagged to be removed*/
    //                               might have some of its elements removed, due to this reason, I do not use vector but rather list*/
    //    bitvector<N> mortonSTL ;  /*!< vector to store geometric points in morton code */
    std::unordered_map<morton<N>, int> derefinelist; /*!< list of elements tagged to be removed, due to 4:1 balance this listi*/

    public:
    bitvector<N> boundarylist; /*!< list of elements of refinelist that are boundary elements*/
    // I wanted seed to each tree but that requires another template parameter M for seed which is normally less than
    // N , i.e. M < N
    //  might not be bad to do this, I have not decided yet

    public:
    Tree( real *length, real *coords, uint nx, uint ny, uint nz ); /*!< constructor*/
    Tree( real *length, real *coords );                            /*!< constructor*/
    Tree() {};
    void construct( real *length, real *coords, uint nx, uint ny, uint nz ); /*!< need to initialize inside forest*/
    virtual void level( morton<N> key, uint *level );                        /*!< obtains the level of the element from morton key */
    void centroid( morton<N> key, real *xyz );   /*!< calculates the centroid of the cube given the morton key of the element*/
    void enclosingBox( morton<N> key, real *X ); /*!< calculates the range that an element occupies in 3D space for a given Element*/
    virtual typename bitmap<N, value>::iterator begin(); /*!< iterator returning the first object*/
    virtual typename bitmap<N, value>::iterator end();   /*!< iterator returning the last object*/
    virtual uint size();                                 /*!< returns the size of the mesh*/
    void reserve( uint *reservedsize );                  /*!<this function reserves the memory given the reservedsize of the mesh */
    void siblings( morton<N> key, uint mylevel, morton<N> *sibkey ); /*!< extracts the siblings from morton code */
    void refine( morton<N> key );                                    /*!< perfomrs refinement for a tagged element given the Morton Key*/
    void derefine( morton<N> key );                                  /*! performs derefinement on a single element given a morton key*/
    void refineRefineList();                                         /*!<performs the refinement */
    void refineRefineList(bitvector<N>&V);                                         /*!<performs the refinement */
    void fourToOneP( uint istart,
                     uint iend ); /*!< imposes 4:1 balance locally for each tree, while loop is eliminated due to parallel implementation*/

    void refineRefineList( uint istart, uint iend ); /*!<performs the refinement */
    void fourToOne(); /*!< imposes 4:1 balance given the list of elments to be refined in the vector refine list*/
    void findFlipLevel( morton<N> key, uint *mylevel, uint *changedirectionlevel,
                        uint *direction ); /*!< detects the flip level, this info used in finding nonlocal neighbors*/
    void flipForNbr( morton<N> *key, uint *mylevel, uint *changedirectionlevel,
                     uint *direction );       /*!< perform the actual operation to identify the nonlocal nbr */
    uint IsInVectorList( morton<N> key );     /*!< checks to see if a given code is already in the list */
    void addToList( morton<N> key );          /*!<adds element to refinelist */
    uint count( morton<N> key );              /*!< counts the number of elements*/
    void addToDerefineList( morton<N> key );  /*!< adds the element to derefinelist */
    void derefineDerefineList( uint nlevel ); /*!< Derefines the mesh*/
                                              // geometry
    uint isInsideSolid( const morton<N> key, const real *geom_xyz,
                        uint n ); /*!< tags the elements if any points of the gemoetry resides in the enclosing box*/
    virtual typename bitmap<N, value>::iterator find( morton<N> key ); /*!< this function is to find a value given the key*/
    void insertKey( morton<N> key );
    void convertStl2Morton( uint geom_size, real *geom_xyz ); /*!<converts stl coordinates to morton and puts them in morton STL*/
                                                              //    template<size_t N1>
    void pushToRefinelist( uint level ); /*!<Note that for dynamic mesh we need to make sure the element exists before adding to this list
                                            */
    bool isBoundary( morton<N> &key );
    void extractBoundary();
    void getDirections( morton<N> &key, vector<uint> &directions );
    //    void getElemNbr(morton<N> &key);
    uint refineListSize();
    void clearRefineList();
    void extractBoundaryP( uint istart, uint iend );
    bool isInMeshList( const morton<N> &key );
    bool isInRefineList( const morton<N> &key );
    void constructHigherLevelNbrs( const morton<N> &key, const uint &keylevel, const uint &direction, morton<N> *nbr );
    void printMesh();
    bool isBoundary( const morton<N> &key, uint direction );
    //  void readRefineList(morton<N>&, uint);
    std::pair<morton<N>, int> readRefineList( typename std::unordered_map<morton<N>, int>::iterator it );
    morton<N> readDerefineList( typename std::unordered_map<morton<N>, int>::iterator it );
    // typename std::unordered_map<morton<N>,int>::iterator  readDerefineList(morton<N> key);
    void getKey( uint i, morton<N> &key );
    void clearMortonSTL();
    // derefine functions
    void retainFourToOne();
    void removeFromDerefineList( typename std::unordered_map<morton<N>, int>::iterator it );
    typename unordered_map<morton<N>, int>::iterator Dbegin();
    typename unordered_map<morton<N>, int>::iterator Dend();
    void derefineDerefineList();
    void clearMesh();
    void pushToDerefinelist( uint nlevel );
    typename unordered_map<morton<N>, int>::iterator Rbegin();
    typename unordered_map<morton<N>, int>::iterator Rend();
    typename std::unordered_map<morton<N>, int>::iterator findInDerefine( morton<N> key ); /*!<  Get the iterator from derefinelist*/
    void mortonSTLclear();
    void flipRefineElemTag( typename std::unordered_map<morton<N>, int>::iterator it );
    void refinelistReset();
    void constructNonlocalHigherLevelNbrs( const morton<N> &key, const uint &keylevel, const uint &direction, morton<N> *nbr );
    void insertSeed( morton<N> &key );

    void insertNbrs( vector<int> &Nbrs );
    void enclosingBoxFixedLevel( morton<N> key, uint mylevel, real *X );
    void centroidFixedLevel( morton<N> key, const uint mylevel, real *xyz );

    virtual void convertCoordToMorton( real *xyz, morton<N> &key ); /*!<converts coordinates of a point to morton code */
    typename std::unordered_map<morton<N>, int>::iterator  findInList( morton<N> key );                                                        

    void setToZero(  );
    void ignoreInactive(morton<N> key );
    void ignoreInactive(int nInactive, real * box);
	void ignoreInactiveVertices(int nInactive, real * box);

     void  getCoords( real *X );

     bool operator==(const  Tree<N, value>  & T0) const;  

//    void derefineDerefineList();
/*

       template <size_t N1, typename Nvalue1, size_t N2>
       friend void convertCoordToMorton(Tree<N1,Nvalue1> T, real *xyz, morton<N2> &key );
    */

    /*
      template <size_t N1, typename Nvalue1, size_t M1, typename Mvalue1>
      friend class Phdf5;
  */

    ~Tree(); /*!< Destructor of the class, it frees the memeory pointed by pointer in the hashmap value  if allocated*/
};
// the reason I am keeping these separate is because of the need for the solution vector in the base class
// template<size_t N>
// class Hdf5XmfV;

/*!    \class Voxel
 * \brief  This Class Generates an unbalancerd Voxel to improve search by geometry partitioning
 *
*/

template <size_t N, typename value>
class Voxel : public Tree<N, value>
{
    private:
    uint maxlevel; /*!< maximum level of refinement*/
    uint numMax;   /*!< number of elements having the highest level*/
    bitmap<N, value> mesh;
    bitvector<N> lookup;

    public:
    Voxel( real *length, real *coords ) : Tree<N, value>( length, coords, 2, 2, 2 ) {}; /*!< constructor */
    void setLevel( uint *l );                                                           /*!< sets the maximum level for refinement*/
    void generateSearchTree( real *geom_xyz, uint n );                                  /*!< generates an initial tree*/
    void distributeGeomToLeaves( real *geom_xyz, uint n ); /*!< distributes geometry to different cells (leaves) */
    uint checkSiblingStatus( morton<N> key, morton<N> *sibkey );
    ;                        /*!<checks to see if siblings include any points and whether they have the same level */
    void derefineGeomTree(); /*!< derefines the tree based on geometry*/
    bool IsInsideSegment( morton<N> key, real *xyz );

    template <size_t N1, typename value1>
    friend class Hdf5XmfV; /*!< this is a friend class to write out in hdf5 format*/
    ~Voxel();              /*!< Destructor of this class*/
                           // add ete <size_t N, typename value>
                           // void Tree<N, value>::refine( morton<N> key )
    // rase because regular refine does not free the memory, though I reallocate in the next step, normally it should not leak
    // void  erase(morton key); /* rather specialized version of erase, since the memory has to be freed as well*/
};

/*!    \class FullTree
 * \brief  This Class is specifically designed for weak analysis to operate on the fulltree topology
 * in this class, level is preset by the user and the level function is redefined to make use of the
 * polymorphism
 *
*/

#if ( 1 )
template <size_t N, typename value>
class FullTree : public Tree<N, value>
{
    private:
    uint fixedlevel; /*!< maximum level of full tree*/
    bitmap<N, value> mesh;
    std::unordered_map<morton<N>, int> refinelist;   /*!<  list of elements tagged to be refined  */
    bitvector<N> mortonSTL;                          /*!< vector to store geometric points in morton code */
    std::unordered_map<morton<N>, int> derefinelist; /*!< list of elements tagged to be removed, due to 4:1 balance this listi*/

    public:
    FullTree() {};
    FullTree( real *length, real *coords ) : Tree<N, value>( length, coords, 2, 2, 2 ) {}; /*!< constructor */
    void setLevel( const uint &l );                                                        /*!< sets the maximum level for refinement*/
    void level( morton<N> key, uint *level );
    void insertKey( morton<N> key );
    uint getLevel();
    // this class has its own mesh as private so it needs ti have its own iterators as well, inheritance does not do this automatically
    void nbrsConstrcut( vector<uint> &Nbrs, uint myrank );
    bool isBoundary( uint &direction, uint myrank );
    typename bitmap<N, value>::iterator find( morton<N> key ); /*!< this function is to find a value given the key*/
    uint size();
    //
    void convertCoordToMorton(real* xyz , morton<N>&key);
    void assignProcs( vector<uint> &Nbrs, uint myrank );
    void findFlipLevel( morton<WSIZE> key,uint fixedlevel, uint *changedirectionlevel, uint *direction ) ;   
    void flipForNbr( morton<WSIZE> *key, uint fixedlevel, uint *changedirectionlevel, uint *direction );

    typename bitmap<N, value>::iterator begin(); /*!< iterator returning the first object*/
    typename bitmap<N, value>::iterator end();
//     void centroid( morton<N> key, real *xyz );   /*!< calculates the centroid of the cube given the morton key of the element*/
//     void enclosingBox( morton<N> key, real *X ); /*!< calculates the range that an element occupies in 3D space for a given Element*/
    ~FullTree() {}; /*!< Destructor of this class*/
    // add erase because regular refine does not free the memory, though I reallocate in the next step, normally it should not leak
    // void  erase(morton key); /* rather specialized version of erase, since the memory has to be freed as well*/
};
#endif

#endif
