#ifndef _REBLAMRFULL_H_
#define _REBLAMRFULL_H_
#include "communicate.h"
#include "definitions.h"
#include "geomSTL.h"
#include "templateForest.h"
#include "templatePhdf5.h"
#include "tree.h"
#include "typedefs.h"
/*!    \class RebleAmr
 * \brief  This is the main class for Rebl-AMR and it uses all other classes
 *developed in this project locally
 *
 *
 */

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
class ReblAmrFull
{
    private:
    real                                                  xyz1[6];
    FullTree<WSIZE,uint>                                  Proc; 
    GeomSTL                                               GMT;          /*!< object to read and manipulate the STL file  */
    int                                                   proclevel;    /*!< Initial refinement level for domain decomposition  */
    int                                                   meshlevel;    /*!< refinement level for each seed */
//    Partition                                             Part;         /*!< object used for partitioning  */
    MpiCom                                                Com;          /*!< communicator of the class */
    TemplateForest<N, Nvalue, M, Mvalue, FullTree<M,Mvalue>> Forest;                     /*!< main forest class  */
    vector<uint> Nbrs;


    public:
    ReblAmrFull<N, Nvalue, M, Mvalue>( int argcs, char *pArgs[], real *length, real *coords, uint nx, uint ny,
                                   uint nz ); /*!< class constructor*/
    void MPIStartUp(); /*!< duplicates and starts up its own Communicator */
#if(1)
    void countPointsinBox();                  /*!< Calculates the weights to be used for Partitioning Algorithms   */
   void generateProcTopology();              /*!<  Generates processor Topology which is a lower level Tree*/
    void setForestParams();          /*!<sets some parameters for forest */
   void createComPattern();         /*!< constructs communication pattern */
    void forestConstruct( int argcs, char *pArgs[], real *length, real *coords, uint nx, uint ny,
                          uint nz ); /*!< constructor of the forest object */
    void moveGeometry( double *xx ); /*<! moves the geometry given the xx coordinates of its center of mass  */
    void getTotalMeshSize();         /*!< return the global mesh size */
    void refineForest();             /*!< performs refinement on forest */
// this is the functions for regulat tree
//if(std::is_same<P, Tree<PROCSIZE,uint>>::value)
// these functions do not make sense for class FullTree and hence their invocation will return an error
 //   void assignPartToProc();                  /*!< Partitions and assigns each element to the corresponding process accordignly  */
//    void asignWeightsForPart();               /*!< assigns wights for each element for weighted partitioning */
//    void distributeTopology();

   void writeMesh(int index);                /*!< writes out mesh in hdf5 format */
   void writeRunInfo();             /*! writes out the log file*/
#endif
    ~ReblAmrFull(); /*!< destructor of the class*/
};

#endif
