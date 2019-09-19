#ifndef _PHDF5_H_
#define _PHDF5_H_
#include "definitions.h"

/*!    \class Phdf5
 *     \brief  This Writes out Tree data in hdf5 format in parallel with *.xmf as metadata suitable for paraview and visit
 *
 */

template <size_t N,typename Nvalue,size_t M,typename Mvalue>
class Forest;


template <size_t N, typename Nvalue,size_t M,typename Mvalue>
class Phdf5
{
   private:
   uint totalnumber;


    public:
    Phdf5(){};	
    void writePolyvertex( Forest<N, Nvalue, M, Mvalue> &F, uint appx ); /*!< writes only the centeroids of the mesh, appx sets the appendix as string for the output file  */
    void xdmfPolyvertex( integer my_rank, uint appx );
    void writeMultiBlock( Forest<N, Nvalue, M, Mvalue> &F,uint appx ); /*!< writes each element as block and mesh is combination of blocks, appx sets the appendix as string for the output file  */
    
    void xdmfMultiBlock( Forest<N, Nvalue, M, Mvalue> &F, integer comsize, integer my_rank,uint offset, uint appx );

#if(0)
    void writeHdf5PolyvertexSerial( Tree<N, value> &T, uint appx );
    void writeP( Tree<N, value> &T, uint iterate ); /*!< write only the centroids of the cubes */

    void writeHdf5MultiBlockSerial( Tree<N, value> &T, uint appx ); /*!<writes out xdmf metadata required to read the hdf5 file */
    void xdmfMultiBlockSerial( Tree<N, value> &T, uint appx );      /*!<writes out the mesh in hdf5  */
    void write( Tree<N, value> &T, uint iterate );                  /*!< write all of the the elements */

    void writeHdf5MultiBlockHighLevel( Tree<N, value> &T, uint appx );
    void xdmfMultiBlockHighLevel( Tree<N, value> &T, uint appx );
    void writeH( Tree<N, value> &T, uint appx ); /*!< write only the elements with highest level*/
#endif
~Phdf5(){};
};


#if(0)
/*!    \class Hdf5XmfV
 *   
 *  *\brief  This class writes out Voxel data in hdf5 format with *.xmf as metadata suitable for paraview and visit,
 *  \details the reason for keeping these separate is because Tree data will have a solution vector soon and
 *    that class will be much more different down the road
 *   */

template <size_t N, typename value>
class Hdf5XmfV
{
    public:
    // Hdf5Xmf;

    void xdmfPolyvertexSerial( Voxel<N, value> &V );
    void writeHdf5PolyvertexSerial( Voxel<N, value> &V );
    void writeP( Voxel<N, value> &V );                    /*!< write only the centroids of the cubes */
    void writeHdf5MultiBlockSerial( Voxel<N, value> &V ); /*!<writes out all the elements */
    void xdmfMultiBlockSerial( Voxel<N, value> &V );
    void write( Voxel<N, value> &V ); /*!<writes out only the elements with highest level in hdf5  */
};

#endif

#endif
