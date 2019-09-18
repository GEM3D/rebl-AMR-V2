#ifndef _TEMPLATEPHDF5_H_
#define _TEMPLATEPHDF5_H_
#include "definitions.h"

/*!    \class TemplatePhdf5
 *     \brief  This Writes out Tree data in hdf5 format in parallel with *.xmf as metadata suitable for paraview and visit
 *
 */

template <size_t N,typename Nvalue,size_t M,typename Mvalue, class T>
class templateForest;


template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
class templatePhdf5
{
   private:
   uint totalnumber;


    public:
    templatePhdf5(){};	
   // void writePolyvertex( templateForest<N, Nvalue, M, Mvalue, T> &F, uint appx ); /*!< writes only the centeroids of the mesh, appx sets the appendix as string for the output file  */
    void writePolyvertex(TemplateForest<N,Nvalue,M,Mvalue,T> &F,  uint appx ); /*!< writes only the centeroids of the mesh, appx sets the appendix as string for the output file  */
//    void writePolyvertex( T &F, uint appx ); /*!< writes only the centeroids of the mesh, appx sets the appendix as string for the output file  */

    void xdmfPolyvertex( integer my_rank, uint appx );

    void writeMultiBlock(TemplateForest<N,Nvalue,M,Mvalue,T> &F, uint appx);

// for the weak scaling where we hve billions of elements, multiblock might overflow due to the growth of the text file 
// beyond buffer limits  
//
   void xdmfMultiBlock(TemplateForest<N, Nvalue, M, Mvalue,T> &F, integer comsize ,integer my_rank, uint offset,uint appx );


#if(0)
 

    void writeMultiBlock( Forest<N, Nvalue, M, Mvalue> &F,uint appx ); /*!< writes each element as block and mesh is combination of blocks, appx sets the appendix as string for the output file  */
    
    void xdmfMultiBlock( Forest<N, Nvalue, M, Mvalue> &F, integer comsize, integer my_rank,uint offset, uint appx );
   void writeHdf5PolyvertexSerial( Tree<N, value> &T, uint appx );
    void writeP( Tree<N, value> &T, uint iterate ); /*!< write only the centroids of the cubes */

    void writeHdf5MultiBlockSerial( Tree<N, value> &T, uint appx ); /*!<writes out xdmf metadata required to read the hdf5 file */
    void xdmfMultiBlockSerial( Tree<N, value> &T, uint appx );      /*!<writes out the mesh in hdf5  */
    void write( Tree<N, value> &T, uint iterate );                  /*!< write all of the the elements */

    void writeHdf5MultiBlockHighLevel( Tree<N, value> &T, uint appx );
    void xdmfMultiBlockHighLevel( Tree<N, value> &T, uint appx );
    void writeH( Tree<N, value> &T, uint appx ); /*!< write only the elements with highest level*/
#endif
~templatePhdf5(){};
};

#endif
