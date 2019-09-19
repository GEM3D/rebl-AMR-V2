#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "definitions.h"
#include "params.h"

template <typename Nvalue>
class solver
{
    private:
    uint nx, ny, nz;
    uint nxg, nyg, nzg;
    char bc[6];
	real Dirichlet[6];
	real Neumann[6];
  	int  tags[3]; /*!< extract the info for exact function construction*/

  
    public:
    solver(){};
    void setGrid( uint n0, uint n1, uint n2 );
    // initialize exact

    void initializeTrigonometric( Nvalue* q, real* XYZ );
    real getGradient( Nvalue* q, real* XYZ );
    real updateInterior( Nvalue* q, real* XYZ );

    void initializeTrigonometricForMMS( Nvalue* q, real* XYZ );
    void setBoundaryConditions( Nvalue* q );
    void exactValueForGhost( Nvalue* point, real* XYZ );

    void setBC(const char* bc , const real* Dirichlet, const real* Neumann );

    // notes to shams
    // please add the following functions
    // so we have  6 faces each face can be for now D or N,
    // like we did with direction and facetag
    // depending on the char 'D' of 'N' write a function to perform  the following
    // as usual we will keep all directions separately and the have a function to call all 3
    // similar to what we did in swapGhostCells()
    void imposeBoundaryXdirection( Nvalue* point, int faceTag , real* XYZ);
    void imposeBoundaryYdirection( Nvalue* point, int faceTag, real* XYZ );
    void imposeBoundaryZdirection( Nvalue* point, int faceTag, real* XYZ );
    // void imposeBoundary();
    // a method that takes in the frequency and location
   
    real exactValue( double omega, double x, int tag ); /*!< this is to automatically switch between sin and cos */
    void extractTag(); /*!< extracts tags (int) from bc (character)*/

    real getError( Nvalue *q, real *XYZ ); /*!< calculate l1-norm */
    real redBlackGS(Nvalue *q, real *XYZ); /*!< red black Gauss Sidel */
	// add the folowing
    // GS
    // RB GS
    // CG
    // GMRES ...
    //
};

#endif
