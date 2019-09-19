#ifndef _INTERPOLATE_H_
#define _INTERPOLATE_H_
#include "definitions.h"
#include "params.h"



/* \class  interpolate
  * \brief  incorporates the interpolation techniques across the boundariesi, methods in this class are static and hence
  *         there is no need for object creation just calli t with :: and you should be good to go 
  *    
  */


typedef void (*func_ptr)(double ,double ,double *);



template<typename Nvalue>
class interpolate
{
/*
private:
func_ptr *f1=NULL;
*/
public:
interpolate(){};


static void quad(const real Q0,const real Q1,const real Q2,const real ksi , real *Qout );

static void restrictFace(const Nvalue *Q,const int npx,const int npy,const int npz, const int direction,const int faceID,Nvalue *Q_out);

static void line(const real Q0,const real Q1,const real ksi , real *Qout );

void setPointers();

static void quad2D( Nvalue *Qin, const real ksi, const real zeta, Nvalue *Qout );
//static void test()

static void quad1D(const double q0, const double q1, const double q2, double const delx,const double x , double *qbc );

static void getInterpolantNodes(const Nvalue *soln ,const int nrowWght,  const int istart, const int jstart, Nvalue *QresultBase);

static void fetchFaceChunks(const Nvalue *point0, const int pxg , const int pyg, const int pzg, const int direction , const int faceTag, Nvalue **pointChunk);

static void fetchFaceF2C(const Nvalue *point0, const int pxg , const int pyg, const int pzg, const int direction , const int faceTag, Nvalue *face, Nvalue* innerFace);

static void interpolateFace(const Nvalue *pointChunk , const int nColumnChunkWGst, const int nRowChunkWGst, Nvalue *Qinterpolated);

static void restrictFace( const Nvalue *solution , const int nColumnChunkWGst, const int nRowChunkWGst, Nvalue *Qrestricted);

static void interpolateGhostCells(Nvalue* point, const int pxg, const int pyg, const int pzg, const int direction);

static void quadraticOneD (Nvalue *point, const int indexGhost, const int index1, const int index2);

//this function does the transverse interpolation on the ghost values received (post-swap interpolation)
void updateRcvdGhstValWithnCoarseBlock(Nvalue *point0, const int pxg, const int pyg, const int pzg, const int direction,const int faceTag);

//this function does the transverse interpolation on the ghost values received (post-swap interpolation)
void updateRcvdGhstValWithnFineBlock(Nvalue *point0, const int pxg, const int pyg, const int pzg, const int direction,const int faceTag);


//this function does the transverse interpolation on the send buffer (pre-swap interpolation) when data exchanged from fine block to coarse one 
void updateSndBufferPreSwap(Nvalue *point0, const int pxg, const int pyg, const int pzg, const int direction,const int faceTag,Nvalue*pointF2C);


};


#endif

