#include "interpolate.h"
#include "params.h"

// function pointers for B--Quadratic interpolation in the AMR mesh
#if ( 0 )
void basis( double ksi, double zeta_Var, double *phi );
void phi0( double ksi, double zeta_Var, double *phi );
void phi1( double ksi, double zeta_Var, double *phi );
void phi2( double ksi, double zeta_Var, double *phi );
void phi3( double ksi, double zeta_Var, double *phi );
void phi4( double ksi, double zeta_Var, double *phi );
void phi5( double ksi, double zeta_Var, double *phi );
void phi6( double ksi, double zeta_Var, double *phi );
void phi7( double ksi, double zeta_Var, double *phi );
void phi8( double ksi, double zeta_Var, double *phi );

void basis( double ksi, double zeta_Var, double *phi )
{
    phi[0] = ksi * zeta_Var * ( ksi - 1.0 ) * ( zeta_Var - 1.0 ) * ( 1.0 / 4.0 );
    phi[1] = ksi * zeta_Var * ( ksi + 1.0 ) * ( zeta_Var - 1.0 ) * ( 1.0 / 4.0 );
    phi[2] = ksi * zeta_Var * ( ksi + 1.0 ) * ( zeta_Var + 1.0 ) * ( 1.0 / 4.0 );
    phi[3] = ksi * zeta_Var * ( ksi - 1.0 ) * ( zeta_Var + 1.0 ) * ( 1.0 / 4.0 );
    phi[4] = -zeta_Var * ( ksi * ( 1.0 / 2.0 ) - 1.0 / 2.0 ) * ( ksi + 1.0 ) * ( zeta_Var - 1.0 );
    phi[5] = ksi * ( ksi + 1.0 ) * ( zeta_Var - 1.0 ) * ( zeta_Var + 1.0 ) * ( -1.0 / 2.0 );
    phi[6] = -zeta_Var * ( ksi * ( 1.0 / 2.0 ) - 1.0 / 2.0 ) * ( ksi + 1.0 ) * ( zeta_Var + 1.0 );
    phi[7] = ksi * ( ksi - 1.0 ) * ( zeta_Var - 1.0 ) * ( zeta_Var + 1.0 ) * ( -1.0 / 2.0 );
    phi[8] = ( ksi - 1.0 ) * ( ksi + 1.0 ) * ( zeta_Var - 1.0 ) * ( zeta_Var + 1.0 );
}

void phi0( double ksi, double zeta_Var, double *phi ) { ( *phi ) = ksi * zeta_Var * ( ksi - 1.0 ) * ( zeta_Var - 1.0 ) * ( 1.0 / 4.0 ); }

void phi1( double ksi, double zeta_Var, double *phi ) { *phi = ksi * zeta_Var * ( ksi + 1.0 ) * ( zeta_Var - 1.0 ) * ( 1.0 / 4.0 ); }

void phi2( double ksi, double zeta_Var, double *phi ) { *phi = ksi * zeta_Var * ( ksi + 1.0 ) * ( zeta_Var + 1.0 ) * ( 1.0 / 4.0 ); }

void phi3( double ksi, double zeta_Var, double *phi ) { *phi = ksi * zeta_Var * ( ksi - 1.0 ) * ( zeta_Var + 1.0 ) * ( 1.0 / 4.0 ); }

void phi4( double ksi, double zeta_Var, double *phi )
{
    *phi = -zeta_Var * ( ksi * ( 1.0 / 2.0 ) - 1.0 / 2.0 ) * ( ksi + 1.0 ) * ( zeta_Var - 1.0 );
}

void phi5( double ksi, double zeta_Var, double *phi )
{
    *phi = ksi * ( ksi + 1.0 ) * ( zeta_Var - 1.0 ) * ( zeta_Var + 1.0 ) * ( -1.0 / 2.0 );
}

void phi6( double ksi, double zeta_Var, double *phi )
{
    *phi = -zeta_Var * ( ksi * ( 1.0 / 2.0 ) - 1.0 / 2.0 ) * ( ksi + 1.0 ) * ( zeta_Var + 1.0 );
}

void phi7( double ksi, double zeta_Var, double *phi )
{
    *phi = ksi * ( ksi - 1.0 ) * ( zeta_Var - 1.0 ) * ( zeta_Var + 1.0 ) * ( -1.0 / 2.0 );
}

void phi8( double ksi, double zeta_Var, double *phi ) { *phi = ( ksi - 1.0 ) * ( ksi + 1.0 ) * ( zeta_Var - 1.0 ) * ( zeta_Var + 1.0 ); }

// func_ptr *f1=NULL;
// f1=(func_ptr*)malloc(9*sizeof(func_ptr));

template <typename Nvalue>
void interpolate<Nvalue>::setPointers()
{
    f1[0] = phi0;
    f1[1] = phi1;
    f1[2] = phi2;
    f1[3] = phi3;
    f1[4] = phi4;
    f1[5] = phi5;
    f1[6] = phi6;
    f1[7] = phi7;
    f1[8] = phi8;
}

#endif

template <typename Nvalue>
void interpolate<Nvalue>::quad( const real Q0, const real Q1, const real Q2, const real ksi, real *Qinterpolated )
{
    // Using Langrange Polynomials
    // these are used in FEM regularly

    double c0 = 0.5 * ( ksi - 1.0 ) * ksi;
    double c1 = ( ksi + 1.0 ) * ( 1.0 - ksi );
    double c2 = 0.5 * ( 1.0 + ksi ) * ksi;

    *Qinterpolated = c0 * Q0 + c1 * Q1 + c2 * Q2;
}

template <typename Nvalue>
void interpolate<Nvalue>::line( const real Q0, const real Q1, const real ksi, real *Qinterpolated )
{
    // Using Langrange Polynomials

    double c0 = ( 1.0 - ksi ) * 0.5;
    double c1 = ( ksi + 1.0 ) * 0.5;

    *Qinterpolated = c0 * Q0 + c1 * Q1;
}

template <typename Nvalue>
void interpolate<Nvalue>::quad1D( const double q0, const double q1, const double q2, double const delx, const double x, double *qbc )
{
    double rhs0 = q1 - q0;
    double rhs1 = q2 - q0;

    double c01 = delx;
    double c00 = c01 * c01;

    double c11 = delx * 2.5;
    double c10 = c11 * c11;

    double det = c00 * c11 - c01 * c10;

    if ( fabs( det ) < 1.e-15 )
    {
        cout << "det =0.0 " << endl;
        exit( 0 );
    }

    double a0 = 1. / det * ( c11 * rhs0 - c01 * rhs1 );
    double a1 = 1. / det * ( -c10 * rhs0 + c00 * rhs1 );

    // printf("det = %lf a0=%lf a1=%lf \n",det,a0,a1);

    *qbc = q0 + a0 * ( x * x ) + a1 * x;
}

template <typename Nvalue>
void interpolate<Nvalue>::quad2D( Nvalue *Qin, const real ksi, const real zeta, Nvalue *Qinterpolated )
{
// Using Langrange Polynomials
// these are used in FEM regularly
#if ( 1 )
    double c0 = 0.5 * ( ksi - 1.0 ) * ksi;
    double c1 = ( ksi + 1.0 ) * ( 1.0 - ksi );
    double c2 = 0.5 * ( 1.0 + ksi ) * ksi;
    double d0 = 0.5 * ( zeta - 1.0 ) * zeta;
    double d1 = ( zeta + 1.0 ) * ( 1.0 - zeta );
    double d2 = 0.5 * ( 1.0 + zeta ) * zeta;

    double c00 = c0 * d0;
    double c10 = c1 * d0;
    double c20 = c2 * d0;
    double c01 = c0 * d1;
    double c11 = c1 * d1;
    double c21 = c2 * d1;
    double c02 = c0 * d2;
    double c12 = c1 * d2;
    double c22 = c2 * d2;
    /*
    Qinterpolated[0].p=c00*Qin[0].p+c10*Qin[1].p+c20*Qin[2].p
         +c01*Qin[3].p+c11*Qin[4].p+c21*Qin[5].p
         +c02*Qin[6].p+c12*Qin[7].p+c22*Qin[8].p;
    */

    Qinterpolated[0]
    = Qin[0] * c00 + Qin[1] * c10 + Qin[2] * c20 + Qin[3] * c01 + Qin[4] * c11 + Qin[5] * c21 + Qin[6] * c02 + Qin[7] * c12 + Qin[8] * c22;

#endif
    /*
    double v=0.0,val1;

         for(int l2=0;l2<9;l2++)
             {
               v=v+Qin[l2]*val1;
    //	   f1[l2](ksi,zeta,&val1);
             }
    *Qinterpolated=v;
    */
}

template <typename Nvalue>
void interpolate<Nvalue>::restrictFace( const Nvalue *solution, const int nColumnChunkWGst, const int nRowChunkWGst, Nvalue *Qrestricted )
{
    int     jstart = 1;
    int     istart = 1;
    Nvalue *Qchunk;
    Nvalue *Qsingle;

    Qchunk  = new Nvalue[9];
    Qsingle = new Nvalue[1];

    real ksi  = -0.5;
    real zeta = -0.5;

    int nRowRestrictedNoGst    = ( nRowChunkWGst - 2 ) / 2;
    int nColumnRestrictedNoGst = ( nColumnChunkWGst - 2 ) / 2;

    for ( int jndex = 0; jndex < nRowRestrictedNoGst; jndex++ )
    {
        for ( int index = 0; index < nColumnRestrictedNoGst; index++ )
        {
            // if chunk is at the bottom boundary
            if ( jstart + 3 > nRowChunkWGst - 1 )
            {
                jstart = nRowChunkWGst - 4;
                zeta   = 0.5;
            }

            // if the chunk is at the right boundary
            if ( ( istart + 3 ) >= nColumnChunkWGst - 1 )
            {
                istart = nColumnChunkWGst - 4;
                ksi    = 0.5;
            }
            getInterpolantNodes( solution, nColumnChunkWGst, istart, jstart, Qchunk ); // extract 9 points
            // print the chunk
            //     cout<<"\n index = "<<index<<", "<<" istart = : "<<istart<<endl;

            quad2D( Qchunk, ksi, zeta, Qsingle );
            Qrestricted[nColumnRestrictedNoGst * jndex + index] = Qsingle[0];
            istart += 2;
        }

        jstart += 2;
        ksi    = -0.5;
        istart = 1;
    }
    delete[] Qchunk;
    delete[] Qsingle;
}

template <typename Nvalue>
void interpolate<Nvalue>::interpolateFace( const Nvalue *pointChunk, const int nColumnChunkWGst, const int nRowChunkWGst,
                                           Nvalue *Qinterpolated )
{
    Nvalue *QresultBase;
    Nvalue *QresultFromBase;
    QresultBase     = new Nvalue[9];
    QresultFromBase = new Nvalue[1];
    real ksi        = -0.75;
    real zeta       = -0.75;

    int nColumnDoubleG = 2 * ( nColumnChunkWGst - 2 ) + 2; // with ghos:
    int nRowDoubleG    = 2 * ( nRowChunkWGst - 2 ) + 2;    // with ghost


    for ( int jndex = 0; jndex < nRowDoubleG; jndex++ )
    {
        for ( int index = 0; index < nColumnDoubleG; index++ )
        {
            int istart = 2 * ( index / 4 );
            int jstart = 2 * ( jndex / 4 );

            if ( istart + 3 >= nColumnChunkWGst )
            {
                istart = nColumnChunkWGst - 3;
                if ( index == nColumnDoubleG - 2 )
                {
                    ksi = 0.25;
                }
            }
            if ( jstart + 3 >= nRowChunkWGst )
            {
                jstart = nRowChunkWGst - 3;
                if ( jndex == nRowDoubleG - 2 )
                {
                    zeta = 0.25;
                }
            }
            if ( ksi > 0.75 )
            {
                ksi = -0.75;
            }

            if ( zeta > 0.75 )
            {
                zeta = -0.75;
            }

            getInterpolantNodes( pointChunk, nColumnChunkWGst, istart, jstart, QresultBase ); // extract 9 points

            // interpolate the new value at ksi and zeta
            quad2D( QresultBase, ksi, zeta, QresultFromBase );
            Qinterpolated[nColumnDoubleG * jndex + index] = QresultFromBase[0];
            // cout<<" interpolated is :"<<QresultFromBase[0].p <<endl;
            ksi += 0.5;
        }
        // next row begins
        ksi = -0.75;
        zeta += 0.5;
    }
    delete[] QresultBase;
    delete[] QresultFromBase;
}

template <typename Nvalue>
void interpolate<Nvalue>::getInterpolantNodes( const Nvalue *soln, const int nrowWght, const int istart, const int jstart,
                                               Nvalue *QresultBase )
{
    int jend = jstart + 3;

    int iend = istart + 3;

    for ( int j = jstart; j < jend; j++ )
    {
        for ( int i = istart; i < iend; i++ )
        {
            QresultBase[3 * ( j - jstart ) + ( i - istart )] = soln[(nrowWght)*j + i];
        }
    }
}



template <typename Nvalue>
void interpolate<Nvalue>::updateRcvdGhstValWithnFineBlock(Nvalue *point0, const int pxg, const int pyg, const int pzg, const int direction,const int faceTag)
{
	/*!< This function interpolates the Ghost value in fine blocks after receiving the Ghost value from the coarse blocks 
 	 *   Using lagrange interpolation
 	 */

real  Xg = 0.5; /*!<coordinate where the ghost needs to be calculated */
real  X1 = 0.0; /*!<coordinate of the first point(ghost)received from the coarse block*/
real  X2 = 1.5; /*!<coordinate of the second point inside the fine block*/
real  X3 = 2.5; /*!<coordinate of the third point inside the fine block*/


real  L1  = (Xg - X2)*(Xg - X3)/((X1 - X2)*( X1 - X3));
real  L2  = (Xg - X1)*(Xg - X3)/((X2 - X1)*( X2 - X3));
real  L3  = (Xg - X1)*(Xg - X2)/((X3 - X1)*( X3 - X2));

// q@xu = 

// x direction 
    if ( direction == 0 )
    {
        switch ( faceTag )
        {
            case 0: 

                for ( int k = 0; k < pzg; k++ )
                {
                    for ( int j = 0; j < pyg; j++ )
                    {
                        int indexGhost = ( pxg ) * (pyg)*k + (pxg)*j  +  0;
						int point1    = indexGhost;
                        int point2    = indexGhost + 1;
                        int point3    = indexGhost + 2;
               
                        point0[indexGhost] =  point0[point1]*L1 + point0[point2]*L2 + point0[point3]*L3;
                
                         
                    }
                }

                break;

            case 1: // right face 

                for ( int k = 0; k < pzg; k++ )
                {
                    for ( int j = 0; j < pyg; j++ )
                    {
                    
				            int indexGhost = ( pxg ) * (pyg)*k + (pxg)*j + pxg-1;
                			int point1    = indexGhost;
                        	int point2    = indexGhost - 1;
                        	int point3    = indexGhost - 2;
               
                        point0[indexGhost] =  point0[point1]*L1 + point0[point2]*L2 + point0[point3]*L3;
                                
                                
                


					}
                }

                break;
        }
    }
    else if ( direction == 1 ) // y direction 
    {
        switch ( faceTag )
        {
            case 0: // left face j = 1;

                for ( int k = 0; k < pzg; k++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {

							int indexGhost = ( pxg ) * (pyg)*k + (pxg)*0 + i;
                 			int point1    = indexGhost;
	           				int point2 = ( pxg ) * (pyg)*k + (pxg)*(1)+ i;
	           				int point3 = ( pxg ) * (pyg)*k + (pxg)*(2)+ i;
 
              
                        point0[indexGhost] =  point0[point1]*L1 + point0[point2]*L2 + point0[point3]*L3;
               




                    }
                }

                break;

            case 1: // right face j = pyg-2

                for ( int k = 0; k < pzg; k++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {
	           				int indexGhost = ( pxg ) * (pyg)*k + (pxg)*(pyg-1)+ i;
               				int point1 = indexGhost;
	           				int point2 = ( pxg ) * (pyg)*k + (pxg)*(pyg-2)+ i;
	           				int point3 = ( pxg ) * (pyg)*k + (pxg)*(pyg-3)+ i;
                             
                        point0[indexGhost] =  point0[point1]*L1 + point0[point2]*L2 + point0[point3]*L3;
               



                   }
                }

                break;
        }
    }
    else // Z direction 
    {
        switch ( faceTag )
        {
            case 0: // left face k = 1;

                for ( int j = 0; j < pyg; j++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {

                        int indexGhost  = ( pxg ) * (pyg)*0 + (pxg)*j + i;
                        int point1      = indexGhost ;
                        int point2		= ( pxg ) * (pyg)*1 + (pxg)*j + i;
                        int point3 		= ( pxg ) * (pyg)*2 + (pxg)*j + i;
                        point0[indexGhost] =  point0[point1]*L1 + point0[point2]*L2 + point0[point3]*L3;


                    }
                }

                break;

            case 1: // right face k = pzg -2;

                for ( int j = 0; j < pyg; j++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {

                        int indexGhost = ( pxg ) * (pyg)*(pzg-1) + (pxg)*j + i;
						int point1     = indexGhost;
                        int point2		= ( pxg ) * (pyg)*(pzg-2) + (pxg)*j + i;
                        int point3 		= ( pxg ) * (pyg)*(pzg-3) + (pxg)*j + i;
                        point0[indexGhost] =  point0[point1]*L1 + point0[point2]*L2 + point0[point3]*L3;


                
                    }
                }

                break;
        }
    }
}


template <typename Nvalue>
void interpolate<Nvalue>::updateRcvdGhstValWithnCoarseBlock(Nvalue *point0, const int pxg, const int pyg, const int pzg, const int direction,const int faceTag)
{
	/*!< This function interpolates the Ghost value in coarse blocks after receiving the Ghost value from the fine blocks 
 	 *   Using lagrange interpolation
 	 */


real  Xg = 1.0; /*!<coordinate where the ghost needs to be calculated */
real  X1 = 0.5; /*!<coordinate of the first point in fine block*/
real  X2 = 1.5; /*!<coordinate of the second point infine block*/
real  X3 = 3.0; /*!<coordinate of the third point inside the coarse block*/

real  L1  = (Xg - X2)*(Xg - X3)/((X1 - X2)*( X1 - X3));
real  L2  = (Xg - X1)*(Xg - X3)/((X2 - X1)*( X2 - X3));
real  L3  = (Xg - X1)*(Xg - X2)/((X3 - X1)*( X3 - X2));



// x direction 
    if ( direction == 0 )
    {
        switch ( faceTag )
        {
            case 0: // left face i = 1

                for ( int k = 0; k < pzg; k++ )
                {
                    for ( int j = 0; j < pyg; j++ )
                    {
                        int i = 0;

                        int indexGhost = ( pxg ) * (pyg)*k + (pxg)*j + i;
                                      
                        int point3     = ( pxg ) * (pyg)*k + (pxg)*j + (i+1);
                                      
                        point0[indexGhost] =  point0[indexGhost] + point0[point3]*L3;
                
                         
                    }
                }

                break;

            case 1: // right face i = pxg-2

                for ( int k = 0; k < pzg; k++ )
                {
                    for ( int j = 0; j < pyg; j++ )
                    {
                        int i = pxg - 1;
                    
				           int indexGhost = ( pxg ) * (pyg)*k + (pxg)*j + i;
                                
                            int point3  = ( pxg ) * (pyg)*k + (pxg)*j + (i-1);
                                
			  			    point0[indexGhost] =  point0[indexGhost] + point0[point3]*L3;
                


					}
                }

                break;
        }
    }
    else if ( direction == 1 ) // y direction 
    {
        switch ( faceTag )
        {
            case 0: // left face j = 1;

                for ( int k = 0; k < pzg; k++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {
                        int j = 0;

						int indexGhost = ( pxg ) * (pyg)*k + (pxg)*j + i;
                                      
                        int point3  = ( pxg ) * (pyg)*k + (pxg)*(j+1) +i;
                                      
                  		point0[indexGhost] =  point0[indexGhost] + point0[point3]*L3;
                




                    }
                }

                break;

            case 1: // right face j = pyg-2

                for ( int k = 0; k < pzg; k++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {
 
					  int j = pyg - 1;
 
	                     int indexGhost = ( pxg ) * (pyg)*k + (pxg)*j + i;

	                     int point3 = ( pxg ) * (pyg)*k + (pxg)*(j-1) + i;
                
  						 point0[indexGhost] =  point0[indexGhost] + point0[point3]*L3;
                

                   }
                }

                break;
        }
    }
    else // Z direction 
    {
        switch ( faceTag )
        {
            case 0: // left face k = 1;

                for ( int j = 0; j < pyg; j++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {
                        int k = 0;

                        int indexGhost = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        int point3 = ( pxg ) * (pyg)*(k+1) + (pxg)*j + i;

						  point0[indexGhost] =  point0[indexGhost] + point0[point3]*L3;
                


                    }
                }

                break;

            case 1: // right face k = pzg -2;

                for ( int j = 0; j < pyg; j++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {
                        int k = pzg - 1;

                        int indexGhost = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        int point3 = ( pxg ) * (pyg)*(k-1) + (pxg)*j + i;
                
                  		  point0[indexGhost] =  point0[indexGhost] + point0[point3]*L3;
                
                    }
                }

                break;
        }
    }
}

template <typename Nvalue>
void interpolate<Nvalue>::fetchFaceF2C( const Nvalue *point0, const int pxg, const int pyg, const int pzg, const int direction,
                                        const int faceTag, Nvalue *face, Nvalue* innerFace )
{
    if ( direction == 0 )
    {    
        switch ( faceTag )
        {    
            case 0: // left face i = 1

                for ( int k = 0; k < pzg; k++ )
                {    
                    for ( int j = 0; j < pyg; j++ )
                    {    
                        int i = 1; 

                        int index      = ( pxg ) * (pyg)*k + (pxg)*j + i; 
                     
                        int innerIndex = ( pxg ) * (pyg)*k + (pxg)*j + (i+1); 
                    
                        // surface values    
                        face[k * pyg + j] = point0[index];
                     
                        // values below the surface
                        innerFace[k * pyg + j] = point0[innerIndex];
                              
                    }    
                }    

                break;

            case 1: // right face i = pxg-2

                for ( int k = 0; k < pzg; k++ )
                {    
                    for ( int j = 0; j < pyg; j++ )
                    {    
                        int i = pxg - 2; 
                         
                        int index      = ( pxg ) * (pyg)*k + (pxg)*j + i; 
                     
                        int innerIndex = ( pxg ) * (pyg)*k + (pxg)*j + (i-1); 
                    
                        // surface values    
                        face[k * pyg + j] = point0[index];
                    
                         // values below the surface
                        innerFace[k * pyg + j] = point0[innerIndex];

                    }
                }

                break;
        }
    }
 else if ( direction == 1 )
    {
        switch ( faceTag )
        {
            case 0: // left face j = 1;

                for ( int k = 0; k < pzg; k++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {
                        int j = 1;

                        // index for surface values
                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        // index for values below the surface
                        
                        int innerIndex = ( pxg ) * (pyg)*k + (pxg)*(j+1) + i;

                        //surface values
                        face[k * pxg + i]     = point0[index];

                        //below surface values
                        innerFace[k * pxg + i] = point0[innerIndex];

                    }
                }

                break;

            case 1: // right face j = pyg-2

                for ( int k = 0; k < pzg; k++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {
 
                      int j = pyg - 2;
 
                         int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                         int innerIndex = ( pxg ) * (pyg)*k + (pxg)*(j-1) + i;

                         face[k * pxg + i] = point0[index];

                         innerFace[k * pxg + i] = point0[innerIndex];
                    }
                }

                break;
        }
    }
                                                                       else
    {
        switch ( faceTag )
        {
            case 0: // left face k = 1;

                for ( int j = 0; j < pyg; j++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {
                        int k = 1;

                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        int innerIndex = ( pxg ) * (pyg)*(k+1) + (pxg)*j + i;

                        face[j * pxg + i]     = point0[index];

                        innerFace[j * pxg + i] = point0[innerIndex];


                    }
                }

                break;

            case 1: // right face k = pzg -2;

                for ( int j = 0; j < pyg; j++ )
                {
                    for ( int i = 0; i < pxg; i++ )
                    {
                        int k = pzg - 2;

                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        int innerIndex = ( pxg ) * (pyg)*(k-1) + (pxg)*j + i;

                        face[j * pxg + i]     = point0[index];

                        innerFace[j * pxg + i] = point0[innerIndex];

                    }
                }

                break;
        }
    }
}


template <typename Nvalue>
void interpolate<Nvalue>::fetchFaceChunks( const Nvalue *point0, const int pxg, const int pyg, const int pzg, const int direction,
                                           const int faceTag, Nvalue **pointChunk )
{
    //  Nvalue *pointChunk[4];
    int npxchunk = ( pxg - 2 ) / 2 + 2;
    int npychunk = ( pyg - 2 ) / 2 + 2;

    if ( direction == 0 )
    {
        switch ( faceTag )
        {
            case 0: // x direction , left face , where i = 1

                // chunk 0

                for ( int k = 0; k < pzg / 2 + 1; k++ )
                {
                    for ( int j = 0; j < pyg / 2 + 1; j++ )
                    {
                        int i = 1;
#if ( !POISSON )
                        i = 0;
#endif

                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[0][k * npychunk + j] = point0[index];
                    }
                }

                // chunk 1

                for ( int k = pzg / 2 - 1; k < pzg; k++ )
                {
                    for ( int j = 0; j < pyg / 2 + 1; j++ )
                    {
                        int i = 1;
#if ( !POISSON )
                        i = 0;
#endif

                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[1][( k - ( pzg / 2 - 1 ) ) * npychunk + j] = point0[index];
                    }
                }

                // chunk2

                for ( int k = 0; k < pzg / 2 + 1; k++ )
                {
                    for ( int j = pyg / 2 - 1; j < pyg; j++ )
                    {
                        int i = 1;
#if ( !POISSON )
                        i = 0;
#endif

                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[2][k * npychunk + j - ( pyg / 2 - 1 )] = point0[index];
                    }
                }

                // chunk 3

                for ( int k = pzg / 2 - 1; k < pzg; k++ )
                {
                    for ( int j = pyg / 2 - 1; j < pyg; j++ )
                    {
                        int i = 1;
#if ( !POISSON )
                        i = 0;
#endif

                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[3][( k - ( pzg / 2 - 1 ) ) * npychunk + j - ( pyg / 2 - 1 )] = point0[index];
                    }
                }

                break;

            case 1: // x direction , right face , where i = pxg-2

                // chunk 0

                for ( int k = 0; k < pzg / 2 + 1; k++ )
                {
                    for ( int j = 0; j < pyg / 2 + 1; j++ )
                    {
                        int i = pxg - 2;

#if ( !POISSON )
                        i = pxg - 1;
#endif

                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[0][k * npychunk + j] = point0[index];
                    }
                }

                // chunk 1

                for ( int k = pzg / 2 - 1; k < pzg; k++ )
                {
                    for ( int j = 0; j < pyg / 2 + 1; j++ )
                    {
                        int i = pxg - 2;

#if ( !POISSON )
                        i = pxg - 1;
#endif

                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[1][( k - ( pzg / 2 - 1 ) ) * npychunk + j] = point0[index];
                    }
                }

                // chunk2

                for ( int k = 0; k < pzg / 2 + 1; k++ )
                {
                    for ( int j = pyg / 2 - 1; j < pyg; j++ )
                    {
                        int i = pxg - 2;

#if ( !POISSON )
                        i = pxg - 1;
#endif

                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[2][k * npychunk + j - ( pyg / 2 - 1 )] = point0[index];
                    }
                }

                // chunk 3

                for ( int k = pzg / 2 - 1; k < pzg; k++ )
                {
                    for ( int j = pyg / 2 - 1; j < pyg; j++ )
                    {
                        int i = pxg - 2;

#if ( !POISSON )
                        i = pxg - 1;
#endif

                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[3][( k - ( pzg / 2 - 1 ) ) * npychunk + j - ( pyg / 2 - 1 )] = point0[index];
                    }
                }

                break;
        } // end of switch statement
    }
    else if ( direction == 1 ) // faces along y axis************************************************************************
    {
        switch ( faceTag )
        {
            case 0: // y direction , left face , where j = 1

                // chunk0

                for ( int k = 0; k < pzg / 2 + 1; k++ )
                {
                    int j = 1;

#if ( !POISSON )
                    j = 0;
#endif

                    for ( int i = 0; i < pxg / 2 + 1; i++ )
                    {
                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[0][k * npxchunk + i] = point0[index];
                    }
                }

                // chunk1
                for ( int k = pzg / 2 - 1; k < pzg; k++ )
                {
                    int j = 1;

#if ( !POISSON )
                    j = 0;
#endif

                    for ( int i = 0; i < pxg / 2 + 1; i++ )
                    {
                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[1][( k - ( pzg / 2 - 1 ) ) * npxchunk + i] = point0[index];
                    }
                }

                // chunk 2

                for ( int k = 0; k < pzg / 2 + 1; k++ )
                {
                    int j = 1;

#if ( !POISSON )
                    j = 0;
#endif

                    for ( int i = pxg / 2 - 1; i < pxg; i++ )
                    {
                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[2][k * npxchunk + i - ( pxg / 2 - 1 )] = point0[index];
                    }
                }

                // chunk 3

                for ( int k = pzg / 2 - 1; k < pzg; k++ )
                {
                    int j = 1;

#if ( !POISSON )
                    j = 0;
#endif

                    for ( int i = pxg / 2 - 1; i < pxg; i++ )
                    {
                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[3][( k - ( pzg / 2 - 1 ) ) * npxchunk + i - ( pxg / 2 - 1 )] = point0[index];
                    }
                }

                break;

            case 1: // y direction , right face , where j = pyg-2

                // chunk0

                for ( int k = 0; k < pzg / 2 + 1; k++ )
                {
                    int j = pyg - 2;

#if ( !POISSON )
                    j = pyg - 1;
#endif

                    for ( int i = 0; i < pxg / 2 + 1; i++ )
                    {
                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[0][k * npxchunk + i] = point0[index];
                    }
                }

                // chunk1
                for ( int k = pzg / 2 - 1; k < pzg; k++ )
                {
                    int j = pyg - 2;

#if ( !POISSON )
                    j = pyg - 1;
#endif

                    for ( int i = 0; i < pxg / 2 + 1; i++ )
                    {
                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[1][( k - ( pzg / 2 - 1 ) ) * npxchunk + i] = point0[index];
                    }
                }

                // chunk 2

                for ( int k = 0; k < pzg / 2 + 1; k++ )
                {
                    int j = pyg - 2;

#if ( !POISSON )
                    j = pyg - 1;
#endif

                    for ( int i = pxg / 2 - 1; i < pxg; i++ )
                    {
                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[2][k * npxchunk + i - ( pxg / 2 - 1 )] = point0[index];
                    }
                }

                // chunk 3

                for ( int k = pzg / 2 - 1; k < pzg; k++ )
                {
                    int j = pyg - 2;

#if ( !POISSON )
                    j = pyg - 1;
#endif

                    for ( int i = pxg / 2 - 1; i < pxg; i++ )
                    {
                        int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                        pointChunk[3][( k - ( pzg / 2 - 1 ) ) * npychunk + i - ( pxg / 2 - 1 )] = point0[index];
                    }
                }

                break;
        } // end of switch statement
    }
    else
    {
        // faces a long z axis ********************************************************************************************

        switch ( faceTag )
        {
            case 0: // z direction , left face , where k = 1;

                // chunk 0
                {
                    int k = 1;

#if ( !POISSON )
                    k = 0;
#endif

                    for ( int j = 0; j < pyg / 2 + 1; j++ )
                    {
                        for ( int i = 0; i < pxg / 2 + 1; i++ )
                        {
                            int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                            pointChunk[0][j * npxchunk + i] = point0[index];
                        }
                    }
                }
                // chunk 1
                {
                    int k = 1;

#if ( !POISSON )
                    k = 0;
#endif

                    for ( int j = pyg / 2 - 1; j < pyg; j++ )
                    {
                        for ( int i = 0; i < pxg / 2 + 1; i++ )
                        {
                            int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                            pointChunk[1][( j - ( pyg / 2 - 1 ) ) * npxchunk + i] = point0[index];
                        }
                    }
                }

                // chunk2

                {
                    int k = 1;

#if ( !POISSON )
                    k = 0;
#endif

                    for ( int j = 0; j < pyg / 2 + 1; j++ )
                    {
                        for ( int i = pxg / 2 - 1; i < pxg; i++ )
                        {
                            int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                            pointChunk[2][j * npxchunk + i - ( pxg / 2 - 1 )] = point0[index];
                        }
                    }
                }

                // chunk 3

                {
                    int k = 1;

#if ( !POISSON )
                    k = 0;
#endif

                    for ( int j = pyg / 2 - 1; j < pyg; j++ )
                    {
                        for ( int i = pxg / 2 - 1; i < pxg; i++ )
                        {
                            int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                            pointChunk[3][( j - ( pyg / 2 - 1 ) ) * npxchunk + i - ( pxg / 2 - 1 )] = point0[index];
                        }
                    }
                }
                break;

            case 1: // z direction , right face , where k = pzg-2
                // chunk 0
                {
                    int k = pzg - 2;

#if ( !POISSON )
                    k = pzg - 1;
#endif

                    for ( int j = 0; j < pyg / 2 + 1; j++ )
                    {
                        for ( int i = 0; i < pxg / 2 + 1; i++ )
                        {
                            int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                            pointChunk[0][j * npxchunk + i] = point0[index];
                        }
                    }
                }

                // chunk1
                {
                    int k = pzg - 2;

#if ( !POISSON )
                    k = pzg - 1;
#endif

                    for ( int j = pyg / 2 - 1; j < pyg; j++ )
                    {
                        for ( int i = 0; i < pxg / 2 + 1; i++ )
                        {
                            int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                            pointChunk[1][( j - ( pyg / 2 - 1 ) ) * npxchunk + i] = point0[index];
                        }
                    }
                }

                // chunk 2
                {
                    int k = pzg - 2;

#if ( !POISSON )
                    k = pzg - 1;
#endif

                    for ( int j = 0; j < pyg / 2 + 1; j++ )
                    {
                        for ( int i = pxg / 2 - 1; i < pxg; i++ )
                        {
                            int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                            pointChunk[2][j * npxchunk + i - ( pxg / 2 - 1 )] = point0[index];
                        }
                    }
                }

                // chunk 3
                {
                    int k = pzg - 2;

#if ( !POISSON )
                    k = pzg - 1;
#endif

                    for ( int j = pyg / 2 - 1; j < pyg; j++ )
                    {
                        for ( int i = pxg / 2 - 1; i < pxg; i++ )
                        {
                            int index = ( pxg ) * (pyg)*k + (pxg)*j + i;

                            pointChunk[3][( j - ( pyg / 2 - 1 ) ) * npxchunk + i - ( pxg / 2 - 1 )] = point0[index];
                        }
                    }
                }
                break;
        }
    }
}

/* this function uses 3 point to construct a quadratic polynomial*/
template <typename Nvalue>
void interpolate<Nvalue>::quadraticOneD( Nvalue *point, const int indexGhost, const int index1, const int index2 )
{
    /* f(x) = a0 + a1*x + a2*x*x */
    /* f1(0 = QGhost, f2 = Q1, f3 = Q2 */
    /* we use lagrange's formula here
     * p2(x) = y0*L0(x) + y1*L1(x) + y2*L2(x);
     * X0 = 1;
     * X1 = 2.5;
     * X2 = 3.5;
     * Xu = 1.5; Location of the Ghost cell in fine block*/

    real X0 = 1.0;
    real X1 = 2.5;
    real X2 = 3.5;
    real Xu = 1.5;

    real L0 = ( Xu - X1 ) * ( Xu - X2 ) / ( ( X0 - X1 ) * ( X0 - X2 ) );

    real L1 = ( Xu - X0 ) * ( Xu - X2 ) / ( ( X1 - X0 ) * ( X1 - X2 ) );

    real L2 = ( Xu - X0 ) * ( Xu - X1 ) / ( ( X2 - X0 ) * ( X2 - X1 ) );

    // Later on (*) operator should be overloaded so that multiplication for other types happen automatically
    // Nvalue which is the solution vector contains other variables( u, v , w) and we don't want to write a sepparate line for each
    point[indexGhost].p = point[indexGhost].p * L0 + point[index1].p * L1 + point[index1].p * L2;

    cout << " QGhost is updated " << endl;
}

// interpolateGhostCells is a function that uses 3 points in the normal direction to the face of a block to
// create a quadratic function which can be used to find the right value for the guard(Ghost) cells.
// Because values send by the coarser blocks are all with double meshsizes (e.g. 2*dx, 2*dy, 2*dz)  wherease values in the finer blocks
// are half of the coarse blocks' mesh sizes (e.ge dx, dy, dz) . Therefore one needs to interpolate those values to find a new
// that matches the mesh size.
//
template <typename Nvalue>
void interpolate<Nvalue>::interpolateGhostCells( Nvalue *point, const int pxg, const int pyg, const int pzg, const int direction )
{
    // (0,0)Ghost Vale <-2dx-> node1 <- 3dx -> node2
    // Ghost Value(@dx) = function(Ghost Value, node1, node2);

    /* for face i = 0 */

    real X0 = 1.0;
    real X1 = 2.5;
    real X2 = 3.5;
    real Xu = 1.5;

    real L0 = ( Xu - X1 ) * ( Xu - X2 ) / ( ( X0 - X1 ) * ( X0 - X2 ) );

    real L1 = ( Xu - X0 ) * ( Xu - X2 ) / ( ( X1 - X0 ) * ( X1 - X2 ) );

    real L2 = ( Xu - X0 ) * ( Xu - X1 ) / ( ( X2 - X0 ) * ( X2 - X1 ) );

    // Later on (*) operator should be overloaded so that multiplication for other types happen automatically
    // // Nvalue which is the solution vector contains other variables( u, v , w) and we don't want to write a sepparate line for each
    // point[indexGhost].p = point[indexGhost].p*L0 + point[index1].p*L1 + point[index1].p*L2;

    switch ( direction )
    {
        // x direction
        case 0:

            for ( int k = 1; k < pzg - 1; k++ )
            {
                for ( int j = 1; j < pyg - 1; j++ )
                {
                    // we need i = 0 , i = 1, i = 2;

                    int indexGhost = pxg * pyg * k + pxg * j + 0;
                    int index1     = pxg * pyg * k + pxg * j + 1;
                    int index2     = pxg * pyg * k + pxg * j + 2;

                    // updates the point[indexGhost];
                    // quadraticOneD(point, indexGhost,index1,index2); instead of calling the function, do it here.
                    point[indexGhost] = point[indexGhost] * L0 + point[index1] * L1 + point[index2] * L2;
                }
            }

            /* for face i = pgx-1 */

            for ( int k = 1; k < pzg - 1; k++ )
            {
                for ( int j = 1; j < pyg - 1; j++ )
                {
                    // we need i = pxg-1, i = pxg-2, i = pxg-3;

                    int indexGhost = pxg * pyg * k + pxg * j + pxg - 1;
                    int index1     = pxg * pyg * k + pxg * j + pxg - 2;
                    int index2     = pxg * pyg * k + pxg * j + pxg - 3;

                    point[indexGhost] = point[indexGhost] * L0 + point[index1] * L1 + point[index2] * L2;
                }
            }

            break;

        // y direction
        case 1:

            /* for face j = 0 */

            for ( int k = 1; k < pzg - 1; k++ )
            {
                for ( int i = 1; i < pxg - 1; i++ )
                {
                    // we need j = 0, i = 1, i = 2;

                    int indexGhost = pxg * pyg * k + pxg * 0 + i;
                    int index1     = pxg * pyg * k + pxg * 1 + i;
                    int index2     = pxg * pyg * k + pxg * 2 + i;

                    point[indexGhost] = point[indexGhost] * L0 + point[index1] * L1 + point[index2] * L2;
                }
            }

            /* for face j = pyg-1 */

            for ( int k = 1; k < pzg - 1; k++ )
            {
                for ( int i = 1; i < pxg - 1; i++ )
                {
                    // we need j = pyg-1, i = pyg-2, i = pyg-3;

                    int indexGhost = pxg * pyg * k + pxg * ( pyg - 1 ) + i;
                    int index1     = pxg * pyg * k + pxg * ( pyg - 2 ) + i;
                    int index2     = pxg * pyg * k + pxg * ( pyg - 3 ) + i;

                    point[indexGhost] = point[indexGhost] * L0 + point[index1] * L1 + point[index2] * L2;
                }
            }

            break;

        // z direction
        case 2:

            /* for face k = 0 */

            for ( int j = 1; j < pyg - 1; j++ )
            {
                for ( int i = 1; i < pxg - 1; i++ )
                {
                    // we need k = 0, k = 1, k = 2;

                    int indexGhost = pxg * pyg * 0 + pxg * j + i;
                    int index1     = pxg * pyg * 1 + pxg * j + i;
                    int index2     = pxg * pyg * 2 + pxg * j + i;

                    point[indexGhost] = point[indexGhost] * L0 + point[index1] * L1 + point[index2] * L2;
                }
            }

            /* for face k = pzg -1 */

            for ( int j = 1; j < pyg - 1; j++ )
            {
                for ( int i = 1; i < pxg - 1; i++ )
                {
                    // we need k = pzg-1, k = pzg-2, k = pzg-3;

                    int indexGhost = pxg * pyg * ( pzg - 1 ) + pxg * j + i;
                    int index1     = pxg * pyg * ( pzg - 2 ) + pxg * j + i;
                    int index2     = pxg * pyg * ( pzg - 3 ) + pxg * j + i;

                    point[indexGhost] = point[indexGhost] * L0 + point[index1] * L1 + point[index2] * L2;
                }
            }

            break;
    }
}

template class interpolate<Q>;