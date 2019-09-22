#include "solver.h"
#include "params.h"

/* These variables are to control the ghost cells of each directions
 * if they are set to 1, the ghost cells of that direction will be set
 * to the exact value of the function */

template <typename Nvalue>
void solver<Nvalue>::setBC(const char *bc, const real* Dirichlet, const real* Neumann )
{
    for ( int i = 0; i < 6; i++ )
    {
        this->bc[i] = bc[i];
		this->Dirichlet[i] = Dirichlet[i];
		this->Neumann[i] = Neumann[i];
    }
}

template <typename Nvalue>
void solver<Nvalue>::initializeTrigonometric( Nvalue *q, real *XYZ )
{
    real xmin, ymin, zmin, xmax, ymax, zmax;
    real x, y, z, delx, dely, delz;

    xmin = XYZ[0];
    xmax = XYZ[1];
    ymin = XYZ[2];
    ymax = XYZ[3];
    zmin = XYZ[4];
    zmax = XYZ[5];

    real freqScl = 0.5;
    // printf(" X[0] = %5.5lf ,X[1] = %5.5lf \n",xmin, xmax);
    // printf(" Y[0] = %5.4lf , Y[1] = %5.4lf \n",ymin, ymax);
    // printf(" Z[0] = %5.4lf , Z[1] = %5.4lf \n",zmin, zmax);

    delx = ( xmax - xmin ) / real( nx );
    dely = ( ymax - ymin ) / real( ny );
    delz = ( zmax - zmin ) / real( nz );

    // printf("dx =  %5.5lf \n",delx);
    // printf("dy =  %5.5lf \n",dely);
    // printf("dz =  %5.5lf \n",delz);

    const real pi = acos( -1.0 );
    // printf("pi = %4.5lf \n",pi);
    // cout<<nx<<" "<<ny<<" "<<nz<<endl;

    int  start = 0;
    int  end   = start;
    uint index;

    for ( uint k = start; k < nzg - end; k++ )
    {
        for ( uint j = start; j < nyg - end; j++ )
        {
            for ( uint i = start; i < nxg - end; i++ )
            {
                x = xmin + delx * real( i - start ) + delx * 0.5;
                y = ymin + dely * real( j - start ) + dely * 0.5;
                z = zmin + delz * real( k - start ) + delz * 0.5;

                index = ( nxg ) * (nyg)*k + (nxg)*j + i;

#if ( !POISSON )
                //		cout<<RED<<" pressure is initialized non-zero"<<RESET<<endl;
            /*    
			        q[index].p = sin( freqScl * OMEGA0 * pi * x - pi / 2 ) * sin( freqScl * OMEGA1 * pi * y - pi / 2 )
                                             * sin( freqScl * OMEGA2 * pi * z - pi / 2 );

        	*/
     
		  //  q[index].p = *XYZ;
			
			
		     q[index].p = i;
			 q[index].f = i+*XYZ;
					

#endif

#if ( POISSON )
                //  cout<<BLUE<<"poisson solution is invoked "<<RESET<<endl;
                q[index].p = 0.0;
			   
                
					q[index].f = cos( freqScl * OMEGA0 * pi * x -pi/2) * sin( freqScl * OMEGA1 * pi * y -pi/2)
                                             * sin( freqScl * OMEGA2 * pi *z -pi/2);
              // q[index].f = 0.0;

#endif
            }
        }
    }
}

/****************** updating the Ghost cells ************/

template <typename Nvalue>
void solver<Nvalue>::exactValueForGhost( Nvalue *q, real *XYZ )
{
    // cout<<"exactValueForGhost is called "<<endl;
    // Exact value of Ghost cells on Faces

    real xmin, ymin, zmin, xmax, ymax, zmax;
    real x, y, z, delx, dely, delz;

    xmin = XYZ[0];
    xmax = XYZ[1];
    ymin = XYZ[2];
    ymax = XYZ[3];
    zmin = XYZ[4];
    zmax = XYZ[5];

    real freqScl = 0.5;

    delx = ( xmax - xmin ) / real( nx );
    dely = ( ymax - ymin ) / real( ny );
    delz = ( zmax - zmin ) / real( nz );

    const real pi = acos( -1.0 );

    int start = 1;
    int end   = start;
    int index;

#if ( EXACTBC_Z && POISSON )

    // cout<<" Zdirection Ghost is set to Exact "<<endl;

    for ( int j = start; j < nyg - end; j++ )
    {
        for ( int i = start; i < nxg - end; i++ )
        {
            int k = 0;

            x = xmin + delx * real( i - start ) + delx * 0.5;
            y = ymin + dely * real( j - start ) + dely * 0.5;
            z = zmin - 0.5 * delz;

            index = ( nxg ) * (nyg)*k + (nxg)*j + i;

            q[index].p = sin( freqScl * OMEGA0 * pi * x - pi / 2 ) * sin( freqScl * OMEGA1 * pi * y - pi / 2 )
                         * sin( freqScl * OMEGA2 * pi * z - pi / 2 );
        }
    }

    for ( int j = start; j < nyg - end; j++ )
    {
        for ( int i = start; i < nxg - end; i++ )
        {
            int k = nzg - 1;

            x = xmin + delx * real( i - start ) + delx * 0.5;
            y = ymin + dely * real( j - start ) + dely * 0.5;
            z = zmax + 0.5 * delz;

            index = ( nxg ) * (nyg)*k + (nxg)*j + i;

            q[index].p = sin( freqScl * OMEGA0 * pi * x - pi / 2 ) * sin( freqScl * OMEGA1 * pi * y - pi / 2 )
                         * sin( freqScl * OMEGA2 * pi * z - pi / 2 );
        }
    }

#endif

#if ( EXACTBC_Y && POISSON )
    // cout<<" Ydirection Ghost is set to Exact "<<endl;

    for ( int k = start; k < nzg - end; k++ )
    {
        for ( int i = start; i < nxg - end; i++ )
        {
            int j = 0;

            x = xmin + delx * real( i - start ) + delx * 0.5;

            y = ymin - 0.5 * dely;

            z = zmin + delz * real( k - start ) + delz * 0.5;

            index = ( nxg ) * (nyg)*k + (nxg)*j + i;

            q[index].p = sin( freqScl * OMEGA0 * pi * x - pi / 2 ) * sin( freqScl * OMEGA1 * pi * y - pi / 2 )
                         * sin( freqScl * OMEGA2 * pi * z - pi / 2 );
        }
    }

    for ( int k = start; k < nzg - end; k++ )
    {
        for ( int i = start; i < nxg - end; i++ )
        {
            int j = nyg - 1;

            x = xmin + delx * real( i - start ) + delx * 0.5;
            y = ymax + 0.5 * dely;
            z = zmin + delz * real( k - start ) + delz * 0.5;

            index = ( nxg ) * (nyg)*k + (nxg)*j + i;

            q[index].p = sin( freqScl * OMEGA0 * pi * x - pi / 2 ) * sin( freqScl * OMEGA1 * pi * y - pi / 2 )
                         * sin( freqScl * OMEGA2 * pi * z - pi / 2 );
        }
    }

#endif

#if ( EXACTBC_X && POISSON )

    // cout<<" Xdirection Ghost is set to Exact "<<endl;

    for ( int k = start; k < nzg - end; k++ )
    {
        for ( int j = start; j < nyg - end; j++ )
        {
            int i = 0;

            x = xmin - 0.5 * delx;

            y = ymin + dely * real( j - start ) + dely * 0.5;

            z = zmin + delz * real( k - start ) + delz * 0.5;

            index = ( nxg ) * (nyg)*k + (nxg)*j + i;

            q[index].p = sin( freqScl * OMEGA0 * pi * x - pi / 2 ) * sin( freqScl * OMEGA1 * pi * y - pi / 2 )
                         * sin( freqScl * OMEGA2 * pi * z - pi / 2 );
        }
    }

    for ( int k = start; k < nzg - end; k++ )
    {
        for ( int j = start; j < nyg - end; j++ )
        {
            int i = nxg - 1;

            x = xmax + 0.5 * delx;

            y = ymin + dely * real( j - start ) + dely * 0.5;

            z = zmin + delz * real( k - start ) + delz * 0.5;

            index = ( nxg ) * (nyg)*k + (nxg)*j + i;

            q[index].p = sin( freqScl * OMEGA0 * pi * x - pi / 2 ) * sin( freqScl * OMEGA1 * pi * y - pi / 2 )
                         * sin( freqScl * OMEGA2 * pi * z - pi / 2 );
        }
    }

#endif
}

template <typename Nvalue>
void solver<Nvalue>::setGrid( uint n0, uint n1, uint n2 )
{
    nx = n0 - 1;
    ny = n1 - 1;
    nz = n2 - 1;

    nxg = n0 + 1;
    nyg = n1 + 1;
    nzg = n2 + 1;
}

template <typename Nvalue>
real solver<Nvalue>::redBlackGS(Nvalue *q, real *XYZ )
{ // red black Guss Siedel
	
	real xmin, ymin, zmin, xmax, ymax, zmax;
    real x, y, z, delx, dely, delz;
    real max = 0;
    xmin     = XYZ[0];
    xmax     = XYZ[1];
    ymin     = XYZ[2];
    ymax     = XYZ[3];
    zmin     = XYZ[4];
    zmax     = XYZ[5];


    delx = ( xmax - xmin ) / real( nx );
    dely = ( ymax - ymin ) / real( ny );
    delz = ( zmax - zmin ) / real( nz );

    real cft = 2.0 / ( delx * delx ) + 2.0 / ( dely * dely ) + 2.0 / ( delz * delz );
    real pi  = acos( -1.0 );

    double qtmp;
	// black nodes
	for ( uint k = 1; k < nzg - 1; k++ )
    {
        for ( uint j = 1; j < nyg - 1; j++ )
        {
            for ( uint i = 1; i < nxg - 1; i++ )
            {
                 
             if((i+j+k)%2 == 0) {
					uint ijk  = ( nxg ) * (nyg)*k + (nxg)*j + i;
                	uint imjk = ( nxg ) * (nyg)*k + (nxg)*j + ( i - 1 );
                	uint ipjk = ( nxg ) * (nyg)*k + (nxg)*j + ( i + 1 );
                	uint ijmk = ( nxg ) * (nyg)*k + ( nxg ) * ( j - 1 ) + i;
                	uint ijpk = ( nxg ) * (nyg)*k + ( nxg ) * ( j + 1 ) + i;
                	uint ijkm = ( nxg ) * ( nyg ) * ( k - 1 ) + (nxg)*j + i;
                	uint ijkp = ( nxg ) * ( nyg ) * ( k + 1 ) + (nxg)*j + i;
	
                	qtmp = q[ijk].p;
                	q[ijk].p = ( ( q[imjk].p + q[ipjk].p ) / ( delx * delx ) + ( q[ijmk].p + q[ijpk].p ) / ( dely * dely )
                             + ( q[ijkm].p + q[ijkp].p ) / ( delz * delz ) - q[ijk].f )
                           / cft;
                	max = fmax( max, fabs( q[ijk].p - qtmp ) );
           		}
	         }
        }
    }

	// red nodes 
	for ( uint k = 1; k < nzg - 1; k++ )
    {
        for ( uint j = 1; j < nyg - 1; j++ )
        {
            for ( uint i = 1; i < nxg - 1; i++ )
            {
                if((i+j+k)%2 == 1) {
					uint ijk  = ( nxg ) * (nyg)*k + (nxg)*j + i;
                	uint imjk = ( nxg ) * (nyg)*k + (nxg)*j + ( i - 1 );
                	uint ipjk = ( nxg ) * (nyg)*k + (nxg)*j + ( i + 1 );
                	uint ijmk = ( nxg ) * (nyg)*k + ( nxg ) * ( j - 1 ) + i;
                	uint ijpk = ( nxg ) * (nyg)*k + ( nxg ) * ( j + 1 ) + i;
                	uint ijkm = ( nxg ) * ( nyg ) * ( k - 1 ) + (nxg)*j + i;
                	uint ijkp = ( nxg ) * ( nyg ) * ( k + 1 ) + (nxg)*j + i;
	
                	qtmp = q[ijk].p;
                	q[ijk].p = ( ( q[imjk].p + q[ipjk].p ) / ( delx * delx ) + ( q[ijmk].p + q[ijpk].p ) / ( dely * dely )
                             + ( q[ijkm].p + q[ijkp].p ) / ( delz * delz ) - q[ijk].f )
                           / cft;
                	max = fmax( max, fabs( q[ijk].p - qtmp ) );
           		}
			 }
        }
    }
    return max;

} 



template <typename Nvalue>
real solver<Nvalue>::updateInterior( Nvalue *q, real *XYZ )
{
    real xmin, ymin, zmin, xmax, ymax, zmax;
    real x, y, z, delx, dely, delz;
    real max = 0;
    xmin     = XYZ[0];
    xmax     = XYZ[1];
    ymin     = XYZ[2];
    ymax     = XYZ[3];
    zmin     = XYZ[4];
    zmax     = XYZ[5];

    // printf(" X[0] = %5.5lf ,X[1] = %5.5lf \n",xmin, xmax);
    // printf(" Y[0] = %5.4lf , Y[1] = %5.4lf \n",ymin, ymax);
    // printf(" Z[0] = %5.4lf , Z[1] = %5.4lf \n",zmin, zmax);

    delx = ( xmax - xmin ) / real( nx );
    dely = ( ymax - ymin ) / real( ny );
    delz = ( zmax - zmin ) / real( nz );

    real cft = 2.0 / ( delx * delx ) + 2.0 / ( dely * dely ) + 2.0 / ( delz * delz );
    real pi  = acos( -1.0 );

    double qtmp;
    // uint nLoops = 1000 ;

    // for (uint l = 1; l < nLoops; l++){
if(GS == 1) 
{
    for ( uint k = 1; k < nzg - 1; k++ )
    {
        for ( uint j = 1; j < nyg - 1; j++ )
        {
            for ( uint i = 1; i < nxg - 1; i++ )
            {
                uint ijk  = ( nxg ) * (nyg)*k + (nxg)*j + i;
                uint imjk = ( nxg ) * (nyg)*k + (nxg)*j + ( i - 1 );
                uint ipjk = ( nxg ) * (nyg)*k + (nxg)*j + ( i + 1 );
                uint ijmk = ( nxg ) * (nyg)*k + ( nxg ) * ( j - 1 ) + i;
                uint ijpk = ( nxg ) * (nyg)*k + ( nxg ) * ( j + 1 ) + i;
                uint ijkm = ( nxg ) * ( nyg ) * ( k - 1 ) + (nxg)*j + i;
                uint ijkp = ( nxg ) * ( nyg ) * ( k + 1 ) + (nxg)*j + i;

                qtmp = q[ijk].p;

                q[ijk].p = ( ( q[imjk].p + q[ipjk].p ) / ( delx * delx ) + ( q[ijmk].p + q[ijpk].p ) / ( dely * dely )
                             + ( q[ijkm].p + q[ijkp].p ) / ( delz * delz ) - q[ijk].f )
                           / cft;
                //	q[ijk].p = 100;
                //      q[ijk].p = q[ijk].f;
                //	cout<<" x "<<x<<" "<<y<<" "<<z <<" "<<q[ijk].p<<endl;
                max = fmax( max, fabs( q[ijk].p - qtmp ) );
            }
        }
    }
}
else 
{
	max = redBlackGS(q,XYZ);
}

    return max;
}


template <typename Nvalue>
real solver<Nvalue>::getGradient( Nvalue *q, real *XYZ )
{
    real xmin, ymin, zmin, xmax, ymax, zmax;
    real delx, dely, delz;
    // calculate gradient and average it out on the entire domain, then normalize it
    real grad = 0.0;
    int  index0;
    int  index1;
    int  index2;
    int  index3;
    int  index4;
    int  index5;
    int  index6;

    xmin = XYZ[0];
    xmax = XYZ[1];
    ymin = XYZ[2];
    ymax = XYZ[3];
    zmin = XYZ[4];
    zmax = XYZ[5];

    delx = ( xmax - xmin ) / real( nx );
    dely = ( ymax - ymin ) / real( ny );
    delz = ( zmax - zmin ) / real( nz );

    for ( uint k = 1; k < nzg - 1; k++ )
    {
        for ( uint j = 1; j < nyg - 1; j++ )
        {
            for ( uint i = 1; i < nxg - 1; i++ )
            {
                //
                index0 = ( nxg ) * (nyg)*k + (nxg)*j + i + 1;

                index1 = ( nxg ) * (nyg)*k + (nxg)*j + ( i + 1 ) + 1;
                index2 = ( nxg ) * (nyg)*k + (nxg)*j + ( i - 1 ) + 1;

                index3 = ( nxg ) * (nyg)*k + ( nxg ) * ( j + 1 ) + i + 1;
                index4 = ( nxg ) * (nyg)*k + ( nxg ) * ( j - 1 ) + i + 1;

                index5 = ( nxg ) * ( nyg ) * ( k + 1 ) + (nxg)*j + i + 1;
                index6 = ( nxg ) * ( nyg ) * ( k - 1 ) + (nxg)*j + i + 1;

                grad += ( fabs( ( ( q[index1].p - q[index0].p ) ) / delx ) + fabs( ( ( q[index3].p - q[index0].p ) ) / dely )
                          + fabs( ( ( q[index5].p - q[index0].p ) ) / delz ) );

                //  cout<<" x "<<x<<" "<<y<<" "<<z <<" "<<q[index]<<endl;
            }
        }
    }

    cout << " grad " << grad / nxg / nyg / nxg << " delxs " << delx << " " << dely << " " << delz << endl;
    return ( grad / ( nxg * nyg * nzg ) );
}

template <typename Nvalue>
void solver<Nvalue>::imposeBoundaryXdirection( Nvalue *point, int faceTag, real* XYZ)
{

// impose X boundary condition
    real xmin, ymin, zmin, xmax, ymax, zmax;
    real x, y, z, delx, dely, delz;

    xmin = XYZ[0];
    xmax = XYZ[1];
    ymin = XYZ[2];
    ymax = XYZ[3];
    zmin = XYZ[4];
    zmax = XYZ[5];

    delx = ( xmax - xmin ) / real( nx );
    dely = ( ymax - ymin ) / real( ny );
    delz = ( zmax - zmin ) / real( nz );

    static const real pi = acos( -1.0 ); // make it a data member
    double omega[3] = {OMEGA0 * pi, OMEGA1 * pi, OMEGA2 * pi}; // make it a data member - data hiding and ADT
    int start = 1;
    int end   = start;


    if ( faceTag == 0 ) // left face
    {
      for ( int k = start; k < nzg - end; k++ )
        {
          for ( int j = start; j < nyg - end; j++ )
           {
		     // check if the left side has Dirichlet boundary condition
		     if(bc[0] == 'D') {
 
                int i = 0;

                x = xmin - 0.5 * delx;

                y = ymin + dely * real( j - start ) + dely * 0.5;

                z = zmin + delz * real( k - start ) + delz * 0.5;

                int indexGhost = ( nxg ) * (nyg)*k + (nxg)*j + i;
			    // initialize the boundaries with the exact equation for calculating the order of convergence using the MMS method
	    	    point[indexGhost].p = exactValue(omega[0], x, tags[0])*exactValue(omega[1],y,tags[1])* exactValue(omega[2],z,tags[2]);

#if(TRNS_INTRP_TEST)
		
				point[indexGhost].p = Dirichlet[0]; 
#endif

		      }  
			  // check if the left face has Neumann boundary condition
		      if(bc[0] == 'N') {

               	 int indexGhost      = nxg * nyg * k + nxg * j + 0;
			  	 int indexInnerFace  = nxg * nyg * k + nxg * j + 2;
			 	 // Neumann Boundary condition
			     point[indexGhost].p  = point[indexInnerFace].p - Neumann[0]*2.0*delx;
	
			  }
           
		    }
        }

      } else if ( faceTag == 1 ) 
	  { // right face
     
		for ( int k = start; k < nzg - end; k++ )
     	{
          for ( int j = start; j < nyg - end; j++ )
       	  {
			// check if the right face Dirichlet boundary condition
	   		if(bc[1] == 'D') {
 
             int i = nxg-1;

             x = xmax + 0.5 * delx; // x coordinate of the most right Ghost cells 

             y = ymin + dely * real( j - start ) + dely * 0.5;

             z = zmin + delz * real( k - start ) + delz * 0.5;

             int indexGhost = ( nxg ) * (nyg)*k + (nxg)*j + i;

			 // initialize the boundaries with the exact equation for calculating the order of convergence using the MMS method
	    	 point[indexGhost].p = exactValue(omega[0], x, tags[0])*exactValue(omega[1],y,tags[1])* exactValue(omega[2],z,tags[2]);

#if(TRNS_INTRP_TEST)				
			 point[indexGhost].p = Dirichlet[1]; 
#endif
	
		   }
		   // check if the right face has Neumann boundary condition
		   if(bc[1] == 'N') {

             int indexGhost      = nxg * nyg * k + nxg * j + nxg-1;
			 int indexInnerFace  = nxg * nyg * k + nxg * j + nxg-3;
			 // Neumann Boundary condition
			 point[indexGhost].p  = point[indexInnerFace].p + Neumann[1]*2.0*delx;
		
			}
           
		   }
        }
	 }

// end 
}

template <typename Nvalue>
void solver<Nvalue>::imposeBoundaryYdirection( Nvalue *point, int faceTag, real*XYZ )
{
// impose Y direction boundary conditions 

    real xmin, ymin, zmin, xmax, ymax, zmax;
    real x, y, z, delx, dely, delz;

    xmin = XYZ[0];
    xmax = XYZ[1];
    ymin = XYZ[2];
    ymax = XYZ[3];
    zmin = XYZ[4];
    zmax = XYZ[5];

    delx = ( xmax - xmin ) / real( nx );
    dely = ( ymax - ymin ) / real( ny );
    delz = ( zmax - zmin ) / real( nz );

    static const real pi = acos( -1.0 ); // make it a data member
    double omega[3] = {OMEGA0 * pi, OMEGA1 * pi, OMEGA2 * pi}; // make it a data member - data hiding and ADT
    int start = 1;
    int end   = start;


    if ( faceTag == 0 ) // left face
    {
      for ( int k = start; k < nzg - end; k++ )
        {
          for ( int i = start; i < nxg - end; i++ )
           {
		     // check if the left side has Dirichlet boundary condition
		     if(bc[2] == 'D') {
 
                int j = 0;
            	x = xmin + delx * real( i - start ) + delx * 0.5;
        	    y = ymin - 0.5 * dely;
     	       	z = zmin + delz * real( k - start ) + delz * 0.5;
                int indexGhost = ( nxg ) * (nyg)*k + (nxg)*j + i;
			    // initialize the boundaries with the exact equation for calculating the order of convergence using the MMS method
	    	    point[indexGhost].p = exactValue(omega[0], x, tags[0])*exactValue(omega[1],y,tags[1])* exactValue(omega[2],z,tags[2]);		
#if(TRNS_INTRP_TEST)		     
			 point[indexGhost].p = Dirichlet[2]; 
#endif

				 }  
			  // check if the left face has Neumann boundary condition
		      if(bc[2] == 'N') {

               	 int indexGhost      = nxg * nyg * k + nxg * 0 + i;
			  	 int indexInnerFace  = nxg * nyg * k + nxg * 2 + i;
			 	 // Neumann Boundary condition
			     point[indexGhost].p  = point[indexInnerFace].p - Neumann[2]*2.0*dely;
		
			  }
           
		    }
        }

      } else if ( faceTag == 1 ) 
	  { // right face
      for ( int k = start; k < nzg - end; k++ )
        {
          for ( int i = start; i < nxg - end; i++ )
           {
		     // check if the left side has Dirichlet boundary condition
		     if(bc[3] == 'D') {
 
                int j = nyg-1;
            	x = xmin + delx * real( i - start ) + delx * 0.5;
        	    y = ymax + 0.5 * dely;
     	       	z = zmin + delz * real( k - start ) + delz * 0.5;
                int indexGhost = ( nxg ) * (nyg)*k + (nxg)*j + i;
			    // initialize the boundaries with the exact equation for calculating the order of convergence using the MMS method
	    	    point[indexGhost].p = exactValue(omega[0], x, tags[0])*exactValue(omega[1],y,tags[1])* exactValue(omega[2],z,tags[2]);
#if(TRNS_INTRP_TEST) 	
			 cout<<" boundary is not for MMS"<<endl;	
			 point[indexGhost].p = Dirichlet[3]; 
#endif

		      }  
			  // check if the left face has Neumann boundary condition
		      if(bc[3] == 'N') {

               	 int indexGhost      = nxg * nyg * k + nxg * (nyg-1) + i;
			  	 int indexInnerFace  = nxg * nyg * k + nxg * (nyg-3) + i;
			 	 // Neumann Boundary condition
			     point[indexGhost].p  = point[indexInnerFace].p + Neumann[3]*2.0*dely;
		
			  }
           
		    }
        }

    }

// end 
}
  


template <typename Nvalue>
void solver<Nvalue>::imposeBoundaryZdirection( Nvalue *point, int faceTag, real*XYZ )
{
// impose Z direction boundary conditions 


    real xmin, ymin, zmin, xmax, ymax, zmax;
    real x, y, z, delx, dely, delz;

    xmin = XYZ[0];
    xmax = XYZ[1];
    ymin = XYZ[2];
    ymax = XYZ[3];
    zmin = XYZ[4];
    zmax = XYZ[5];

    delx = ( xmax - xmin ) / real( nx );
    dely = ( ymax - ymin ) / real( ny );
    delz = ( zmax - zmin ) / real( nz );

    static const real pi = acos( -1.0 ); // make it a data member
    double omega[3] = {OMEGA0 * pi, OMEGA1 * pi, OMEGA2 * pi}; // make it a data member - data hiding and ADT
    int start = 1;
    int end   = start;


    if ( faceTag == 0 ) // left face
    {
 
     	for ( int j = start; j < nyg - end; j++ )
    	{
        	for ( int i = start; i < nxg - end; i++ )
        	{
   
		     // check if the left side has Dirichlet boundary condition
		     if(bc[4] == 'D') {

	         	int k = 0;
            	x = xmin + delx * real( i - start ) + delx * 0.5;
            	y = ymin + dely * real( j - start ) + dely * 0.5;
            	z = zmin - 0.5 * delz;
				int indexGhost = ( nxg ) * (nyg)*k + (nxg)*j + i;
			    // initialize the boundaries with the exact equation for calculating the order of convergence using the MMS method
	    	    point[indexGhost].p = exactValue(omega[0], x, tags[0])*exactValue(omega[1],y,tags[1])* exactValue(omega[2],z,tags[2]);
			 }
			  // check if the left face has Neumann boundary condition
		      if(bc[4] == 'N') {

               	 int indexGhost      = nxg * nyg * 0 + nxg * j + i;
			  	 int indexInnerFace  = nxg * nyg * 2 + nxg * j + i;
			 	 // Neumann Boundary condition
			     point[indexGhost].p  = point[indexInnerFace].p - Neumann[4]*2.0*delz;
		
			  }

        	}
    	}

	} else if (faceTag == 1 ) 
	
	{
    	for ( int j = start; j < nyg - end; j++ )
    	{
        	for ( int i = start; i < nxg - end; i++ )
       		
			 {
   
		     // check if the left side has Dirichlet boundary condition
		     if(bc[5] == 'D') {

	         	int k = nzg-1;
            	x = xmin + delx * real( i - start ) + delx * 0.5;
            	y = ymin + dely * real( j - start ) + dely * 0.5;
            	z = zmax + 0.5 * delz;
				int indexGhost = ( nxg ) * (nyg)*k + (nxg)*j + i;
			    // initialize the boundaries with the exact equation for calculating the order of convergence using the MMS method
	    	    point[indexGhost].p = exactValue(omega[0], x, tags[0])*exactValue(omega[1],y,tags[1])* exactValue(omega[2],z,tags[2]);
			 }
			  // check if the left face has Neumann boundary condition
		      if(bc[5] == 'N') {

               	 int indexGhost      = nxg * nyg * (nzg-1) + nxg * j + i;
			  	 int indexInnerFace  = nxg * nyg * (nzg-3) + nxg * j + i;
			 	 // Neumann Boundary condition
			     point[indexGhost].p  = point[indexInnerFace].p + Neumann[5]*2.0*delz;
	

			  }  
     
		 	}
	   	
		} 	
	}

// end 
}


template <typename Nvalue>
void solver<Nvalue>::initializeTrigonometricForMMS( Nvalue *q, real *XYZ )
{
    real xmin, ymin, zmin, xmax, ymax, zmax;
    real x, y, z, delx, dely, delz;

    static const real pi       = acos( -1.0 );
    double            omega[3] = {OMEGA0 * pi, OMEGA1 * pi, OMEGA2 * pi};

    xmin = XYZ[0];
    xmax = XYZ[1];
    ymin = XYZ[2];
    ymax = XYZ[3];
    zmin = XYZ[4];
    zmax = XYZ[5];

    // printf(" X[0] = %5.5lf ,X[1] = %5.5lf \n",xmin, xmax);
    // printf(" Y[0] = %5.4lf , Y[1] = %5.4lf \n",ymin, ymax);
    // printf(" Z[0] = %5.4lf , Z[1] = %5.4lf \n",zmin, zmax);

    delx = ( xmax - xmin ) / real( nx );
    dely = ( ymax - ymin ) / real( ny );
    delz = ( zmax - zmin ) / real( nz );

    // printf("dx =  %5.5lf \n",delx);
    // printf("dy =  %5.5lf \n",dely);
    // printf("dz =  %5.5lf \n",delz);

    // printf("pi = %4.5lf \n",pi);
    // cout<<nx<<" "<<ny<<" "<<nz<<endl;

	// set the tags based on the boundary conditions
	extractTag();

    int  start = 1;
    int  end   = start;
    uint index;

    for ( uint k = start; k < nzg - end; k++ )
    {
        for ( uint j = start; j < nyg - end; j++ )
        {
            for ( uint i = start; i < nxg - end; i++ )
            {
                x = xmin + delx * real( i - start ) + delx * 0.5;
                y = ymin + dely * real( j - start ) + dely * 0.5;
                z = zmin + delz * real( k - start ) + delz * 0.5;

                index = ( nxg ) * (nyg)*k + (nxg)*j + i;

                q[index].f = ( omega[0] * omega[0] ) * exactValue( omega[0], x, tags[0] ) * exactValue( omega[1], y, tags[1] )
                             * exactValue( omega[2], z, tags[2] )
                             + ( omega[1] * omega[1] ) * exactValue( omega[0], x, tags[0] ) * exactValue( omega[1], y, tags[1] )
                               * exactValue( omega[2], z, tags[2] )
                             + ( omega[2] * omega[2] ) * exactValue( omega[0], x, tags[0] ) * exactValue( omega[1], y, tags[1] )
                               * exactValue( omega[2], z, tags[2] );
                               
                q[index].f = -q[index].f;
				q[index].p = 0.0;
            }
        }
    }
}
// Dirichlet for poisson is rather theoretical for fluid flow, we are more interested in
// Neumann and Periodic, for now we are not focusing on periodic
// write a method here
// to assign the Dirichlet Boundary condition,
// if tags[i]==0 (for sin function where we have Dirichlet BC)
// assign boundary values
// 
//  q[index].p = exactValue( omega[0], x, tags[0] ) * exactValue( omega[1], y, tags[1] ) * exactValue( omega[2], z, tags[2] );
// where the values are evaluated at the ghost cell's center 


// added for MMS
template <typename Nvalue>
real solver<Nvalue>::exactValue( double omega, double x, int tag )
{
    double result = 0.0;

    if ( tag == 0 )
    {
        result = sin( omega * x );
    }
    else if ( tag == 1 )
    {
        result = cos( omega * x );
    }
    else
    {
     //   cout << " only Drichlet and Neumann are considered " << endl;
        exit( 0 );
    }

    return ( result );
}

template <typename Nvalue>
void solver<Nvalue>::extractTag()
{
    for ( int i = 0; i < 3; i++ )
    {
        if ( bc[2 * i] == 'D' && bc[2 * i + 1] == 'D' )
        {
            tags[i] = 0;
        }
        else if ( bc[2 * i] == 'N' && bc[2 * i + 1] == 'N' )
        {
            tags[i] = 1;
        }
    }
    //std::cout << "tags " << tags[0] << " " << tags[1] << " " << tags[2] << std::endl;
}

template <typename Nvalue>
real solver<Nvalue>::getError( Nvalue *q, real *XYZ )
{
    real xmin, ymin, zmin, xmax, ymax, zmax;
    real x, y, z, delx, dely, delz;

    static const real pi       = acos( -1.0 );
    double            omega[3] = {OMEGA0 * pi, OMEGA1 * pi, OMEGA2 * pi};

    xmin = XYZ[0];
    xmax = XYZ[1];
    ymin = XYZ[2];
    ymax = XYZ[3];
    zmin = XYZ[4];
    zmax = XYZ[5];

    // printf(" X[0] = %5.5lf ,X[1] = %5.5lf \n",xmin, xmax);
    // printf(" Y[0] = %5.4lf , Y[1] = %5.4lf \n",ymin, ymax);
    // printf(" Z[0] = %5.4lf , Z[1] = %5.4lf \n",zmin, zmax);

    delx = ( xmax - xmin ) / real( nx );
    dely = ( ymax - ymin ) / real( ny );
    delz = ( zmax - zmin ) / real( nz );

    // printf("dx =  %5.5lf \n",delx);
    // printf("dy =  %5.5lf \n",dely);
    // printf("dz =  %5.5lf \n",delz);

    // printf("pi = %4.5lf \n",pi);
    // cout<<nx<<" "<<ny<<" "<<nz<<endl;

    int  start = 1;
    int  end   = start;
    uint index;
    real max=0.0;
    real diff;

    for ( uint k = start; k < nzg - end; k++ )
    {
        for ( uint j = start; j < nyg - end; j++ )
        {
            for ( uint i = start; i < nxg - end; i++ )
            {
                x = xmin + delx * real( i - start ) + delx * 0.5;
                y = ymin + dely * real( j - start ) + dely * 0.5;
                z = zmin + delz * real( k - start ) + delz * 0.5;

                index = ( nxg ) * (nyg)*k + (nxg)*j + i;

                diff=q[index].p - exactValue( omega[0], x, tags[0] ) * exactValue( omega[1], y, tags[1] )* exactValue( omega[2], z, tags[2] );
// l-inf norm              
                max=fmax(fabs(diff),max);

            }
        }
    }

return(max);

}



template class solver<Q>;
