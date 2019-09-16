#include "ReblAmr.h"
#include "ReblAmrFull.h"

/** \mainpage
 *
 *
 *  <STRONG>
 *   "rebl-AMR: An Open-Source Software for Binarized-Octree Mesh Generation around Immersed Geometries"
 * </STRONG>
 *
 * \details
 * This software generates forest of octrees for an immersed geometry. <br>
 * reble-AMR uses red-black tree data strcture to store the nodes in a Z-curve compliant fashion. <br>
 * This project is part of the <STRONG>NSF-GEM3D </STRONG>Award No. 1440638. <br>
 *
 *	Required Libraries:
 *
 *      (1)  CMAKE
 *      (2)  MPI (MPI 3.0 standard compliant version)
 *      (3)  Zoltan
 *      (4)  ParMetis
 *      (5)  HDF5
 *
 *      Usage:
 *      mpirun -np <number of processes> progName input/myGeometry.stl <params.txt
 *
 * @authors     Jaber J. Hasbestan and Inanc Senocak (PI) <br>
 *
 * @date        Feb 2018
 *
 * @copyright  (c)
 * Department of Mechanical Engineering and Materials Science <br>
 * University of Pittsburgh, Pittsburgh, PA <br>
 * <br>
 *
 * Please refere to either of the JCP papers \cite Hasbestan2018 and \cite Hasbestan2017 listed below if you are using this code
 *
 * @ref  "[1] Jaber J. Hasbestan and Inanc Senocak. Binarized-octree generation for cartesian adaptive mesh refinement around immersed geometries. Journal of Computational Physics, 368:179–195, sep 2018"
 *
 * @ref  "[2] Jaber J. Hasbestan and Inanc Senocak. A short note on the use of the red textendash black tree in cartesian adaptive mesh refinement algorithms. Journal of Computational Physics, 351:473–477, dec 2017" 
 *  
 *
 *  \image html f16_bunny.png
 *  
 *   
 *  <STRONG>
 *	Subject to GPL 3.0 license
 *  </STRONG>
 *
*/

int main( int argcs, char *pArgs[] )
{
 
   int proclevel = atoi (pArgs[2]);
   int meshlevel = atoi (pArgs[3]);

   if(argcs!=4)
   {
    cout<< " need three input vars " <<endl;
    cout<< " STL geometry " <<endl;
    cout<< " proclevel " <<endl;
    cout<< " mesh level " <<endl;
    }

    int nSTL=argcs-3;

    unsigned int i, j, k, l;
    int my_rank, comsize;

    real ancestorlength[3] = {4.f, 4.f, 4.f};
    real ancestorcoords[3] = {0.f,0.0, 0.f};

    real *init_loc=new real[3* (nSTL)];

// set the initial locations of the geometries. 
// 3*i+j  
//

    init_loc[0]=1.0;
    init_loc[1]=0.0;
    init_loc[2]=0.0;

if(nSTL>1)
{ 
    init_loc[3]=0.0;
    init_loc[4]=1.0;
    init_loc[5]=0.0;
}

if(nSTL>2)
{
    init_loc[6]=0.0;
    init_loc[7]=0.5;
    init_loc[8]=0.0;
}

if(nSTL>3)
{
    init_loc[9]=0.0;
    init_loc[10]=-0.5;
    init_loc[11]=0.0;
}

if(nSTL>4)
{
    init_loc[12]=0.0;
    init_loc[13]=0.0;
    init_loc[14]=0.5;
}

if(nSTL>5)
{
    init_loc[15]=0.0;
    init_loc[16]=0.0;
    init_loc[17]=-0.5;
}


    ReblAmr<TREESIZE, real, PROCSIZE, uint> AMR( argcs, pArgs, ancestorlength, ancestorcoords, npx, npy, npz,init_loc);

    AMR.distributeTopology( );

    AMR.forestConstruct(argcs,pArgs,ancestorlength, ancestorcoords, npx,npy, npz );

    AMR.createComPattern( );   

//    AMR.writeMesh(0);
    int index  		= 0;

    real *xx=new real[3* nSTL];
 

    memset(xx,0.0,3*(nSTL)*sizeof(double));

    cout<< nSTL <<endl;

    int activeList[nSTL];
    
    memset(activeList,1,(nSTL)*sizeof(int));

    AMR.moveGeometry(xx,activeList);
    AMR.refineForest();
    AMR.writeMesh(index);
    AMR.writeRunInfo();

	xx[4] = -1.0;
    index = 1;
	activeList[0] = 0;
	AMR.derefineMesh(activeList); 
    AMR.derefineProcTopology(activeList);
   
	AMR.moveGeometry(xx,activeList);
	AMR.distributeTopology();
	AMR.generateProcTopology();
	AMR.updateSeedsAndTrees();
	AMR.updateGeom(xx, activeList); 

    AMR.refineForest();	

	AMR.writeMesh(index);
    AMR.writeRunInfo();





#if(0)
 	real dy 		 = -0.25;
    xx[1]			 = dy;
	real dt 		 = 1.0;	 			 	 // time stepsize
	real Tsimulation = 30.0;	 			 // simulation time 
   	real radius      = 0.5;  
    real zeta 		 = abs(dy);
	AMR.bounceOffGround(dt,xx,ancestorlength,ancestorcoords,Tsimulation,radius,zeta);

#endif

 delete[] init_loc;

 delete[] xx;

}

