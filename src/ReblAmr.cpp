#include "ReblAmr.h"
#include "definitions.h"
#include <type_traits>

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
ReblAmr<N, Nvalue, M, Mvalue>::ReblAmr( int argcs, char *pArgs[], real *length, real *coords, uint nx, uint ny,uint nz, real  *init_loc )
{
    // corrds is the center
    //
    // find the number of geometry files in the input and 
    // 
    proclevel = atoi( pArgs[argcs-2] );
    meshlevel = atoi( pArgs[argcs-1] );
  
    nSTL=argcs-3;

    cout<<" program name is "<<pArgs[0] <<endl; 
    cout<<" number of inputs "<<argcs <<endl; 
    
    cout<<" procLevel "<<proclevel <<endl; 
    cout<<" meshLevel "<<meshlevel <<endl; 
    cout<<" NumebrOfSTLFiles "<<nSTL <<endl; 

    xyz1[0] = coords[0] - length[0] * 0.5;
    xyz1[1] = coords[0] + length[0] * 0.5;

    xyz1[2] = coords[1] - length[1] * 0.5;
    xyz1[3] = coords[1] + length[1] * 0.5;

    xyz1[4] = coords[2] - length[2] * 0.5;
    xyz1[5] = coords[2] + length[2] * 0.5;

    cout << "proc level =" << proclevel << endl;
    cout << "mesh level =" << meshlevel << endl;


    int initFlag = 0;

    if ( MPI_Initialized( &initFlag ) != MPI_SUCCESS && Com.myrank == 0 )
    {
        cout << " Exit Code : " << ReblAmrGetErrorEnum( MPI_INIT_CHECK_FAIL ) << endl;
        exit( 1 );
    }

    if ( initFlag == 0 )
    {
        if ( MPI_Init( &argcs, &pArgs ) != MPI_SUCCESS && Com.myrank == 0 )
        {
            cout << " Exit Code : " << ReblAmrGetErrorEnum( MPI_INIT_FAIL ) << endl;
            exit( 1 );
        }
    }

    MPIStartUp();
//
// Multi  Geometrry is being added 
// Create Geometry object with the same size
//
  
   GMT0=new GeomSTL[nSTL];    

   real tmp[3];
 
   for(int i=0;i<nSTL;i++)  
   { 
     GMT0[i].construct(nSTL, xyz1 );
     tmp[0]=init_loc[i*3+0];
     tmp[1]=init_loc[i*3+1];
     tmp[2]=init_loc[i*3+2];
     GMT0[i].readSTLGeom( pArgs, i+1, xyz1,tmp );
    }

#if(1)

// old for single geometry
//   GMT.construct( xyz1 );
//   GMT.readSTLGeom( pArgs,1, xyz1 );

    Proc.construct( length, coords, nx, ny, nz );

    generateProcTopology();

    ID = new unsigned int[Proc.size()];

  if(ZOLTAN_ON)
  {
    Part.construct( argcs, pArgs, PART_METHOD, Proc.size() );
  }
 #endif 
  // uses distributed proc topology, without success in zoltan part this will fail
    //
}


template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::MPIStartUp()
{
    if ( MPI_SUCCESS != MPI_Comm_dup( MPI_COMM_WORLD, &Com.mpicom ) )
    {
        cout << BLUE << " Rank(" << Com.myrank << ") > Exit Code : " << MPI_DUP_FAIL << RESET << endl;
        cout << BLUE << ReblAmrGetErrorEnum( MPI_DUP_FAIL ) << RESET << endl;
        exit( 1 );
    }

    if ( MPI_SUCCESS != MPI_Comm_rank( Com.mpicom, &Com.myrank ) )
    {
        cout << BLUE << " Rank(" << Com.myrank << ") > Exit Code : " << MPI_GET_RANK_FAIL << RESET << endl;
        cout << BLUE << ReblAmrGetErrorEnum( MPI_GET_RANK_FAIL ) << RESET << endl;
        exit( 1 );
    }

    if ( MPI_SUCCESS != MPI_Comm_size( Com.mpicom, &Com.comsize ) )
    {
        cout << BLUE << " Rank(" << Com.myrank << ") > Exit Code : " << MPI_COMSIZE_FAIL << RESET << endl;
        cout << BLUE << ReblAmrGetErrorEnum( MPI_COMSIZE_FAIL ) << RESET << endl;
        exit( 1 );
    }
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
ReblAmr<N, Nvalue, M, Mvalue>::~ReblAmr()
{

 delete [] GMT0;

 if(ID!=nullptr)
  {
    delete[] ID;
  }
}


#if(1)
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::generateProcTopology()
{
    //**************************************************************
    //
    //                 generate processor topology
    //
    //**************************************************************
    double t1;
    MPI_Barrier( MPI_COMM_WORLD );
    t1 = MPI_Wtime();

  cout<<" procesor topology size "<<Proc.size()<<" proclevel  "<<proclevel<<endl;

  if(Proc.size()==1)
  {
    Proc.refine( 0 );
   }

    // get collective size
    //  


    for(int i=0;i<nSTL;i++)
    {
    Proc.convertStl2Morton( GMT0[i].geom_nn, GMT0[i].geom_xyz);
    }


    for ( uint j = 0; j < proclevel - 1; j++ )
    {
        Proc.pushToRefinelist( j + 1 );
        Proc.fourToOne();
        Proc.refineRefineList();
       //        cout<<j<<" : "<<Proc.size()<<endl;
    }

    double t2;
    MPI_Barrier( MPI_COMM_WORLD );
    t2 = MPI_Wtime();

cout<<" proc size "<<Proc.size()<<endl;
 
}
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::countPointsinBox()
{
    real      xyz[6];
    uint      count;
    morton<M> key;

    for ( auto it = Proc.begin(); it != Proc.end(); it++ )
    {
        // count is set to 1 such that an empty box has a weighting of one
        count = 1;

        key = it->first;
        Proc.enclosingBox( key, xyz );
        it->second = (uint *)malloc( sizeof( uint ) );

     for(uint k=0;k<nSTL;k++)
     {  
        for ( uint j = 0; j < GMT0[k].geom_nn; j++ )
        {
            if ( GMT0[k].geom_xyz[3 * j + 0] >= xyz[0] && GMT0[k].geom_xyz[3 * j + 0] <= xyz[1] )
            {
                if ( GMT0[k].geom_xyz[3 * j + 1] >= xyz[2] && GMT0[k].geom_xyz[3 * j + 1] <= xyz[3] )
                {
                    if ( GMT0[k].geom_xyz[3 * j + 2] >= xyz[4] && GMT0[k].geom_xyz[3 * j + 2] <= xyz[5] )
                    {
                        count++;
                    }
                }
            }
        }
     }
        it->second[0] = count;
    }

    /*
       for(auto it=proc.begin();it!=proc.end();it++)
       {
       count=0;
       key=it->first;
       cout<<"size "<<it->second[0]<<endl;
       }
      */
    // ???????????????????????????????????????????????????????????????????????????????????????
    // to be integrated soon
    /*????????????????????????????????????????????????????????????????????????????????????????

       if ( Proc.size() < comsize )
        {
            // cout<<proclevel<<" "<<treesize/3<<endl;
            throw std::runtime_error( RED "please increase topology level : Topolgy level is not distributable among the processes\n" RESET
    );
        }
    ??????????????????????????????????????????????????????????????????????????????????????????
    */
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::forestConstruct( int argcs, char *pArgs[], real *length, real *coords, uint nx, uint ny, uint nz )
{
    Forest.construct( argcs, pArgs, Proc, length, coords, nx, ny, nz );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::moveGeometry( double *xx)
{
   real tmp[3];
 
   for(int i=0;i<nSTL;i++)  
   { 
     tmp[0]=xx[i*3+0];
     tmp[1]=xx[i*3+1];
     tmp[2]=xx[i*3+2];
     Forest.moveGeom( Proc, GMT0[i].geom_xyz, GMT0[i].geom_nn, tmp );
  } 

}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::moveGeometry( double *xx,int *activeList)
{
    Forest.moveGeom( Proc, GMT0, nSTL, xx, activeList );
}

/*
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::moveGeometryDebug( double *xx )
{
    //Forest.moveGeomDebug( Proc, GMT.geom_xyz, GMT.geom_nn, xx );
}
*/

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::getTotalMeshSize()
{
    Forest.getTotalMeshSize();
}
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::refineForest()
{   //cout<< " ReblAmr.cpp refineForest -> Proc size : " <<Proc.size() <<endl;
    
    
	Forest.refineForestBalanced( meshlevel, Proc );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::assignPartToProc()
{
    //if ( std::is_same<P, Tree<N, Nvalue>>::value )
        int myvalue = Proc.size();
        int offset, totalvalue;
        totalvalue = Proc.size();

        // need to do comesize and comerank

        Zoltan_Out zoltan_out;

        Part.zoltanGeometricPartitionerSerial( myvalue, totalvalue, offset, Com.comsize, &zoltan_out );

        //        Part.zoltanGeometricPartitionerSerial( myvalue, totalvalue, offset, 1, zz, &zoltan_out );
        //    if ( my_rank == 0 )
        // non scalable part for zoltan partitioning

         // realloc the size 
         ID=(uint*)realloc(ID,myvalue*sizeof(uint));

        for ( uint i = 0; i < Proc.size(); i++ )
        {
            ID[i] = 0;
        }

        // elements to be exported are given,

        for ( int i = 0; i < zoltan_out.numExport; i++ )
        {
            ID[zoltan_out.exportLocalGids[i]] = zoltan_out.exportToPart[i];
            //      cout<<" dist  " << zoltan_out.exportLocalGids[i]<<endl;
        }

        int co = 0;

        for ( auto it = Proc.begin(); it != Proc.end(); it++ )
        {
            it->second[0] = ID[co];
            //                      cout <<"tag  "<< it->second[0] << endl;
            co++;
        }
   }

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::assignPartToProcDebug()
{
    //if ( std::is_same<P, Tree<N, Nvalue>>::value )
    
        int myvalue = Proc.size();
        int offset, totalvalue;
        totalvalue = Proc.size();

        // need to do comesize and comerank

        Zoltan_Out zoltan_out;

        Part.zoltanGeometricPartitionerSerial( myvalue, totalvalue, offset, Com.comsize, &zoltan_out );

        //  Part.zoltanGeometricPartitionerSerial( myvalue, totalvalue, offset, 1, zz, &zoltan_out );
        //  if ( my_rank == 0 )
        //  non scalable part for zoltan partitioning
        
         ID=(uint*)realloc(ID,myvalue*sizeof(uint));


        for ( uint i = 0; i < Proc.size(); i++ )
        {
            ID[i] = 0;
        }


        // elements to be exported are given,

        for ( int i = 0; i < zoltan_out.numExport; i++ )
        {
            ID[zoltan_out.exportLocalGids[i]] = zoltan_out.exportToPart[i];
            //      cout<<" dist  " << zoltan_out.exportLocalGids[i]<<endl;
        }

        int co = 0;

        for ( auto it = Proc.begin(); it != Proc.end(); it++ )
        {
            it->second[0] = ID[co];
            //                      cout <<"tag  "<< it->second[0] << endl;
            co++;
        }

}


template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::assignWeightsForPartDebug()
{
        real xyzc[3];
        uint co = 0;

        Part.reSize(Proc.size());

        for ( auto it = Proc.begin(); it != Proc.end(); it++ )
        {
            Proc.centroid( it->first, xyzc );
            Part.XYZ.push_back( CenterCoords() );
            Part.XYZ.at( co ).x = xyzc[0];
            Part.XYZ.at( co ).y = xyzc[1];
            Part.XYZ.at( co ).z = xyzc[2];

            // cout << " corrds " << Part.XYZ.at( co ).x << " " << Part.XYZ.at( co ).y << " " << Part.XYZ.at( co ).z << endl;
            if ( WEIGHT )
            {
                Part.weight[co] = it->second[0];
            }
            else
            {
                Part.weight[co] = 1.0;
            }
            // cout<<co<<"\t" <<Part.weight[co]<<endl;
            // ipecify rank
            it->second[0] = 0;
            co++;
        }
}


template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::assignWeightsForPart()
{
        real xyzc[3];
        uint co = 0;

// do not forget to resize the array for the new proc size
        Part.reSize(Proc.size());

        for ( auto it = Proc.begin(); it != Proc.end(); it++ )
        {
            Proc.centroid( it->first, xyzc );
            Part.XYZ.push_back( CenterCoords() );
            Part.XYZ.at( co ).x = xyzc[0];
            Part.XYZ.at( co ).y = xyzc[1];
            Part.XYZ.at( co ).z = xyzc[2];

            // cout << " corrds " << Part.XYZ.at( co ).x << " " << Part.XYZ.at( co ).y << " " << Part.XYZ.at( co ).z << endl;
            if ( WEIGHT )
            {
                Part.weight[co] = it->second[0];
            }
            else
            {
                Part.weight[co] = 1.0;
            }
            // cout<<co<<"\t" <<Part.weight[co]<<endl;
            // ipecify rank
            it->second[0] = 0;
            co++;
        }
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::distributeTopology()
{
  // cout<< std::is_same<P, Tree<N, Nvalue>>::value<<endl; 
   // if( std::is_same<P, Tree<N, Nvalue>>::value )

// cout<<" distributing " <<endl;
        countPointsinBox();
        assignWeightsForPart();
        assignPartToProc();
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::setForestParams()
{
    Forest.getMaxSeedsLevel( Proc );
    Forest.setMaxProcLevel( proclevel );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::createComPattern()
{
    setForestParams();

    Forest.comPatternConstruct( Proc );

    Forest.createCommGraph( 0 );

    Forest.checkGraphConsistency();
}
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::writeRunInfo()
{
    Forest.runInfo();
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::writeMesh( int index )
{
    //if ( std::is_same<P, Tree<N, uint>>::value )
    {
        templatePhdf5<N, Nvalue, M, Mvalue, Tree<M, Mvalue>> IO;

        if ( WR == 1 )
        {
            IO.writePolyvertex( Forest, index );
        }
        else if ( WR == 2 )
        {
            IO.writeMultiBlock( Forest, index );
        }
    }

}
#endif

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void  ReblAmr<N, Nvalue, M, Mvalue>::derefineMesh(int *activeList){
/* meshlevel is derefined from meshlevel-1 inside the pushToDerefine list */

 int nInactive=1;
#if(1) 
    real *parpd = new real[6*nInactive];
 
   int count = 0;

    for (int i = 0 ; i < nSTL; i++){
                    
        if(activeList[i]==0)
        {   
        for (int j = 0; j < 6; j++){

           parpd[count*6 + j] = GMT0[i].fitParlpd[j]; 

        }   
        
        count++;
       }   
    }   

 cout<<"parpd[0] = "<<parpd[0]<<" parpd[1] = "<<parpd[1] << " parpd[2] = "<<parpd[2]<<endl;
 cout<<"parpd[3] = "<<parpd[3]<<" parpd[4] = "<<parpd[4] << " parpd[5] = "<<parpd[5]<<endl;


#endif

 //   real parpd[6]={0,2.,-2.,2.,-2.,2.};


 
    for ( uint j = meshlevel ; j > 1 ; j-- )
    	{
 		//cout<<"ReblAmr.cpp  derefined mesh level is  = "<< j<<endl;
		//cout<<"ReblAmr.cpp  Proc size is  = "<<Proc.size()<<endl; 
		Forest.pushToDerefineEachTree(j,Proc,nInactive,parpd);
	}
	
	//delete [] parpd;

}


template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void  ReblAmr<N, Nvalue, M, Mvalue>::derefineMeshTopology(){
/* meshlevel is derefined from meshlevel-1 inside the pushToDerefine list */
 
    for ( uint j = meshlevel ; j > 1 ; j-- )
    	{
 		//cout<<"ReblAmr.cpp  derefined mesh level is  = "<< j<<endl;
		//cout<<"ReblAmr.cpp  Proc size is  = "<<Proc.size()<<endl; 
		Forest.pushToDerefineEachTree(j, Proc);
	}

}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void  ReblAmr<N, Nvalue, M, Mvalue>::updateGeom(real *x, int *activeList){

	//Forest.updateGeom(Proc, GMT.geom_xyz, GMT.geom_nn);

	Forest.updateGeom(Proc, GMT0, nSTL,x, activeList);

}
 
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void  ReblAmr<N, Nvalue, M, Mvalue>::getMeshLevel(){
 
	//	cout<<" Rebelamr.Cpp " << " Meshlevel = "<< meshlevel<<endl;
}
/*
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::selectKeysToIgnore(int *activeList)
{

//Forest.selectKeysToIgnore(Proc,activeList);

}
*/
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::derefineProcTopology(int *activeList)
{
//**************************************************************
//
//                 generate processor topology
//
//**************************************************************
//   added for multiple geometry 
//  define number of inactive lists in nInactive
//  allocate box with the size of 6 * nInactive 
//
    int nInactive=1;

	real *box = new real[6*nInactive];
#if(1)  
   int count = 0;

    for (int i = 0 ; i < nSTL; i++){
          		
        if(activeList[i]==0)
        {
		for (int j = 0; j < 6; j++){

	       box[count*6 + j] = GMT0[i].fitParlpd[j];	

		}
	
		count++;
       }
	}

 cout<<"box[0] = "<<box[0]<<" box[1] = "<<box[1] << " box[2] = "<<box[2]<<endl;
 cout<<"box[3] = "<<box[3]<<" box[4] = "<<box[4] << " box[5] = "<<box[5]<<endl;


#endif

   //real  box[6]={0,2.,-2.,2.,-2.,2.};

    double t1;
    MPI_Barrier( MPI_COMM_WORLD );
    t1 = MPI_Wtime();
    Forest.getMaxSeedsLevel( Proc );

   uint start =Forest.getMaxSeedLevel();
  cout<<" ReblAmr.cpp derefineProcTopology " << " max proc level = "<<start<<endl;
 	
	for ( uint j = start ;j > 1 ; j-- )
    {
 // uint j  = start;     
        Proc.pushToDerefinelist( j );

// 
// added to disable the derefinement for inactive list
//
      	//Proc.ignoreInactive(nInactive, box);

        Proc.ignoreInactiveVertices(nInactive, box);
       
        cout<<"level = "<<j<<"  Proc size before derefinement: "<<Proc.size()<<endl;
        Proc.retainFourToOne();
        Proc.derefineDerefineList();
		//cout<<" Proc is derefined  " << j<<" level"<<endl;
        cout<<"level = "<<j<<"  Proc size after derefinement: "<<Proc.size()<<endl;
  }

    cout<<" Final Proc Size "<<Proc.size()<<endl;

/* Here there is a problem with updating the list */
 
    Forest.updateSeeds(Proc);
    distributeTopology( );

    createComPattern();

    Forest.updateSeedsAndTrees(Proc);

    double t2;
    MPI_Barrier( MPI_COMM_WORLD );
    t2 = MPI_Wtime();

  for(auto it=Proc.begin();it!=Proc.end();it++)
  {
    //cout<<it->first<<" "<<it->second[0] << endl;
  }
delete [] box;
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::updateSeedsAndTrees()
{
    Forest.updateSeeds(Proc);

    distributeTopology( );

// for debug 
/*
     countPointsinBox();
     assignWeightsForPart();
     assignPartToProc();
*/
    createComPattern();
    Forest.updateSeedsAndTrees(Proc);

}


#if(0)
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::animateGeom(real *xyz , int index)
{
   	derefineMeshTopology(); 
    derefineProcTopology(); 
    moveGeometry(xyz); 
    
    distributeTopology();
    generateProcTopology();
    updateSeedsAndTrees();
    
    updateGeom();
    refineForest();
    writeMesh(index); 


}


template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::bounceOffGround(real dt, real *xx,real *length, real *coords, real Tfinal,real R_sphere,real zeta)
{

    int index        = 0;
    //real xx[3] = {0.0, 0.0, 0.0};
#if(1)

    real Time        = 0.0;                  // time
    //real dt          = 1.0;                  // time stepsize
    real radius      =  0.5;  
    real y_min       =  - length[1] * 0.5 + R_sphere;
    cout             << " y_min " << y_min << endl;
    real dy          = - 0.5;
    real y_top       = - coords[1];
    //real zeta        =   abs(dy);
    real y           = y_top;   
    cout<< " y_top = "<<y_top<<endl;
    cout<< " zeta  = "<< zeta<<endl;
	cout<<"y_min = "<<y_min<<endl;

while (Time < Tfinal)
{
    while( y > y_min && y <= y_top)
    {
    xx[1]       = dy;
    index       = index+1;  
    animateGeom(xx,index);
    cout << "y = "<< y << endl;
    y           = y + dy;   
}

#if(0)
   AMR.getTotalMeshSize();
   AMR.writeRunInfo();
#endif

 cout << "********* new **********"<<endl;

    if (dy < 0){
        y_top   = y_top - zeta;
        y       = y + abs(dy);
    }
    else {

         y       = y - abs(dy);
    }

    if( y_top <= y_min/2.0) {
        dy = dy*0.5;
	    y_min = y_min + abs(2*dy);
    }


    if(y_top <= y_min){

    break;                          
   }

    dy              = -dy;
    Time            = Time + dt;
    cout<<"y_top = " <<y_top<<endl;
	cout<<"y_min = "<< y_min<<endl;
    cout<<"dy = "<<dy<<endl; 
}

#endif 
}

#endif

template class ReblAmr<TREESIZE, real, PROCSIZE, uint>;
/*
template class ReblAmr<TREESIZE, real, PROCSIZE, uint, Tree<PROCSIZE, uint>>;
template class ReblAmr<TREESIZE, real, WSIZE, uint, FullTree<WSIZE, uint>>;
*/
