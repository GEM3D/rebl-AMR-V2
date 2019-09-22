#include "ReblAmr.h"
#include "definitions.h"
#include <type_traits>

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
ReblAmr<N, Nvalue, M, Mvalue>::ReblAmr( int argcs, char *pArgs[], real *length, real *coords, uint nx, uint ny, uint nz )
{
    // corrds is the center
    //

    proclevel = atoi( pArgs[2] );
    meshlevel = atoi( pArgs[3] );

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

    GMT.construct( xyz1 );

    GMT.readSTLGeom( pArgs, xyz1 );

    Proc.construct( length, coords, nx, ny, nz );

    //  FullProc.construct(length,coords, nx, ny, nz);


    generateProcTopology();

    ID = new unsigned int[Proc.size()];

   if(ZOLTAN_ON)
   {
    Part.construct( argcs, pArgs, PART_METHOD, Proc.size() );
   }
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
    delete[] ID;
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

    Proc.refine( 0 );

    Proc.convertStl2Morton( GMT.geom_nn, GMT.geom_xyz );

    for ( uint j = 0; j < proclevel - 1; j++ )
    {
        Proc.pushToRefinelist( j + 1 );
        Proc.fourToOne();
        Proc.refineRefineList();
        //       cout<<j<<" : "<<Proc.size()<<endl;
    }

    double t2;
    MPI_Barrier( MPI_COMM_WORLD );
    t2 = MPI_Wtime();
}

//==================================================================
//
//
//==================================================================
/*
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::derefineProcTopology()
{
// **************************************************************
//
//                 generate processor topology
//
// **************************************************************
   
    double t1;
    MPI_Barrier( MPI_COMM_WORLD );
    t1 = MPI_Wtime();

//    Proc.refine( 0 );

//    Proc.convertStl2Morton( GMT.geom_nn, GMT.geom_xyz );

    for ( uint j = 0; j < proclevel - 1; j++ )
    {
        Proc.pushToDerefinelist( j + 1 );
        Proc.retainFourToOne();
        Proc.derefineDerefineList();
        //       cout<<j<<" : "<<Proc.size()<<endl;
    }

    double t2;
    MPI_Barrier( MPI_COMM_WORLD );
    t2 = MPI_Wtime();
}
*/

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

        for ( uint j = 0; j < GMT.geom_nn; j++ )
        {
            if ( GMT.geom_xyz[3 * j + 0] >= xyz[0] && GMT.geom_xyz[3 * j + 0] <= xyz[1] )
            {
                if ( GMT.geom_xyz[3 * j + 1] >= xyz[2] && GMT.geom_xyz[3 * j + 1] <= xyz[3] )
                {
                    if ( GMT.geom_xyz[3 * j + 2] >= xyz[4] && GMT.geom_xyz[3 * j + 2] <= xyz[5] )
                    {
                        count++;
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
void ReblAmr<N, Nvalue, M, Mvalue>::moveGeometry( double *xx )
{
    Forest.moveGeom( Proc, GMT.geom_xyz, GMT.geom_nn, xx );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::getTotalMeshSize()
{
    Forest.getTotalMeshSize();
}
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::refineForest()
{
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
void ReblAmr<N, Nvalue, M, Mvalue>::asignWeightsForPart()
{
        real xyzc[3];
        uint co = 0;

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
        countPointsinBox();
        asignWeightsForPart();
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

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::allocatePointers(  )
{
 Forest.allocatePointersforEachBlock(Proc);
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::tagElemGradBased(  )
{

  Forest.addToRefineLisGradBased();

}


template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::solvePoisson(const real &Epsilon, const char*bc, const real* Dirichlet, const real* Neumann)      
{
Forest.solvePoisson(Proc,Epsilon,bc,Dirichlet,Neumann);

}


template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::orderAnalysis(const real &Epsilon, const char*bc, const real* Dirichlet, const real* Neumann)      
{
Forest.orderAnalysis(Proc,Epsilon,bc,Dirichlet,Neumann);

}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmr<N, Nvalue, M, Mvalue>::transInterpTest()
{
// this is a test function to verify transverse interpolations
	//Forest.preSwaptestTransInterpolation(Proc);

	Forest.postSwaptestTransInterpolation(Proc);

}





#endif





//template class ReblAmr<TREESIZE, real, PROCSIZE, uint>;
template class ReblAmr<TREESIZE, Q, PROCSIZE, uint>;
/*
template class ReblAmr<TREESIZE, real, PROCSIZE, uint, Tree<PROCSIZE, uint>>;
template class ReblAmr<TREESIZE, real, WSIZE, uint, FullTree<WSIZE, uint>>;
*/
