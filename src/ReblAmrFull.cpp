#include "ReblAmrFull.h"
#include "definitions.h"
#include <type_traits>

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
ReblAmrFull<N, Nvalue, M, Mvalue>::ReblAmrFull( int argcs, char *pArgs[], real *length, real *coords, uint nx, uint ny, uint nz )
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

/*
    GMT.construct( xyz1 );

    GMT.readSTLGeom( pArgs, xyz1 );


    Proc.construct( length, coords, nx, ny, nz );

    generateProcTopology();
*/
    //  FullProc.construct(length,coords, nx, ny, nz);

    // uses distributed proc topology, without success in zoltan part this will fail
    //
}


template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmrFull<N, Nvalue, M, Mvalue>::MPIStartUp()
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
void ReblAmrFull<N, Nvalue, M, Mvalue>::generateProcTopology()
{
    const uint proclevel = WSIZE / 3;
    Proc.setLevel( proclevel );
    morton<proclevel * 3> seedKey( Com.myrank );
    cout << seedKey << endl;
    morton<WSIZE> seed;
    for (int i = 0; i < 3 * proclevel; i++ )
    {
       // seed[TREESIZE - 1 - i] = seedKey[3 * proclevel - i - 1];
        seed[WSIZE - 1 - i] = seedKey[3 * proclevel - i - 1];
    }
    cout <<" seed "<< seed << endl;
    Proc.insertKey( seed );
    /*
       auto it=proc.begin();
       it->second=new uint[1];
       it->second[0]=(uint)Com.myrank;
       cout<<GREEN<<Com.myrank<<" "<<it->first<<" "<<it->second[0]<<RESET<<endl;
    */
    Proc.nbrsConstrcut( Nbrs, Com.myrank );
    //  std::sort(Nbrs.begin(),Nbrs.end());
    Proc.assignProcs( Nbrs, Com.myrank );
   
    Proc.convertStl2Morton( GMT.geom_nn, GMT.geom_xyz );


}


template <size_t N, typename Nvalue, size_t M, typename Mvalue>
ReblAmrFull<N, Nvalue, M, Mvalue>::~ReblAmrFull()
{

}


#if(1)
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmrFull<N, Nvalue, M, Mvalue>::countPointsinBox()
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
        cout<<" count = "<<count<<endl;
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
void ReblAmrFull<N, Nvalue, M, Mvalue>::forestConstruct( int argcs, char *pArgs[], real *length, real *coords, uint nx, uint ny, uint nz )
{
    Forest.construct( argcs, pArgs, Proc, length, coords, nx, ny, nz );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmrFull<N, Nvalue, M, Mvalue>::moveGeometry( double *xx )
{
    Forest.moveGeom( Proc, GMT.geom_xyz, GMT.geom_nn, xx );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmrFull<N, Nvalue, M, Mvalue>::getTotalMeshSize()
{
    Forest.getTotalMeshSize();
}
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmrFull<N, Nvalue, M, Mvalue>::refineForest()
{
    Forest.refineForestBalanced( meshlevel, Proc );
}
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmrFull<N, Nvalue, M, Mvalue>::setForestParams()
{
    Forest.getMaxSeedsLevel( Proc );
    Forest.setMaxProcLevel( proclevel );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmrFull<N, Nvalue, M, Mvalue>::createComPattern()
{
    setForestParams();
    Forest.comPatternConstruct( Proc, Nbrs );
}
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmrFull<N, Nvalue, M, Mvalue>::writeRunInfo()
{
    Forest.runInfo();
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void ReblAmrFull<N, Nvalue, M, Mvalue>::writeMesh( int index )
{
    //if ( std::is_same<P, Tree<N, uint>>::value )
    {
        templatePhdf5<N, Nvalue, M, Mvalue, FullTree<M, Mvalue>> IO;

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

//template class ReblAmrFull<TREESIZE, real, WSIZE, uint>;
/*
template class ReblAmrFull<TREESIZE, real, PROCSIZE, uint, Tree<PROCSIZE, uint>>;
template class ReblAmrFull<TREESIZE, real, WSIZE, uint, FullTree<WSIZE, uint>>;
*/
