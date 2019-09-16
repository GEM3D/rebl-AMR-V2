#include "scale.h"

template <size_t N>
FullOctreeTop<N>::FullOctreeTop()
{
    if ( N % 3 != 0 )
    {
        throw std::runtime_error( "please adjust your full tree processortopology to a number divisible by 3\n" );
    }
    else
    {
        fullOctreeLevel = N / 3;
    }
}

template <size_t N>
void FullOctreeTop<N>::convertRank2Bits( int myrank )
{
    morton<N> key( myrank );
    rootKey = key;
    //  cout<<"key "<<key<<endl;
}

template <size_t N>
void FullOctreeTop<N>::constructNbrProcs()
{
    if ( N % 3 != 0 )
    {
        throw std::runtime_error( "please adjust your full tree processortopology to a number divisible by 3\n" );
    }

    // push back first the siblings
    // this is a full tree all we need to do is to find neighbors and find nonlocal neighbors
    morton<N> ktemp;

    ktemp = rootKey;

    // in x, y and z-direction
    ktemp.flip( 0 );
    Nbr.push_back( ktemp.to_ulong() );
    // cout<<ktemp<<endl;
    ktemp.flip( 0 );
    ktemp.flip( 1 );
    Nbr.push_back( ktemp.to_ulong() );

    // cout<<ktemp<<endl;
    ktemp.flip( 1 );
    ktemp.flip( 2 );
    Nbr.push_back( ktemp.to_ulong() );

    // cout<<ktemp<<endl;
    //

    // push back nonlocal neighbors

    for ( uint i = 0; i < 3; i++ )
    {
        if ( isBoundary( i ) )
        {
            continue;
        }
        else
        {
            ktemp = rootKey;
            ktemp.flip( N - 1 - 3 * ( fullOctreeLevel - 1 ) - i );
            ktemp.flip( N - 1 - 3 * ( fullOctreeLevel - 2 ) - i );
            Nbr.push_back( ktemp.to_ulong() );
        }
    }
}

// returns false if element is not a boundary element and true if it is a boundary element at a given direction

template <size_t N>
bool FullOctreeTop<N>::isBoundary( uint &direction )
{
    uint mylevel;
    bool bol1, bol2;
    morton<N> kt1, kt2;
    bol2 = true;

    uint j = direction;

    bol1 = rootKey[N - 1 - j];

    kt1 = 0;

    for ( uint k = 0; k < fullOctreeLevel; k++ )
    {
        if ( rootKey[N - 1 - 3 * k - j] != bol1 )
        {
            bol2 = false;
            break;
        }
    }

    return ( bol2 );
}

template <size_t N>
void FullOctreeTop<N>::checkGraphConsistency( int myrank )
/*!< check symmetry of the communication graph, if not symmetric, due to the use of MPI_Probe it will deadlock */
{
    int size;
    MPI_Status status, status1;
    MPI_Request request[Nbr.size()], request1;
    uint dest[Nbr.size()];

    for ( uint i = 0; i < Nbr.size(); i++ )
    {
        dest[i] = Nbr.at( i );
    }

    for ( uint i = 0; i < Nbr.size(); i++ )
    {
        MPI_Isend( dest, Nbr.size(), MPI_UNSIGNED, Nbr.at( i ), myrank, MPI_COMM_WORLD, &request1 );
    }

    uint *recvbuff[Nbr.size()];

    uint val;

    for ( uint i = 0; i < Nbr.size(); i++ )
    {
    //    cout << RED << "myrank " << myrank << " " << Nbr.at( i ) << RESET << endl;
        MPI_Probe( Nbr.at( i ), Nbr.at( i ), MPI_COMM_WORLD, &status );

        MPI_Get_count( &status, MPI_UNSIGNED, &size );

        recvbuff[i] = new uint[size];

        MPI_Irecv( recvbuff[i], size, MPI_UNSIGNED, Nbr.at( i ), Nbr.at( i ), MPI_COMM_WORLD, &request[i] );

        MPI_Wait( &request[i], &status1 );

        if ( std::find( recvbuff[i], recvbuff[i] + size, myrank ) == recvbuff[i] + size )
        {
            cout << "inconsistency occured at rank " << myrank << endl;
        }
    }

    for ( uint i = 0; i < Nbr.size(); i++ )
    {
        delete[] recvbuff[i];
    }
}

/*****************************************/

template <size_t N>
void FullOctreeTop<N>::readRoot( morton<N> &key )
{

    key = rootKey;
}

template <size_t N>
void FullOctreeTop<N>::readNbrs( vector<int> &Nbrs )
{

    for ( uint i = 0; i < Nbr.size(); i++ )
    {
        Nbrs.push_back( Nbr.at( i ) );
    }
}

/*

template <size_t N>
void   FullOctreeTop<N>::translateAndScaleGeom(int geom_nn, real* geom_xyz )
{



}
*/

template class FullOctreeTop<3>;
template class FullOctreeTop<15>;
template class FullOctreeTop<6>;
template class FullOctreeTop<9>;
