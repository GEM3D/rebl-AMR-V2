#include "tree.h"
#include "definitions.h"
#include "typedefs.h"

template <size_t N, typename value>
void FullTree<N, value>::setLevel( const uint &l )
{

    fixedlevel = l;
    // cout<<"fixedlevel"<<l<<endl;
}

template <size_t N, typename value>
void FullTree<N, value>::level( morton<N> key, uint *level )
{

    *level = fixedlevel;
}

template <size_t N, typename value>
typename bitmap<N, value>::iterator FullTree<N, value>::begin()
{
    return ( mesh.begin() );
}

template <size_t N, typename value>
typename bitmap<N, value>::iterator FullTree<N, value>::end()
{
    return ( mesh.end() );
}

template <size_t N, typename value>
void FullTree<N, value>::insertKey( morton<N> key )
{
    mesh.erase( 0 );
    mesh.insert( {key, nullptr} );
}

template <size_t N, typename value>
uint FullTree<N, value>::getLevel()
{

    return ( fixedlevel );
}

template <size_t N, typename value>
typename bitmap<N, value>::iterator FullTree<N, value>::find( morton<N> key ) /**!< this function is to find a value given the key*/
{
    typename bitmap<N, value>::iterator temp;

    temp = mesh.find( key );
    /*
    if(temp==end())
    {
    cout<<"key not found"<<endl;
    exit(0);
    }
    else
    {
    return(temp);
    }
    */
    return ( temp );
}

template <size_t N, typename value>
void FullTree<N, value>::nbrsConstrcut( vector<uint> &Nbrs, uint myrank )
{
    const uint size = WSIZE;

    Nbrs.clear();

    morton<size> ktemp;

    //    uint fullTreeLevel = getLevel();
    //  add siblings as neighbors

    morton<size> rootKey( myrank );
    ktemp = rootKey;
    // in x, y and z-direction
    ktemp.flip( size - 1 - 3 * ( fixedlevel - 1 ) );
    Nbrs.push_back( ktemp.to_ulong() );
    // cout<<ktemp<<endl;
    ktemp.flip( size - 1 - 3 * ( fixedlevel - 1 ) );
    ktemp.flip( size - 1 - 3 * ( fixedlevel - 1 ) - 1 );
    Nbrs.push_back( ktemp.to_ulong() );

    // cout<<ktemp<<endl;
    ktemp.flip( size - 1 - 3 * ( fixedlevel - 1 ) - 1 );
    ktemp.flip( size - 1 - 3 * ( fixedlevel - 1 ) - 2 );
    Nbrs.push_back( ktemp.to_ulong() );

    // cout<<ktemp<<endl;
    //

    // push back nonlocal neighbors
    // 3 stands for directions in 3D

    uint changedirectionlevel = 0;

    for ( uint i = 0; i < 3; i++ )
    {
        if ( isBoundary( i, myrank ) == 0 && fixedlevel > 1 )
        {

            ktemp = rootKey;
            findFlipLevel( ktemp, fixedlevel, &changedirectionlevel, &i );

            //           ktemp.flip( size - 1 - 3 * ( fixedlevel - 1 ) - i );
            //           ktemp.flip( size - 1 - 3 * ( fixedlevel - 2 ) - i );

            flipForNbr( &ktemp, fixedlevel, &changedirectionlevel, &i );

            Nbrs.push_back( ktemp.to_ulong() );
            //  cout<<" inside "<<ktemp<<endl;
        }
    }

    morton<N> ktmp2;

    for ( uint i = 0; i < Nbrs.size(); i++ )
    {
        morton<size> ktmp( Nbrs.at( i ) );
        ktmp2 = 0;

        for ( uint j = 0; j < 3 * fixedlevel; j++ )
        {
            ktmp2[N - j - 1] = ktmp[size - 1 - j];
        }
        mesh.insert( {ktmp2, nullptr} );
        //        cout << GREEN << Nbrs.at( i ) << " " << ktmp2 << RESET << endl;
    }

    //  cout << BLUE<<myrank<<" " << Nbrs.size() << RESET << endl;
    /*
        for ( uint i = 0; i < Nbrs.size(); i++ )
        {
            cout<<myrank<<" "<<Nbrs.at(i)<<endl;
        }
    */
}

template <size_t N, typename value>
void FullTree<N, value>::findFlipLevel( morton<WSIZE> key, uint fixedlevel, uint *changedirectionlevel, uint *direction )
{
    bool bol;

    bol = key[WSIZE - 3 * ( ( fixedlevel ) - 1 ) - ( *direction ) - 1];

    // cout<<"index is = "<<bit-3*((*mylevel)-1)-(*direction)-1<<endl;
    // cout<<"bol = "<<bol<<endl;

    *changedirectionlevel = 0;

    for ( uint i = ( fixedlevel ) - 1; i > 0; i-- )
    {
        // cout<<"i="<<i <<key[bit-3*(i-1)-(*direction)-1]<<endl;

        if ( key[WSIZE - 3 * ( i - 1 ) - ( *direction ) - 1] != bol )
        {
            *changedirectionlevel = i;
            break;
        }
    }
}

template <size_t N, typename value>
void FullTree<N, value>::flipForNbr( morton<WSIZE> *key, uint fixedlevel, uint *changedirectionlevel, uint *direction )
{
    // cout<<"mylevel, changelevel and direction = "<<*mylevel<<"\t"<<*changedirectionlevel<<"\t"<<*direction <<endl;
    // no change in level defines a boundary element assumping that the levek is bigger than 1, just return if this is the case
    if ( *changedirectionlevel == 0 )
    {
        return;
    }

    for ( uint i = ( *changedirectionlevel ); i <= ( fixedlevel ); i++ )
    {
        // cout<<(bit-3*(i-1)-(*direction)-1)<<endl;
        ( *key ).flip( WSIZE - 3 * ( i - 1 ) - ( *direction ) - 1 );
        // cout<<"in here"<<(*key)<<endl;
    }

    // cout<<"exiting flip NBR ..."<<endl;
}

template <size_t N, typename value>
bool FullTree<N, value>::isBoundary( uint &direction, uint myrank )
{
    uint mylevel;
    bool bol1, bol2;
    const uint size = WSIZE;

    morton<size> kt1, kt2, rootKey;
    bol2 = true;
    morton<size> key( myrank );
    rootKey = key;
    uint j = direction;
    bol1 = rootKey[size - 1 - j];
    kt1 = 0;

    for ( uint k = 0; k < fixedlevel; k++ )
    {
        if ( rootKey[size - 1 - 3 * k - j] != bol1 )
        {
            bol2 = false;
            break;
        }
    }

    return ( bol2 );
}

template <size_t N, typename value>
void FullTree<N, value>::assignProcs( vector<uint> &Nbrs, uint myrank )
{

    for ( uint i = 0; i < Nbrs.size(); i++ )
    {
        cout << RED << myrank << " " << Nbrs.at( i ) << RESET << endl;
    }

    const uint size = WSIZE;

    morton<size> ktmp;

    // hash table, order is modified

    for ( auto it = begin(); it != end(); it++ )
    {
        {
            it->second = new value[1];

            for ( uint j = 0; j < 3 * fixedlevel; j++ )
            {
                ktmp[size - 1 - j] = it->first[N - j - 1];
            }

            it->second[0] = ( value )( ktmp.to_ulong() );
            //          cout << RED << myrank << " " << Nbrs.at( co ) << it->first << " " << it->second[0] << RESET << endl;
        }
    }

    int co = 0;
    for ( auto it = mesh.begin(); it != mesh.end(); it++ )
    {

        cout << RED << myrank << " " << it->first << " " << it->second[0] << RESET << endl;
        co++;
    }

    // cout<<RED<<"co ="<<co<<RESET<<endl;
}

/*
template <size_t N, typename value>
void FullTree<N>::convertRank2Bits( int myrank )
{
    morton<3*fixedlevel> key( myrank );

   cout<<key<<endl;

}
*/

#if ( 0 )
template <size_t N, typename value>
void FullTree<N, value>::centroid( morton<N> key, real *xyz )
{
    real result = 0.0, sign = 0.0, dx = 0.0;
    // uint k;

    uint mylevel = fixedlevel;
    ;
    //    level( key, &mylevel );

    // cout<<"key="<<key<<endl;
    // cout<<"level"<<level<<endl;
    for ( uint j = 0; j < 3; j++ )
    {
        dx = 0.0;
        for ( uint i = 0; i < mylevel; i++ )
        {
            if ( key[N - 3 * i - j - 1] == false )
            {
                sign = -1.0;
            }
            else
            {
                sign = 1.0;
            }

            i++;
            TwoPowN( i, &result );
            // cout<<"result:"<<result<<endl;
            dx = dx + sign * 1. / result;
            i--;
        }
        xyz[j] = Tree<N, value>::ancestorcoords[j] + dx * ( Tree<N, value>::ancestorlength[j] ) * 0.5;
    }
}

template <size_t N, typename value>
void FullTree<N, value>::enclosingBox( morton<N> key, real *X )
{
    real result = 0.0, sign = 0.0, dx = 0.0;
    real xyz[3];
    real dy, dz;

    uint mylevel = fixedlevel;

    //    level( key, &mylevel );

    // cout<<key<<" mylevel is "<<mylevel<<endl;
    // get coords

    for ( uint j = 0; j < 3; j++ )
    {
        dx = 0.0;
        for ( uint i = 0; i < mylevel; i++ )
        {
            if ( key[N - 3 * i - j - 1] == false )
            {
                sign = -1.0;
            }
            else
            {
                sign = 1.0;
            }

            i++;
            TwoPowN( i, &result );
            // cout<<"result:"<<result<<endl;
            dx = dx + sign * 1. / result;
            i--;
        }
        xyz[j] = Tree<N, value>::ancestorcoords[j] + dx * ( Tree<N, value>::ancestorlength[j] ) * 0.5;
    }
    real idenum;
    TwoPowN( mylevel, &idenum );

    // cout<<"idenum"<<idenum<<endl;
    idenum = 1. / idenum;
    dx = Tree<N, value>::ancestorlength[0] * idenum * 0.5;
    X[0] = xyz[0] - dx;
    X[1] = xyz[0] + dx;

    // denum=pow(2,level);
    dy = Tree<N, value>::ancestorlength[1] * idenum * 0.5;
    X[2] = xyz[1] - dy;
    X[3] = xyz[1] + dy;

    // denum=pow(2,level);
    dz = Tree<N, value>::ancestorlength[2] * idenum * 0.5;
    X[4] = xyz[2] - dz;
    X[5] = xyz[2] + dz;
}

#endif

template <size_t N, typename value>
uint FullTree<N, value>::size()
{
    return ( mesh.size() );
}

template <size_t N, typename value>
void FullTree<N, value>::convertCoordToMorton( real *xyz, morton<N> &key ) /**!< this function is to find a value given the key*/
{
    real x, y, z, xmid, ymid, zmid;

    static const real half = 0.5;

    real xmin = Tree<N, value>::ancestorcoords[0] - half * Tree<N, value>::ancestorlength[0];
    real xmax = Tree<N, value>::ancestorcoords[0] + half * Tree<N, value>::ancestorlength[0];

    real ymin = Tree<N, value>::ancestorcoords[1] - half * Tree<N, value>::ancestorlength[1];
    real ymax = Tree<N, value>::ancestorcoords[1] + half * Tree<N, value>::ancestorlength[1];

    real zmin = Tree<N, value>::ancestorcoords[2] - half * Tree<N, value>::ancestorlength[2];
    real zmax = Tree<N, value>::ancestorcoords[2] + half * Tree<N, value>::ancestorlength[2];

    x = xyz[0];
    y = xyz[1];
    z = xyz[2];
    key = 0;

    for ( uint j = 0; j < N / 3; j++ )
    {
        xmid = half * ( xmin + xmax );

        ymid = half * ( ymin + ymax );

        zmid = half * ( zmin + zmax );

        if ( x > half * ( xmin + xmax ) )
        {
            key.flip( N - 1 - 3 * j );
            xmin = xmid;
        }
        else
        {
            xmax = xmid;
        }

        if ( y > half * ( ymin + ymax ) )
        {
            key.flip( N - 1 - 3 * j - 1 );
            ymin = ymid;
        }
        else
        {
            ymax = ymid;
        }
        if ( z > half * ( zmin + zmax ) )
        {
            key.flip( N - 1 - 3 * j - 2 );
            zmin = zmid;
        }
        else
        {
            zmax = zmid;
        }
    }

    mortonSTL.push_back( key );

    //cout << "=================== " << mortonSTL.size() << endl;
}

// template class FullTree<3, real>;
// template class FullTree<3, uint>;
// template class FullTree<6, uint>;
// template class FullTree<9, uint>;
template class FullTree<PROCSIZE, uint>;
template class FullTree<WSIZE, uint>;
template class FullTree<WSIZE, real>;
template class FullTree<32, real>;
//template class FullTree<9, real>;
