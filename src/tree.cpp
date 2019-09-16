#include "tree.h"
#include "definitions.h"
#include "typedefs.h"

//====================================================================================
//
//                  Constructor
//
//====================================================================================

template <size_t N, typename value>
Tree<N, value>::Tree( real *length, real *coords, uint nx, uint ny, uint nz )
{
    for ( uint i = 0; i < 3; i++ )
    {
        ancestorlength[i] = length[i];
        ancestorcoords[i] = coords[i];
    }
    npx = nx;
    npy = ny;
    npz = nz;

    mesh.insert( {ancestorkey, nullptr} );

    // cout<<ancestorkey<<endl;
}

template <size_t N, typename value>
Tree<N, value>::Tree( real *length, real *coords )
{

    for ( uint i = 0; i < 3; i++ )
    {
        ancestorlength[i] = length[i];
        ancestorcoords[i] = coords[i];
    }
    npx = 2;
    npy = 2;
    npz = 2;

    mesh.insert( {ancestorkey, nullptr} );

    // cout<<ancestorkey<<endl;
}

template <size_t N, typename value>
Tree<N, value>::~Tree()
{
    for ( auto it = begin(); it != end(); it++ )
    {
        if ( it->second != nullptr )
        {
            delete[] it->second;
        }
    }
}

//====================================================================================
//
//                  Calculates the Level of the ELement from Morton Code
//
//====================================================================================

template <size_t N, typename value>
void Tree<N, value>::level( morton<N> key, uint *level )
{
    *level = N / 3;
    uint rem = N % 3;
    // cout<<"rem="<<rem<<endl; /*!< This remaining will affect the number of bits left over, this offset will be used in the following
    // loop*/
    uint iend = N / 3;
    morton<N> keytemp;

    // cout<<"rem "<<rem<<" iend "<<iend<<endl;

    /*! to prevent unnecesary bit operation, the morton code is placed from starting from  left hand side*/

    for ( uint i = 0; i < iend; i++ )
    {
        if ( key[3 * i + rem] == false && key[3 * i + rem + 1] == false && key[3 * i + rem + 2] == false )
        {
            /*! now look and see if any siblings exist*/
            //
            // cout<<"i=:"<<i <<endl;
            // cout<<key[3*i+rem]<<" "<<key[3*i+rem]<<" "<<key[3*i+rem]<<endl;
            keytemp = key;
            //  cout<<key<<endl;

            if ( mesh.count( keytemp.flip( 3 * i + rem ) ) == 0 )
            {
                // cout<<" i= "<<i<<"flipped key="<<keytemp<<endl;

                *level = *level - 1;

                // cout<<"here"<<endl;
            }
        }
        else
        {
            break;
        }
    }
}

//=====================================================================
//
//
// This routine calculated coordinates based on morton code
//  Tis is bsically a series of the form (-1)^boolean*2^-i
//
//=====================================================================

template <size_t N, typename value>
void Tree<N, value>::centroid( morton<N> key, real *xyz )
{
    real result = 0.0, sign = 0.0, dx = 0.0;
    // uint k;

    uint mylevel;
    level( key, &mylevel );

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
        xyz[j] = ancestorcoords[j] + dx * ( ancestorlength[j] ) * 0.5;
    }
}

template <size_t N, typename value>
void Tree<N, value>::centroidFixedLevel( morton<N> key, const uint mylevel, real *xyz )
{
    real result = 0.0, sign = 0.0, dx = 0.0;
    // uint k;

    //    uint mylevel;
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
        xyz[j] = ancestorcoords[j] + dx * ( ancestorlength[j] ) * 0.5;
    }
}

//====================================================================
//
//                             Enclosing Box
//
//====================================================================
template <size_t N, typename value>
void Tree<N, value>::enclosingBox( morton<N> key, real *X )
{
    real result = 0.0, sign = 0.0, dx = 0.0;
    real xyz[3];
    real dy, dz;

    uint mylevel;

    level( key, &mylevel );

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
        xyz[j] = ancestorcoords[j] + dx * ( ancestorlength[j] ) * 0.5;
    }

    real idenum;
    TwoPowN( mylevel, &idenum );

    // cout<<"idenum"<<idenum<<endl;
    idenum = 1. / idenum;
    dx = ancestorlength[0] * idenum * 0.5;
    X[0] = xyz[0] - dx;
    X[1] = xyz[0] + dx;

    // denum=pow(2,level);
    dy = ancestorlength[1] * idenum * 0.5;
    X[2] = xyz[1] - dy;
    X[3] = xyz[1] + dy;

    // denum=pow(2,level);
    dz = ancestorlength[2] * idenum * 0.5;
    X[4] = xyz[2] - dz;
    X[5] = xyz[2] + dz;
}

// enclosing box with given fixed level

template <size_t N, typename value>
void Tree<N, value>::enclosingBoxFixedLevel( morton<N> key, uint mylevel, real *X )
{
    real result = 0.0, sign = 0.0, dx = 0.0;
    real xyz[3];
    real dy, dz;

    //    uint mylevel;
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
        xyz[j] = ancestorcoords[j] + dx * ( ancestorlength[j] ) * 0.5;
    }

    real idenum;
    TwoPowN( mylevel, &idenum );

    // cout<<"idenum"<<idenum<<endl;
    idenum = 1. / idenum;
    dx = ancestorlength[0] * idenum * 0.5;
    X[0] = xyz[0] - dx;
    X[1] = xyz[0] + dx;

    // denum=pow(2,level);
    dy = ancestorlength[1] * idenum * 0.5;
    X[2] = xyz[1] - dy;
    X[3] = xyz[1] + dy;

    // denum=pow(2,level);
    dz = ancestorlength[2] * idenum * 0.5;
    X[4] = xyz[2] - dz;
    X[5] = xyz[2] + dz;
}

//======================================================================
//
//                  these two functions return theiterator object for
//                      begining and end of the unordered_map
//
//======================================================================

template <size_t N, typename value>
typename bitmap<N, value>::iterator Tree<N, value>::begin()
{
    return ( mesh.begin() );
}

template <size_t N, typename value>
typename bitmap<N, value>::iterator Tree<N, value>::end()
{
    return ( mesh.end() );
}
//======================================================================
//
//    return the Size of mesh and reserve the amount you want,
//          reserve is important to prevent rehashing
//
//======================================================================

template <size_t N, typename value>
uint Tree<N, value>::size()
{
    return ( mesh.size() );
}

/* \brief this initial guess is important to prevent from rehashing*/
template <size_t N, typename value>
void Tree<N, value>::reserve( uint *reservedsize )
{
    mesh.reserve( *reservedsize );
}

//=========================================================
//
//                       Get Siblings
//
//==========================================================

template <size_t N, typename value>
void Tree<N, value>::siblings( morton<N> key, uint mylevel, morton<N> *sibkey )
{
    morton<N> kt = key;

    mylevel = mylevel - 1;

    sibkey[0] = kt.flip( N - 3 * ( mylevel ) - 1 );
    // cout<<"sibkeyi[0] "<<sibkey[0]<<endl;

    kt = key;
    sibkey[1] = kt.flip( N - 3 * ( mylevel ) - 2 );
    // cout<<"sibkeyi[1] "<<sibkey[1]<<endl;

    kt = key;
    sibkey[2] = kt.flip( N - 3 * ( mylevel ) - 3 );
    // cout<<"sibkeyi[2] "<<sibkey[2]<<endl;

    kt = key;
    kt.flip( N - 3 * ( mylevel ) - 1 );
    sibkey[3] = kt.flip( N - 3 * ( mylevel ) - 2 );
    // cout<<"sibkeyi[3] "<<sibkey[3]<<endl;

    kt = key;
    kt.flip( N - 3 * ( mylevel ) - 1 );
    sibkey[4] = kt.flip( N - 3 * ( mylevel ) - 3 );
    // cout<<"sibkeyi[4] "<<sibkey[4]<<endl;

    kt = key;
    kt.flip( N - 3 * ( mylevel ) - 3 );
    sibkey[5] = kt.flip( N - 3 * ( mylevel ) - 2 );
    // cout<<"sibkeyi[3] "<<sibkey[5]<<endl;

    kt = key;
    kt.flip( N - 3 * ( mylevel ) - 1 );
    kt.flip( N - 3 * ( mylevel ) - 2 );

    sibkey[6] = kt.flip( N - 3 * ( mylevel ) - 3 );
}

//==================================================================
//
//                  Refines a single elmenet
//
//==================================================================

template <size_t N, typename value>
void Tree<N, value>::refine( morton<N> key )
{
    uint mylevel;
    level( key, &mylevel );
    morton<N> temp;
    temp = key;
    /* \brief if the morton code does not exist in mesh, refinement is not permitted (refining a nonexsiting element not permitted)*/
    /* \brief for now I have assigned the values as NULL, assign it to what you like for solving your particular problem*/

/*
    if ( mesh.count( key ) == 0 )
    {
        cout << key << endl;
        //           throw std::runtime_error( RED"Erro in refining key does not exist"RESET );
    }
*/
    // cout<<"N= "<<N<<endl;
    // cout<<"mylevel "<<mylevel<<endl;
    // cout<<"key "<<key<<endl;
    // cout<<"index "<<N-3*mylevel-1<<endl;

    // too many type conversions but the asnwer will always be positive so should not be a problem

    temp = key;
    temp.flip( N - 3 * mylevel - 1 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 2 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 3 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 1 );
    temp.flip( N - 3 * mylevel - 2 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 1 );
    temp.flip( N - 3 * mylevel - 3 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 2 );
    temp.flip( N - 3 * mylevel - 3 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 1 );
    temp.flip( N - 3 * mylevel - 2 );
    temp.flip( N - 3 * mylevel - 3 );
    mesh.insert( {temp, nullptr} );
}

//===================================================================
// \brief
//@return perfomes 4:1 consistent refinement for the elements
// in the refinelist, this list is empty after return
//
//
//===================================================================

template <size_t N, typename value>
void Tree<N, value>::refineRefineList()
{
    morton<N> key;
    for ( auto it = refinelist.begin(); it != refinelist.end(); it++ )
    {
        key = it->first;

        refine( key );
    }
    refinelist.clear();
}

template <size_t N, typename value>
void Tree<N, value>::refineRefineList( bitvector<N> &V )
{
    for ( uint i = 0; i < V.size(); i++ )
    {
        auto it = refinelist.find( V.at( i ) );
        refine( it->first );
        refinelist.erase( it->first );
    }
}

template <size_t N, typename value>
void Tree<N, value>::refineRefineList( uint istart, uint iend )
{
    morton<N> key;

    for ( uint i = istart; i < iend; i++ )
    {
        key = refinelist.at( i );
        refine( key );
        //       refinelist.pop_back();
    }
}
template <size_t N, typename value>
void Tree<N, value>::clearRefineList()
{
    refinelist.clear();
}

//====================================================================
//
//             derefine derefinelist
//
// ===================================================================
/*
template <size_t N, typename value>
uint Tree<N, value>::( const morton<N> key, const real *geom_xyz, uint n )
{





}
*/
//=============================================================
//
//
//            Check and see if that is inside solid
//
//============================================================
template <size_t N, typename value>
uint Tree<N, value>::isInsideSolid( const morton<N> key, const real *geom_xyz, uint n )
{
    uint j;
    uint a, b, c;
    real xyz[6];
    uint bol;

    bol = 0;
    // get bounding box

    enclosingBox( key, xyz );

    for ( j = 0; j < n; j++ )
    {
        a = geom_xyz[3 * j] > xyz[0] && geom_xyz[3 * j] < xyz[1];
        b = geom_xyz[3 * j + 1] > xyz[2] && geom_xyz[3 * j + 1] < xyz[3];
        c = geom_xyz[3 * j + 2] > xyz[4] && geom_xyz[3 * j + 2] < xyz[5];

        if ( a && b && c )
        {
            bol = 1;
            break;
        }
    }

    return ( bol );
}

// ================================================================
//
//                       Four to One balance
//
// ================================================================

template <size_t N, typename value>
uint Tree<N, value>::IsInVectorList( morton<N> key )
{
    uint bol = 0;

    for ( uint i = 0; i < refinelist.size(); i++ )
    {
        if ( key == refinelist.at( i ) )
        {
            bol = 1;
            break;
        }
    }
    return ( bol );
}

//====================================================================
//
/*!< all we are interested is the nonlocal neighbors, i.e. the neighbors of the parents as siblings will have same level*/

//====================================================================
#if ( 0 )
template <size_t N, typename value>
void Tree<N, value>::fourToOne() /*!< imposes 4:1 balance given the list of elments to be refined in the vector refine list*/
{
    morton<N> kt;
    uint mylevel, nbrlevel, a, changedirectionlevel;
    morton<N> mykey;
    uint istart = 0;
    uint iend = refinelist.size();

    a = 1;
    // the so-called ripple effect

    while ( a == 1 )
    {
        a = 0;

        //      printf("istart %d iend %d\n",istart,iend);

        for ( uint i = istart; i < iend; i++ )
        {
            mykey = refinelist.at( i );
            level( mykey, &mylevel );

            if ( mylevel > 1 )
            {
                for ( uint j = 0; j < 3; j++ )
                {
                    //      uint j=1;

                    kt = mykey;

                    findFlipLevel( kt, &mylevel, &changedirectionlevel, &j );
                    //      cout<<"mykey\t"<<kt<<endl;
                    // if the change in signe does not happen that is a boundary cube
                    //       cout<<"changedirectionlevel"<<changedirectionlevel<<endl;
                    if ( changedirectionlevel != 0 )
                    {
                        flipForNbr( &kt, &mylevel, &changedirectionlevel, &j );
                        mesh.count( kt );

                        // if this element exists, the level of nbr>=level of tagged element

                        //      cout<<"kt="<<mesh.count(kt)<<endl;
                        if ( mesh.count( kt ) == 0 )
                        {
                            kt[N - 3 * ( mylevel - 1 ) - 1] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 2] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 3] = 0;
                        }
                        //      cout<<"modified key\t"<<kt<<endl;
                        level( kt, &nbrlevel );

                        // cout<<"nbr level"<<nbrlevel<<"my level"<<mylevel<<endl;

                        // if its is not in the list add to list
                        //                   if ( mylevel > nbrlevel && !IsInVectorList( kt ) )
                        // use the new version which is std::find is much faster than my search
                        if ( mylevel > nbrlevel && std::find( refinelist.begin(), refinelist.end(), kt ) == refinelist.end() )

                        {
                            refinelist.push_back( kt );
                            a = 1;
                        }
                    }
                }
            }
        }

        if ( a == 1 )
        {
            istart = iend;
            iend = refinelist.size();
            //  cout<<"istart=\t"<<istart<<"iend\t"<<iend<<endl;
        }
    }

    /*!< this approach eliminates search algorithm  as now we do not have the restrictions on cutting the cube that we had in the previous
     * approach*/
    // std::sort (refine_list.begin(), refine_list.end(), compare_level);
    // cout<<"============================================================"<<endl;
    // cout<<           "EXISTING 4:1 BALANCE CHECK" <<endl;
    // cout<<"============================================================"<<endl;
}

//=================================================================================================
//
//    4:1 balance without a while loop; i.e. parallel 4:1
//    original older version
//================================================================================================

template <size_t N, typename value>
void Tree<N, value>::fourToOneP( uint istart,
                                 uint iend ) /*!< imposes 4:1 balance given the list of elments to be refined in the vector refine list*/
{
    morton<N> kt;
    uint mylevel, nbrlevel, a, changedirectionlevel;
    morton<N> mykey;
    //    uint      istart = 0;
    //    uint      iend   = refinelist.size();
    //      printf("istart %d iend %d\n",istart,iend);

    for ( uint i = istart; i < iend; i++ )
    {
        mykey = refinelist.at( i );
        level( mykey, &mylevel );

        if ( mylevel > 1 )
        {
            for ( uint j = 0; j < 3; j++ )
            {
                //      uint j=1;

                kt = mykey;

                //      cout<<"mykey\t"<<kt<<endl;
                findFlipLevel( kt, &mylevel, &changedirectionlevel, &j );

                // if the change in signe does not happen that is a boundary cube
                //       cout<<"changedirectionlevel"<<changedirectionlevel<<endl;
                if ( changedirectionlevel != 0 )
                {
                    flipForNbr( &kt, &mylevel, &changedirectionlevel, &j );
                    mesh.count( kt );

                    // if this element exists, the level of nbr>=level of tagged element

                    //      cout<<"kt="<<mesh.count(kt)<<endl;
                    if ( mesh.count( kt ) == 0 )
                    {
                        kt[N - 3 * ( mylevel - 1 ) - 1] = 0;
                        kt[N - 3 * ( mylevel - 1 ) - 2] = 0;
                        kt[N - 3 * ( mylevel - 1 ) - 3] = 0;
                    }
                    //      cout<<"modified key\t"<<kt<<endl;
                    level( kt, &nbrlevel );

                    // cout<<"nbr level"<<nbrlevel<<"my level"<<mylevel<<endl;

                    // if its is not in the list add to list
                    //   if ( mylevel > nbrlevel && !IsInVectorList( kt ) )
                    // use the new version which is std::find is much faster than my search
                    if ( mylevel > nbrlevel && std::find( refinelist.begin(), refinelist.end(), kt ) == refinelist.end() )

                    {
                        refinelist.push_back( kt );
                        //                    a = 1;
                    }
                }
            }
        }
    }
}

#endif
//==========================================================================================
// eliminated search
// =========================================================================================

template <size_t N, typename value>
void Tree<N, value>::fourToOne() /*!< imposes 4:1 balance given the list of elments to be refined in the vector refine list*/
{
    morton<N> kt;
    uint mylevel, nbrlevel, a, changedirectionlevel;
    morton<N> mykey;

    a = 1;
    // the so-called ripple effect

    while ( a == 1 )
    {
        a = 0;

        //      printf("istart %d iend %d\n",istart,iend);

        for ( auto it = refinelist.begin(); it != refinelist.end(); it++ )
        {
            mykey = it->first;
            level( mykey, &mylevel );

            if ( mylevel > 1 )
            {
                for ( uint j = 0; j < 3; j++ )
                {
                    //      uint j=1;

                    kt = mykey;

                    //      cout<<"mykey\t"<<kt<<endl;
                    findFlipLevel( kt, &mylevel, &changedirectionlevel, &j );

                    // if the change in signe does not happen that is a boundary cube
                    //       cout<<"changedirectionlevel"<<changedirectionlevel<<endl;
                    if ( changedirectionlevel != 0 )
                    {
                        flipForNbr( &kt, &mylevel, &changedirectionlevel, &j );
                        mesh.count( kt );

                        // if this element exists, the level of nbr>=level of tagged element

                        //      cout<<"kt="<<mesh.count(kt)<<endl;
                        if ( mesh.count( kt ) == 0 )
                        {
                            kt[N - 3 * ( mylevel - 1 ) - 1] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 2] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 3] = 0;
                        }

                        level( kt, &nbrlevel );

                        // cout<<"nbr level"<<nbrlevel<<"my level"<<mylevel<<endl;

                        // if its is not in the list add to list
                        // if ( mylevel > nbrlevel && !IsInVectorList( kt ) )
                        // use std find is faster

                        // edit

                        if ( mylevel > nbrlevel && it->second == 0 )
                        {
                            refinelist.insert( {kt, 0} );
                            //
                            a = 1;
                        }

                        /*
                              else
                                               {
                                                  refine(mykey);
                                                  refinelist1.erase(mykey);
                                                }
                                               */
                    }
                }
                it->second = 1;
                // refine(mykey);
                // refinelist1.erase(mykey);
            }
        }
    }

    // need to reset this for parallel applications

    /*!< this approach eliminates search algorithm  as now we do not have the restrictions on cutting the cube that we had in the previous
     * approach*/
    // std::sort (refine_list.begin(), refine_list.end(), compare_level);
    // cout<<"============================================================"<<endl;
    // cout<<           "EXISTING 4:1 BALANCE CHECK" <<endl;
    // cout<<"============================================================"<<endl;
}

//
//===============================================================================================
// Added this due to need to access this from within forest
//

template <size_t N, typename value>
uint Tree<N, value>::refineListSize()
{
    return ( refinelist.size() );
}

template <size_t N, typename value>
void Tree<N, value>::refinelistReset()
{

    if ( refinelist.size() != 0 )
    {
        for ( auto it = refinelist.begin(); it != refinelist.end(); it++ )

        {
            it->second = 0;
        }
    }
}

//=================================================================================================
//=================================================================================================
// gets the first change of direction
//
template <size_t N, typename value>
void Tree<N, value>::flipForNbr( morton<N> *key, uint *mylevel, uint *changedirectionlevel, uint *direction )
{
    // cout<<"mylevel, changelevel and direction = "<<*mylevel<<"\t"<<*changedirectionlevel<<"\t"<<*direction <<endl;
    // no change in level defines a boundary element assumping that the levek is bigger than 1, just return if this is the case
    if ( *changedirectionlevel == 0 )
    {
        return;
    }

    for ( uint i = ( *changedirectionlevel ); i <= ( *mylevel ); i++ )
    {
        // cout<<(bit-3*(i-1)-(*direction)-1)<<endl;
        ( *key ).flip( N - 3 * ( i - 1 ) - ( *direction ) - 1 );
        // cout<<"in here"<<(*key)<<endl;
    }

    // cout<<"exiting flip NBR ..."<<endl;
}

template <size_t N, typename value>
void Tree<N, value>::addToList( morton<N> key )
{
    refinelist.insert( {key, 0} );

    // cout<<"refinesize"<<refinelist.size()<<endl;
}

//
//   Given the Key find the level to be traversed to find the nbr
//   \brief
//   a) there is not change in sign in a given direction for the entire code, this corresponds to boundary element
//   b) otherwise, find the level at which the change in sign appears
//   c) siblings normally exist
//   direction gets the values 1, 2,3 for x,y, azd z directions
//=============================================================================================

template <size_t N, typename value>
void Tree<N, value>::findFlipLevel( morton<N> key, uint *mylevel, uint *changedirectionlevel, uint *direction )
{
    bool bol;

    bol = key[N - 3 * ( ( *mylevel ) - 1 ) - ( *direction ) - 1];

    // cout<<"index is = "<<bit-3*((*mylevel)-1)-(*direction)-1<<endl;
    // cout<<"bol = "<<bol<<endl;

    *changedirectionlevel = 0;

    for ( uint i = ( *mylevel ) - 1; i > 0; i-- )
    {
        // cout<<"i="<<i <<key[bit-3*(i-1)-(*direction)-1]<<endl;

        if ( key[N - 3 * ( i - 1 ) - ( *direction ) - 1] != bol )
        {
            *changedirectionlevel = i;
            break;
        }
    }

    // cout<<*changedirectionlevel<<endl;

    if ( *changedirectionlevel == 0 )
    {
        // cout<<"element is a boundary element = "<<endl;
    }

    // cout<<"exiting flip level ..."<<endl;
}
//=================================================
template <size_t N, typename value>
uint Tree<N, value>::count( morton<N> key )
{
    return ( mesh.count( key ) );
}

#if ( 0 )
template <size_t N, typename value>
void Tree<N, value>::derefineDerefineList()
{
    morton<N> key;
    bool bol;
    bitvector<N> temp;
    morton<N> sibkey[8];

    // std::sort(derefinelist.begin(),derefinelist.end(), compare);

    uint mylevel;
    uint maxlevel = 1;

    // get max level

    for ( uint i = 0; i < derefinelist.size(); i++ )
    {
        key = derefinelist.at( i );
        level( key, &mylevel );
        if ( mylevel > maxlevel )
        {
            maxlevel = mylevel;
        }
    }

    for ( uint j = maxlevel; j > 0; j-- )
    {
        for ( uint i = 0; i < derefinelist.size(); i++ )
        {
            key = derefinelist.at( i );
            level( key, &mylevel );

            if ( mylevel == j )
            {
                temp.push_back( key );
            }
        }
    }
    /*
    for(uint i=0;i<temp.size();i++)
    {
    cout<<"temp"<<temp.at(i)<<endl;
    }
    */

    for ( uint i = 0; i < temp.size(); i++ )
    {
        key = temp.at( i );
        level( key, &mylevel );
      //  cout << "temp level= " << mylevel << endl;
    }

    if ( temp.size() != derefinelist.size() )
    {
        cout << RED "missing element after sorting" RESET << endl;
    }

    uint changedirectionlevel, nbrlevel;

    morton<N> kt;

    for ( uint i = 0; i < temp.size(); i++ )
    // while(temp.size()!=0)
    {
        key = temp.at( i );
        //
        //
        // key=temp.front();
        // cout<<"tagged elem"<<key<<endl;

        bol = true;

        sibkey[7] = key;

        level( key, &mylevel );

        siblings( key, mylevel, sibkey );

        /*
        for(uint l=0;l<8;l++)
        {
        cout<<sibkey[l]<<endl;
        }
        */
        for ( uint j = 0; j < 8; j++ )
        {
            /*
            cout<<"==========================================\n"<<endl;
            cout<<"key is " <<key<<endl;
            cout<<"level is " <<mylevel<<endl;
            */
            for ( uint k = 0; k < 3; k++ )
            {
                key = sibkey[j];
                findFlipLevel( key, &mylevel, &changedirectionlevel, &k );
                // printf("k= %d\n",k);
                // cout<<"chanege dir "<<changedirectionlevel<<endl;
                if ( changedirectionlevel != 0 )
                {
                    flipForNbr( &key, &mylevel, &changedirectionlevel, &k );
                    mesh.count( key );
                    // cout<<"nbr key "<<key<<endl;
                    // if this element exists, the level of nbr>=level of tagged element

                    //      cout<<"kt="<<mesh.count(kt)<<endl;
                    if ( mesh.count( key ) == 0 )
                    {
                        key[N - 3 * ( mylevel - 1 ) - 1] = 0;
                        key[N - 3 * ( mylevel - 1 ) - 2] = 0;
                        key[N - 3 * ( mylevel - 1 ) - 3] = 0;
                    }
                    //      cout<<"modified key\t"<<kt<<endl;
                    level( key, &nbrlevel );

                    // cout<<"nbr level =  "<<nbrlevel<<" my level = "<<mylevel<<endl;

                    // if its is not in the list add to list
                    if ( mylevel < nbrlevel )
                    {
                        bol = false;
                        break;
                    }
                }
            }
            // cout<<"==========================================\n"<<endl;
        }

        // cout<<"bol= "<<bol<<endl;

        if ( bol == true )
        {
            key = sibkey[7];
            derefine( key );
        }

        // cout<<"*************************************************************"<<endl;
        // printKey();
        // cout<<"*************************************************************"<<endl;
        // cout<<RED "Mesh Size" RESET<<mesh.size()<<endl;
    }
    derefinelist.clear();
    derefinelist.shrink_to_fit();
}

#endif

template <size_t N, typename value>
typename bitmap<N, value>::iterator Tree<N, value>::find( morton<N> key ) /**!< this function is to find a value given the key*/
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
// typename bitmap<N, int>::iterator Tree<N, value>::findInList( morton<N> key ) /**!< this function is to find a value given the key*/
typename std::unordered_map<morton<N>, int>::iterator Tree
<N, value>::findInList( morton<N> key ) /**!< this function is to find a value given the key*/
{
    //    typename bitmap<N, int>::iterator temp;

    typename std::unordered_map<morton<N>, int>::iterator temp;

    temp = refinelist.find( key );
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

//
// This is not going tobe used in Forest, I have another loop for this
//

template <size_t N, typename value>
void Tree<N, value>::convertStl2Morton( uint geom_size, real *geom_xyz ) /**!< this function is to find a value given the key*/
{

//     mortonSTL.clear();

    real x, y, z, xmid, ymid, zmid;
    morton<N> key = 0;

    real xyz[3];

    // std::cout<<"xmax "<<xmax<<" ymax "<<ymax<<" zmax "<<zmax<<endl;

    // std::cout<<"xmin "<<xmin<<" ymin "<<ymin<<" zmin "<<zmin<<endl;
    std::cout << GREEN << "Morton Code Construction for STL Successful it's size is before accumulation " << mortonSTL.size() <<RESET << endl;

    for ( uint i = 0; i < geom_size; i++ )
    {
        xyz[0] = geom_xyz[3 * i];
        xyz[1] = geom_xyz[3 * i + 1];
        xyz[2] = geom_xyz[3 * i + 2];
        convertCoordToMorton( xyz, key );

// pushing back inside the above function
        //       mortonSTL.push_back( key );
    }

// sort this in the future
// std::sort(mortonSTL.begin(),mortonSTL.end());

    real X[6];

#if ( 0 )
    bool bol1, bol2, bol3;

    for ( uint i = 0; i < mortonSTL.size(); i++ )
    {
        x = geom_xyz[3 * i];
        y = geom_xyz[3 * i + 1];
        z = geom_xyz[3 * i + 2];

        key = mortonSTL.at( i );
        enclosingBox( key, X );

        bol1 = X[0] <= x && x <= X[1];
        bol2 = X[2] <= y && y <= X[3];
        bol3 = X[4] <= z && z <= X[5];

        // cout<<bol1<< " " <<bol2 << " "<<bol3 <<endl;
        //      cout<< x << " " <<y<< " " <<z <<endl;
        //    cout<<"  "<<X[0]<<" "<<X[1]<<" " <<X[2]<< " "<<X[3]<<" "<<X[4]<<" "<<X[5]<<endl;
        //      cout<< "key "<<key<<endl;
        if ( !( bol1 && bol2 && bol3 ) )
        {
       //     throw ::runtime_error( RED "error in generating morton code for STL geometry" RESET );
      //  std::cout << YELLOW << "MWarning : the key does not belong to any element in encoded Geometry STL " << mortonSTL.size() <<RESET << endl;
        
        }
    }


#endif
       std::cout << GREEN << "Morton Code Construction for STL Successful it's size is " << mortonSTL.size() <<RESET << endl;
    /*
    std::cout<<"x "<<x<<"y "<<y<<"z "<<z<<endl;
    std::cout<<key<<endl;
    std::cout<<" xmin "<<X[0]<<" xmax "<<X[1]<<endl;

    std::cout<<" ymin "<<X[2]<<" zmax "<<X[3]<<endl;
    std::cout<<" zmin "<<X[4]<<" zmax "<<X[5]<<endl;
    */
}

template <size_t N, typename value>
void Tree<N, value>::convertCoordToMorton( real *xyz, morton<N> &key ) /**!< this function is to find a value given the key*/
{
    real x, y, z, xmid, ymid, zmid;

    static const real half = 0.5;

    real xmin = ancestorcoords[0] - half * ancestorlength[0];
    real xmax = ancestorcoords[0] + half * ancestorlength[0];

    real ymin = ancestorcoords[1] - half * ancestorlength[1];
    real ymax = ancestorcoords[1] + half * ancestorlength[1];

    real zmin = ancestorcoords[2] - half * ancestorlength[2];
    real zmax = ancestorcoords[2] + half * ancestorlength[2];

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
}

template <size_t N, typename value>
void Tree<N, value>::pushToRefinelist( uint nlevel ) /**!< this function is to find a value given the key*/
{
    morton<N> key, key1;

    key = 0;
    //  cout<<BLUE"mortonsize "<<mortonSTL.size()<<RESET<<endl;
    uint mylevel = 0;

    // cout<<" nlevel "<<nlevel<<endl;

    for ( uint i = 0; i < mortonSTL.size(); i++ )
    {
        key = mortonSTL.at( i );

        for ( uint j = 0; j < 3 * nlevel; j++ )
        {
            key1[N - j - 1] = key[N - j - 1];
        }
        // cout<<key1<<" "<<key1<<endl;
        // considering dynamic mesh, we need to check if that element actually exists in the mesh
        level( key1, &mylevel );
        //      cout<<key1<<"mylevel "<<mylevel<<"nlevel "<<nlevel<<endl;
        if ( refinelist.count( key1 ) == 0 && mesh.count( key1 ) != 0 && mylevel == nlevel )
        {
            //      addToList( {key1,0} );
            refinelist.insert( {key1, 0} );
        }
    }
  //  cout << refinelist.size() << endl;
  /*  
    for(uint i=0;i<refinelist.size();i++)
    {
    std::cout<<refinelist.at(i)<<endl;
    }
   */
}

template <size_t N, typename value>
void Tree<N, value>::construct( real *length, real *coords, uint nx, uint ny, uint nz )
{
    for ( uint i = 0; i < 3; i++ )
    {
        ancestorlength[i] = length[i];
        ancestorcoords[i] = coords[i];
    }
    npx = nx;
    npy = ny;
    npz = nz;

    mesh.insert( {ancestorkey, nullptr} );
}

/*
template <size_t N, typename value>
void Tree<N, value>::setToZero(  )
{

    for ( auto it = mesh.begin(); it != mesh.end() ; it++ )
    {
        it->second=new value[1];        
        it->second[0]= ;        
    }

}
*/
/*
template <size_t N, typename value>
void Tree<N, value>::ignoreInactive(morton<N> key )
{
uint geomlevel;

level(key, &geomlevel );

vector<morton<N>> temp;

morton<N> tmp=0,mykey;

for(int i=0;i<3*geomlevel;i++)
{
   tmp.flip(N-i-1);
}

for(auto it = derefinelist.begin();it!= derefinelist.end();it++)
{
   mykey= it->first;
   mykey |= tmp;

   if(mykey==key)
   { 
   temp.push_back(mykey);
   }
}
    
while(temp.size()!=0)
{
  derefinelist.erase(temp.at(temp.size()-1));
}

}
*/

template <size_t N, typename value>
void Tree<N, value>::ignoreInactiveVertices(int nInactive, real * box)
{

real encBox[6];
bool bol0,bol1,bol2;

vector<morton<N>> temp;

cout<<" size of derefinelist before = " <<derefinelist.size()  <<endl;

for(int i=0;i<nInactive;i++)
{
for(auto it = derefinelist.begin();it!= derefinelist.end();it++)
{

   enclosingBox(it->first,encBox);

 bol0 =  (encBox[0] < box[6*i+1] && encBox[0] > box[6*i+0]) || (encBox[1] < box[6*i+1] && encBox[1] > box[6*i+0]);

 bol1 = (encBox[2] < box[6*i+3] && encBox[2] > box[6*i+2]) || (encBox[3] < box[6*i+3] && encBox[3] > box[6*i+2]);
   
 bol2 = (encBox[4] < box[6*i+5] && encBox[4] > box[6*i+4]) || (encBox[5] < box[6*i+5] && encBox[5] > box[6*i+4]);


 
/*
   bol0= (center[0] < box[6*i+1] && center[0]> box[6*i+0]); 
   bol1= (center[1] < box[6*i+3] && center[1]> box[6*i+2]); 
   bol2= (center[2] < box[6*i+5] && center[2]> box[6*i+4]); 
*/
 //cout<<" size of ignorance " <<temp.size()  <<endl;

  if(bol0 && bol1 && bol2 )                     
  {
   temp.push_back(it->first);
  }

 }
}

 //cout<<"Tree.cpp ????????????????? size of derifinelist after " <<temp.size()  <<endl;
    
while(temp.size()!=0)
{
  derefinelist.erase(temp.at(temp.size()-1));
  temp.pop_back();
}
cout<<" size of derefinelist after = " <<derefinelist.size()  <<endl;


}




template <size_t N, typename value>
void Tree<N, value>::ignoreInactive(int nInactive, real * box)
{

real center[6];
bool bol0,bol1,bol2;

vector<morton<N>> temp;

cout<<" size of derefinelist " <<derefinelist.size()  <<endl;

for(int i=0;i<nInactive;i++)
{
for(auto it = derefinelist.begin();it!= derefinelist.end();it++)
{

   centroid(it->first,center); 

   bol0= (center[0] < box[6*i+1] && center[0]> box[6*i+0]); 
   bol1= (center[1] < box[6*i+3] && center[1]> box[6*i+2]); 
   bol2= (center[2] < box[6*i+5] && center[2]> box[6*i+4]);
//  cout<<"center[0] = "<<center[0]<<", center[1] = " << center[1] << ", center[2] = " << center[2]<<endl; 

// cout<<" size of ignorance " <<temp.size()  <<endl;

  if(bol0 && bol1 && bol2 )                     
  {
   temp.push_back(it->first);

  cout<<"center[0] = "<<center[0]<<", center[1] = " << center[1] << ", center[2] = " << center[2]<<endl; 


  }

 }
}

 cout<<"Tree.cpp ????????????????? size of derifinelist after " <<temp.size()  <<endl;
    
while(temp.size()!=0)
{
  derefinelist.erase(temp.at(temp.size()-1));
  temp.pop_back();
}

}



template <size_t N, typename value>
void Tree<N, value>::insertKey( morton<N> key )
{
    mesh.insert( {key, nullptr} );
}

// you can eliminate level calculations and save 100100100... up to 64 bits for x-direction and 010010010... up to 64 for y-direction and
// ...
template <size_t N, typename value>
bool Tree<N, value>::isBoundary( morton<N> &key )
{
    uint mylevel;
    bool bol1, bol2;
    morton<N> kt1, kt2;
    level( key, &mylevel );
    bol2 = false;
    // cout<<"key= "<<key<<"level "<<mylevel<<endl;
    // need to chec for all three directions for every element
    for ( uint j = 0; j < 3; j++ )
    {
        bol1 = key[N - 1 - j];
        // cout<<"j direction"<<bol1<<endl;

        kt1 = 0;
        // if bol1=false basically the generated auxillary code will
        if ( bol1 == true )
        {
            for ( uint k = 0; k < mylevel; k++ )
            {
                kt1[N - 1 - 3 * k - j] = bol1;
            }
        }
        kt2 = kt1;

        if ( ( kt1 & key ) == kt2 )
        {
            bol2 = true;
            //     cout << bol2 << endl;
            break;
        }
    }

    return ( bol2 );
}

template <size_t N, typename value>
std::pair<morton<N>, int> Tree<N, value>::readRefineList( typename std::unordered_map<morton<N>, int>::iterator it )
{
    /*
        if ( i >= refinelist.size() )
        {
            throw std::runtime_error( RED "wrong index for refine list" RESET );
        }
        key = refinelist.at( i );
    */

    return ( std::make_pair( it->first, it->second ) );
}

template <size_t N, typename value>
void Tree<N, value>::flipRefineElemTag( typename std::unordered_map<morton<N>, int>::iterator it )
{

    it->second = 1;
}
/*
template <size_t N, typename value>
typename std::unordered_map<morton<N>,int>::iterator  Tree<N, value>::readDerefineList( morton<N> key)
{

    return(derefinelist.find(key) );
}
*/
template <size_t N, typename value>
morton<N> Tree<N, value>::readDerefineList( typename std::unordered_map<morton<N>, int>::iterator it )
{

    return ( it->first );
}

template <size_t N, typename value>
bool Tree<N, value>::isBoundary( const morton<N> &key, uint direction )
{
    uint mylevel;
    bool bol1, bol2;
    morton<N> kt1, kt2;
    level( key, &mylevel );
    bol2 = true;
    // cout<<"key= "<<key<<"level "<<mylevel<<endl;
    // need to chec for all three directions for every element

    uint j = direction;

    bol1 = key[N - 1 - j];
    // cout<<"j direction"<<bol1<<endl;
    for ( uint i = 1; i < mylevel; i++ )
    {
        if ( key[N - 3 * i - 1 - j] != bol1 )
        {
            bol2 = false;
            break;
        }
    }

    return ( bol2 );
}

// This is to identify the corner elements to find out if a boundary is a corner element
//  or side element, bottom line number of directions that an element has a nonlocal nbr
template <size_t N, typename value>
void Tree<N, value>::getDirections( morton<N> &key, vector<uint> &directions )
{
    directions.clear();
    uint mylevel;
    bool bol1, bol2;
    morton<N> kt1, kt2;
    level( key, &mylevel );
    bol2 = false;
    // excluding when the level is zero
    // since this is an exception as the other elements are siblings not nonlocal neighbor
    // cout<<"key= "<<key<<"level "<<mylevel<<endl;
    // need to chec for all three directions for every element
    if ( mylevel != 1 )
    {
        for ( uint j = 0; j < 3; j++ )
        {
            bol1 = key[N - 1 - j];
            // cout<<"j direction"<<bol1<<endl;
            kt1 = 0;
            // if bol1=false basically the generated auxillary code will be zero
            if ( bol1 == true )
            {
                for ( uint k = 0; k < mylevel; k++ )
                {
                    kt1[N - 1 - 3 * k - j] = bol1;
                }
            }
            kt2 = kt1;
            if ( ( kt1 & key ) == kt2 )
            {
                bol2 = true;
                //         cout << bol2 << endl;
                directions.push_back( j );
                // break;
            }
        }
    }
}

template <size_t N, typename value>
void Tree<N, value>::extractBoundary()
{
    morton<N> key, kt1, kt2;

    // empty the list before assigning

    // boundaryElem.clear();

    bool bol1, bol2;
    uint mylevel = 0;

    // cout<<"list size"<<refinelist.size()<<endl;

    for ( uint i = 0; i < refinelist.size(); i++ )
    {
        key = refinelist.at( i );

        // cout << "bol " << bol2 << endl;
        if ( isBoundary( key ) )
        {
            //    cout << key << endl;
            boundarylist.push_back( key );
        }
    }
}

template <size_t N, typename value>
void Tree<N, value>::extractBoundaryP( uint istart, uint iend )
{
    morton<N> key, kt1, kt2;

    // empty the list before assigning

    // boundaryElem.clear();

    bool bol1, bol2;
    uint mylevel = 0;

    // cout<<"list size"<<refinelist.size()<<endl;

    //  cout << "istart " << istart << "iend " << iend << endl;

    for ( uint i = istart; i < iend; i++ )

    {
        key = refinelist.at( i );

        //  cout << "key " << key << endl;
        if ( isBoundary( key ) )
        {
            //   cout << "boundary key " << key << endl;
            boundarylist.push_back( key );
        }
    }
}

template <size_t N, typename value>
bool Tree<N, value>::isInMeshList( const morton<N> &key )
{
    bool bol;
    bol = false;

    if ( mesh.find( key ) != mesh.end() )
    {
        bol = true;
    }

    return ( bol );
}

template <size_t N, typename value>
bool Tree<N, value>::isInRefineList( const morton<N> &key )
{
    bool bol;
    bol = false;
    if ( refinelist.count( key ) != 0 )
    {
        bol = true;
    }
    return ( bol );
}

template <size_t N, typename value>
void Tree<N, value>::constructHigherLevelNbrs( const morton<N> &key, const uint &keylevel, const uint &direction, morton<N> *nbr )
{
    static const uint index[12] = {0, 1, 2, 1, 1, 0, 2, 0, 2, 0, 1, 0};

    if ( N < 3 * ( keylevel + 1 ) )
    {
        throw std::runtime_error( "level is too big tp fit in the given morton code" );
    }

    morton<N> kt;
    kt = key;
    // generate one level higher in the tree

    // repeat the last three bits
    for ( uint i = 0; i < 3; i++ )
    {
        kt[N - 1 - 3 * ( keylevel ) - i] = kt[N - 1 - 3 * ( keylevel - 1 ) - i];
    }
    // cout<<key<<endl;
    // cout<<direction<<endl;

    // careful implementation eliminates one flip
    for ( uint i = 0; i < 4; i++ )
    {
        kt.flip( N - 1 - 3 * ( keylevel ) - index[4 * direction + i] );
        nbr[i] = kt;
        // cout<<GREEN<<kt<<RESET<<endl;
    }
}

template <size_t N, typename value>
void Tree<N, value>::constructNonlocalHigherLevelNbrs( const morton<N> &key, const uint &keylevel, const uint &direction, morton<N> *nbr )
{
    static const uint index[9] = {1, 2, 1, 0, 2, 0, 0, 1, 0};

    if ( N < 3 * ( keylevel + 1 ) )
    {
        throw std::runtime_error( "level is too big tp fit in the given morton code" );
    }

    morton<N> kt;
    kt = key;
    // generate one level higher in the tree

    // repeat the last three bits
    for ( uint i = 0; i < 3; i++ )
    {
        kt[N - 1 - 3 * ( keylevel ) - i] = kt[N - 1 - 3 * ( keylevel - 1 ) - i];
        nbr[0] = kt;
    }
    // cout<<key<<endl;
    // cout<<direction<<endl;

    // careful implementation eliminates one flip
    for ( uint i = 0; i < 3; i++ )
    {
        kt.flip( N - 1 - 3 * ( keylevel ) - index[3 * direction + i] );
        nbr[i + 1] = kt;
        // cout<<GREEN<<kt<<RESET<<endl;
    }
}

template <size_t N, typename value>
void Tree<N, value>::printMesh()
{

    // bitmap<N,typename value*>::hasher fn = mesh.hash_function();

    for ( auto it = mesh.begin(); it != mesh.end(); it++ )
    {
       // cout << it->first << "   " << it->first.to_ulong() << endl;
        // cout<<fn(it->first)<<endl;
    }

    for ( auto &x : mesh )
    {
        //std::cout << "Element [" << x.first << ":" << x.second << "]";
        //std::cout << " is in bucket #" << mesh.bucket( x.first ) << "  " << x.first.to_ulong() << std::endl;
    }
}

template <size_t N, typename value>
void Tree<N, value>::getKey( uint i, morton<N> &key )
{

    auto it = std::next( mesh.begin(), i );

    key = it->first;
}

template <size_t N, typename value>
void Tree<N, value>::clearMortonSTL()
{
    mortonSTL.clear();
}

template <size_t N, typename value>
typename unordered_map<morton<N>, int>::iterator Tree<N, value>::Rbegin()
{
    return ( refinelist.begin() );
}

template <size_t N, typename value>
typename unordered_map<morton<N>, int>::iterator Tree<N, value>::Rend()
{
    return ( refinelist.end() );
}

//*************************************************************************
//
//                            Derefinement Routines
//
//
//*************************************************************************
// derefinelist is private I need to access it by forest

template <size_t N, typename value>
typename unordered_map<morton<N>, int>::iterator Tree<N, value>::Dbegin()
{
    return ( derefinelist.begin() );
}

template <size_t N, typename value>
typename unordered_map<morton<N>, int>::iterator Tree<N, value>::Dend()
{
    return ( derefinelist.end() );
}

//======================================================================================
/**
*  \brief If any of the siblings are listed in the dereffinement do not add to the list
*  as derefining one child means removing all the siblings
*
*
* */

#if ( 1 )
template <size_t N, typename value>
void Tree<N, value>::addToDerefineList( morton<N> key )
{
    bool bol = true;

    uint mylevel, klevel;
    morton<N> sibkey[7], kt;
    level( key, &mylevel );
    siblings( key, mylevel, sibkey );

    // if any of the siblings has a higher level we can not remove that element
    //
    for ( uint i = 0; i < 7; i++ )
    {
        level( sibkey[i], &klevel );
        if ( klevel > mylevel )
        {
            bol = false;
            break;
        }
    }

    if ( bol && derefinelist.count( key ) == 0 )
    {
        derefinelist.insert( {key, 0} );

        for ( uint i = 0; i < 7; i++ )
        {
            derefinelist.insert( {sibkey[i], 0} );
        }
    }
}
#endif

#if ( 1 )
template <size_t N, typename value>
void Tree<N, value>::pushToDerefinelist( uint nlevel ) /**!< this function is to find a value given the key*/
{
    morton<N> key, key1, kt;

    key = 0;
    // cout<<mortonSTL.size()<<endl;
    uint mylevel;
    morton<N> sibkey[7];
    std::unordered_set<morton<N>> temp;

    // first clear the derefinelist
    //
    //   derefinelist.clear();

    // temp is to speed up the search, mortonSTL is the new geometry encoded based on the new location (moving body)
    /*
    save a temporary set of the mortonSTL (temp) with one level lower than the actual nlevel
    if that element is not inside mortonSTL, insert it to temp

    */
    // look at the parent element, thats why I do (nlevel-1)
    //
    //
    //   cout << "stlsize" << mortonSTL.size() << endl;
    for ( uint i = 0; i < mortonSTL.size(); i++ )
    {
        key = mortonSTL.at( i );

        for ( uint j = 0; j < 3 * ( nlevel - 1 ); j++ )
        {
            key1[N - j - 1] = key[N - j - 1];
        }
        if ( temp.count( key1 ) == 0 )
        {
            temp.insert( key1 );
            //cout<<key1<<endl;
        }
    }

    for ( auto it = mesh.begin(); it != mesh.end(); it++ )
    {

        key = it->first;
        level( key, &mylevel );
        kt = key;

        kt[N - 3 * ( mylevel - 2 ) - 1] = 0;
        kt[N - 3 * ( mylevel - 2 ) - 2] = 0;
        kt[N - 3 * ( mylevel - 2 ) - 3] = 0;

        if ( mylevel == nlevel && temp.count( kt ) == 0 )
        {
            // insert this element and all its siblings

            addToDerefineList( key );
        }
    }
    temp.clear();

//    cout << derefinelist.size() << endl;

//    cout << RED "stl size " << mortonSTL.size() << RESET << endl;
}
#endif

template <size_t N, typename value>
void Tree<N, value>::retainFourToOne() /*!< imposes 4:1 balance given the list of elments to be refined in the vector refine list*/
{
    morton<N> kt;
    uint mylevel, nbrlevel, a, changedirectionlevel;
    morton<N> mykey, sibkey[7];
    //    uint      istart = 0;
    //    uint      iend   = refinelist.size();
    //      printf("istart %d iend %d\n",istart,iend);

    //
    // check the balance for nonlocal neighbors
    //
    for ( auto it = derefinelist.begin(); it != derefinelist.end(); it++ )
    {
        mykey = ( it->first );
        level( mykey, &mylevel );

        if ( mylevel > 1 )
        {
            for ( uint j = 0; j < 3; j++ )
            {
                //      uint j=1;

                kt = mykey;

                //      cout<<"mykey\t"<<kt<<endl;
                findFlipLevel( kt, &mylevel, &changedirectionlevel, &j );

                // if the change in signe does not happen that is a boundary cube
                //       cout<<"changedirectionlevel"<<changedirectionlevel<<endl;
                if ( changedirectionlevel != 0 )
                {
                    flipForNbr( &kt, &mylevel, &changedirectionlevel, &j );
                    mesh.count( kt );

                    // if this element exists, the level of nbr>=level of tagged element

                    //      cout<<"kt="<<mesh.count(kt)<<endl;
                    if ( mesh.count( kt ) == 0 )
                    {
                        kt[N - 3 * ( mylevel - 1 ) - 1] = 0;
                        kt[N - 3 * ( mylevel - 1 ) - 2] = 0;
                        kt[N - 3 * ( mylevel - 1 ) - 3] = 0;
                    }
                    //      cout<<"modified key\t"<<kt<<endl;
                    level( kt, &nbrlevel );

                    // cout<<"nbr level"<<nbrlevel<<"my level"<<mylevel<<endl;

                    // if its level is lower than the neighbor,remove it from the derefine list
                    if ( mylevel < nbrlevel )
                    {
                        // remove the element and all its siblings from the list
                        //                        cout<<mylevel<<" "<<nbrlevel<<endl;
                        removeFromDerefineList( it );
                    }
                }
            }
        }
    }

    cout << derefinelist.size() << endl;
}
template <size_t N, typename value>
typename std::unordered_map<morton<N>, int>::iterator Tree<N, value>::findInDerefine( morton<N> key ) /*!< E;iminate from derefinelist*/
{
    return ( derefinelist.find( key ) );
}

template <size_t N, typename value>
void Tree<N, value>::removeFromDerefineList( typename std::unordered_map<morton<N>, int>::iterator it ) /*!< E;iminate from derefinelist*/
{
    // this is to remove from the list so removes all 8 element
    uint mylevel;
    morton<N> key, sibkey[8];
    key = it->first;
    level( key, &mylevel );
    siblings( key, mylevel, sibkey );
    sibkey[7] = key;
    // find the key and its siblings and put int to 1

    for ( uint i = 0; i < 8; i++ )
    {
        auto it = derefinelist.find( sibkey[i] );
        if ( derefinelist.count( sibkey[i] ) != 0 )
        {
            it->second = 1;
        }
    }
}

//===========================================================
//
//         derefine: Eliminate An Element From Mesh
//
//===========================================================
template <size_t N, typename value>
void Tree<N, value>::derefine( morton<N> key )
{
    /*! \brief if the morton code does not exist in mesh, refinement is not permitted (derefining a nonexsiting element not permitted)
     *  Also, if any of the siblings have a higher level of refinement, derefinement is ignored*/

    // cout<<"key is :"<<key<<endl;

    /*!< \brief if the key does not exist simply igonre doing anything*/
    morton<N> sibkey[7], kt;
    uint mylevel;

    level( key, &mylevel );

    siblings( key, mylevel, sibkey );

    if ( mesh.count( key ) != 0 )
    { // if the siblings

        for ( uint i = 0; i < 7; i++ )
        {
            mesh.erase( sibkey[i] );
            //            cout<<BLUE<<sibkey[i]<<endl;
        }

        mesh.erase( key );
        kt = key;
        kt[N - 3 * ( mylevel - 1 ) - 1] = 0;
        kt[N - 3 * ( mylevel - 1 ) - 2] = 0;
        kt[N - 3 * ( mylevel - 1 ) - 3] = 0;
        mesh.insert( {kt, nullptr} );
    }
    else
    {
        cout << RED "derefinement ignored, mesh not in the list" RESET << endl;
    }
}

template <size_t N, typename value>
void Tree<N, value>::derefineDerefineList()
{
    morton<N> key;
    bool bol;
    bitvector<N> temp;
    morton<N> sibkey[8];

    uint mylevel;
    //    uint maxlevel = nlevel;
    uint changedirectionlevel, nbrlevel;
    morton<N> kt;

    //cout << " derefine list " << derefinelist.size() << endl;

    retainFourToOne();

    for ( auto i = derefinelist.begin(); i != derefinelist.end(); i++ )
    {
        key = i->first;
        //        cout<<key<<" "<<i->second<<endl;
        if ( i->second == 0 )
        {
            derefine( key );
            // cout<<RED<<key<<RESET<<endl;
            removeFromDerefineList( i );
        }
    }
    derefinelist.clear();
}

template <size_t N, typename value>
void Tree<N, value>::clearMesh() /*!< Eliminate from derefinelist*/
{
    mesh.clear();
}

template <size_t N, typename value>
void Tree<N, value>::mortonSTLclear() /*!< Eliminate from derefinelist*/
{

    mortonSTL.clear();
    mortonSTL.shrink_to_fit();
}

// insert single element, this fucntion is to developed for weak scaling analysis
// where each processor has one seed

template <size_t N, typename value>
void Tree<N, value>::insertSeed( morton<N> &key ) /*!< Eliminate from derefinelist*/
{

    // not that constructing the object will put the root as 0
    // for weak scaling we only need one seed given by key as input to this function

    mesh.erase( 0 );
    mesh.insert( {key, nullptr} );
    // cout<<key<<"   "<<mesh.size()<<endl;
    // auto it=mesh.begin();
}

template <size_t N, typename value>
void Tree<N, value>::insertNbrs( vector<int> &Nbrs ) /*!< add neighbors*/
{

    for ( uint i = 0; i < Nbrs.size(); i++ )
    {
        morton<N> key( Nbrs.at( i ) );

        mesh.insert( {key, nullptr} );
    }
}

template <size_t N, typename value>
void Tree<N, value>::getCoords( real *X ) /*!< add neighbors*/
{

      X[0]=ancestorcoords[0]-0.5*ancestorlength[0];
      X[1]=ancestorcoords[0]+0.5*ancestorlength[0];

      X[2]=ancestorcoords[1]-0.5*ancestorlength[1];
      X[3]=ancestorcoords[1]+0.5*ancestorlength[1];
 
      X[4]=ancestorcoords[2]-0.5*ancestorlength[2];
      X[5]=ancestorcoords[2]+0.5*ancestorlength[2];

}

template <size_t N, typename value>
bool Tree<N, value>::operator==(const  Tree<N, value>  & T0) const /*!< add neighbors*/
{

       
       bool bol0=fabs(T0.ancestorcoords[0]-this->ancestorcoords[0])<1e-6;
       bool bol1=fabs(T0.ancestorcoords[1]-this->ancestorcoords[1])<1e-6;
       bool bol2=fabs(T0.ancestorcoords[2]-this->ancestorcoords[2])<1e-6;
     
       bool bol3=fabs(T0.ancestorlength[0]-this->ancestorlength[0])<1e-6;
       bool bol4=fabs(T0.ancestorlength[1]-this->ancestorlength[1])<1e-6;
       bool bol5=fabs(T0.ancestorlength[2]-this->ancestorlength[2])<1e-6;

       bol0=bol0 && bol1 && bol2 && bol3 && bol4 && bol5;

      return(bol0);
}


#if ( 0 )
template <size_t N, typename Nvalue, size_t N1>
void convertCoordToMorton( Tree<N, Nvalue> T, real *xyz, morton<N> &key ) /**!< this function is to find a value given the key*/
{
    real x, y, z, xmid, ymid, zmid;

    static const real half = 0.5;

    real xmin = T.ancestorcoords[0] - half * T.ancestorlength[0];
    real xmax = T.ancestorcoords[0] + half * T.ancestorlength[0];

    real ymin = T.ancestorcoords[1] - half * T.ancestorlength[1];
    real ymax = T.ancestorcoords[1] + half * T.ancestorlength[1];

    real zmin = T.ancestorcoords[2] - half * T.ancestorlength[2];
    real zmax = T.ancestorcoords[2] + half * T.ancestorlength[2];

    x = xyz[0];
    y = xyz[1];
    z = xyz[2];

    key = 0;
#if ( 0 )
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
#endif
}

#endif

template class Tree<TREESIZE, uint>;
template class Tree<TREESIZE, real>;


template class Tree<64, uint>;
template class Tree<64, real>;

template class Tree<16, uint>;
template class Tree<16, real>;
template class Tree<3, uint>;
template class Tree<3, real>;
template class Tree<6, uint>;
template class Tree<9, uint>;
template class Tree<6, real>;
template class Tree<9, real>;

/*
template class Tree<PROCSIZE, uint>;
template class Tree<PROCSIZE, real>;
*/
//TheTemplate<conditional<is_same<Type2, Type4>::value, dummy<4>, Type4>::type>



#if(1)
//template class Tree<conditional<!is_same< Tree<PROCSIZE,uint>, Tree<TREESIZE,uint> > ::value,Tree<5,uint> >;
//template class Tree<!is_same< Tree<PROCSIZE,uint>, Tree<TREESIZE,uint> > ::value,Tree<> >;
//template class  Tree<conditional<(!is_same< Tree<PROCSIZE,uint>, Tree<TREESIZE,uint> > ::value) , Tree<10,uint>,Tree<PROCSIZE,uint>>::type>;

//template class Tree<enable_if<!is_same<TREESIZE,unit, PROCSIZE,uint>::value>,Type3>::type>, uint>;
//template class Tree<TREESIZE, real>;
#endif

//template class Tree<PROCSIZE, real>;


/*
template class Tree<48, real>;
template class Tree<128, real>;
*/
// template class Tree<17, uint>;
