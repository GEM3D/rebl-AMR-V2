#include "tree.h"
//======================================================================
//
//               initializa geometryVoxel object
//
//=====================================================================
template <size_t N, typename value>
void Voxel<N, value>::setLevel( uint *l )
{

    maxlevel = *l;
    real v;
    morton<N> ancestorkey = 0;
    mesh.insert( {ancestorkey, &v} );

    // cout<<maxlevel<<endl;
    if ( (uint)N / maxlevel < 3 )
    {
        cout << RED "Not enough bits for specified level" RESET << endl;
        exit( 0 );
    }

    morton<N> kt = 0;

    for ( uint i = 0; i < maxlevel; i++ )
    {

        kt.flip( N - 1 - 3 * i );
        kt.flip( N - 1 - 3 * i - 1 );
        kt.flip( N - 1 - 3 * i - 2 );

        lookup.push_back( kt );
        // cout<<kt<<endl;
    }
}

template <size_t N, typename value>
Voxel<N, value>::~Voxel() /* !< class Destructor*/
{

    // typename bitmap<N,value>::iterator it;

    for ( auto it = Tree<N, value>::begin(); it != Tree<N, value>::end(); it++ )
    {
        /*
        if(it->second!=nullptr)
        {
        free(it->second);
        }
        */
    }
}

//
//
// generate search tree for refinement and derefinement
//
//

template <size_t N, typename value>
void Voxel<N, value>::generateSearchTree( real *geom_xyz, uint geom_nn )
{
    bitvector<N> mylist;
    uint mylevel;

    // typename bitmap<N,value>::iterator local_it;

    for ( uint j = 0; j < maxlevel; j++ )
    {

        for ( auto local_it = Tree<N, value>::begin(); local_it != Tree<N, value>::end(); local_it++ )
        {

            Tree<N, value>::level( local_it->first, &mylevel );
            // cout<<mylevel<<endl;
            if ( Tree<N, value>::isInsideSolid( local_it->first, geom_xyz, geom_nn ) == 1 && mylevel < maxlevel )
            {
                mylist.push_back( local_it->first );
            }
        }

        cout << RED "list size= " RESET << mylist.size() << endl;

        if ( mylist.size() == 0 )
        {
            break;
        }

        for ( auto it = 0; it < mylist.size(); it++ )
        {
            Tree<N, value>::refine( mylist.at( it ) );
        }

        mylist.clear();
        mylist.shrink_to_fit();
    }
}

//
//
// \brief this function allocates the values using the tree
//
//

template <size_t N, typename value>
void Voxel<N, value>::distributeGeomToLeaves( real *geom_xyz, uint n )
{
    typename bitmap<N, real>::iterator kt;
    real xyz[6];
    morton<N> key;

    uint count = 0;

    real *ppointer;

    // typename bitmap<N,value>::iterator it;

    for ( auto it = Tree<N, value>::begin(); it != Tree<N, value>::end(); it++ )
    {
        // cout<<" inside the func" <<endl;

        key = it->first;

        Tree<N, value>::enclosingBox( key, xyz );

        // delete [] it->second;

        // it->second=(real*)realloc(it->second,sizeof(real));

        // I intentionally do not use vectors, due to inehritance and other complications that comes with it
        //  use classical real pointer that is more general, can template bitmap to implement vectors in the future

        count = 0;

        for ( uint j = 0; j < n; j++ )
        {
            if ( geom_xyz[3 * j + 0] >= xyz[0] && geom_xyz[3 * j + 0] <= xyz[1] )
            {
                if ( geom_xyz[3 * j + 1] >= xyz[2] && geom_xyz[3 * j + 1] <= xyz[3] )
                {
                    if ( geom_xyz[3 * j + 2] >= xyz[4] && geom_xyz[3 * j + 2] <= xyz[5] )
                    {
                        count++;
                    }
                }
            }
        }

        // save the size of the array in the first location

        // it->second=new real[3*count+1];
        if ( count > 0 )
        {
            it->second = (real *)realloc( it->second, ( 3 * count + 1 ) * sizeof( real ) );
            // cout<<(it->second)<<endl;

            it->second[0] = count;

            // cout<<count<<endl;
            // cout<<"count"<<(it->second)[0]<<endl;

            count = 0;

            for ( uint j = 0; j < n; j++ )
            {
                if ( geom_xyz[3 * j + 0] >= xyz[0] && geom_xyz[3 * j + 0] <= xyz[1] )
                {
                    if ( geom_xyz[3 * j + 1] >= xyz[2] && geom_xyz[3 * j + 1] <= xyz[3] )
                    {
                        if ( geom_xyz[3 * j + 2] >= xyz[4] && geom_xyz[3 * j + 2] <= xyz[5] )
                        {
                            it->second[3 * count + 1] = geom_xyz[3 * j + 0];
                            it->second[3 * count + 2] = geom_xyz[3 * j + 1];
                            it->second[3 * count + 3] = geom_xyz[3 * j + 2];

                            count++;
                        }
                    }
                }
            }
        }
    }
    // cout<<sizeof(mesh)<<endl;
}

//
// this routine is used to shorten the derefine geometry
//
// it basically checks to see if foir a given key, are the siblings
// from the same level and also if any of them include any geometry points
template <size_t N, typename value>
uint Voxel
<N, value>::checkSiblingStatus( morton<N> key,
                                morton<N> *sibkey ) /*!< checks to see if siblings include any points and whether they have the same level*/
{
    uint mylevel, siblevel;
    Tree<N, value>::level( key, &mylevel );
    Tree<N, value>::siblings( key, mylevel, sibkey );
    uint bol = 1;
    typename bitmap<N, value>::iterator kt;
    uint size;

    for ( uint m = 0; m < 7; m++ )
    {
        kt = Tree<N, value>::find( sibkey[m] );
        //
        // modify this level check only for 0000000000000 sibling
        //
        if ( kt != Tree<N, value>::end() )
        {
            Tree<N, value>::level( sibkey[m], &siblevel );
            // size=kt->second[0];
            // cout<<size<<endl;

            if ( siblevel != mylevel || kt->second != nullptr )
            {
                bol = 0;
                // break;
            }
            // cout<<"sibs\t"<<sibkey[m]<<endl;
            // cout<<"sib level\t"<<siblevel<<"size\t"<<size<<endl;
        }
    }
    return ( bol );
}

//============================================================
template <size_t N, typename value>
void Voxel<N, value>::derefineGeomTree()
{
    bitvector<N> mylist;
    uint mylevel, siblevel;
    morton<N> key, sibkey[7];

    typename bitmap<N, value>::iterator kt, it;
    uint bol = 0;
    /*
    cout<<"----------------------------\n"<<endl;
    for(auto it=begin();it!=end();it++)
    {
    cout<<it->first<<":\t"<<it->second[0]<<endl;
    }

    cout<<"---------------------------\n"<<endl;

    */
    uint size;
    uint count = 0;
    for ( uint i = maxlevel; i > 1; i-- )
    {
        // uint i=maxlevel;
        for ( it = Tree<N, value>::begin(); it != Tree<N, value>::end(); it++ )
        {
            // auto it=begin();
            key = it->first;
            Tree<N, value>::level( key, &mylevel );
            bol = 0;
            Tree<N, value>::siblings( key, mylevel, sibkey );

            if ( mylevel == i && it->second == nullptr && checkSiblingStatus( key, sibkey ) )
            {
                // cout<<"null "<<it->second<<endl;
                bol = 1;

                for ( uint k = 0; k < mylist.size(); k++ )
                {
                    for ( uint l = 0; l < mylist.size(); l++ )
                    {
                        if ( mylist.at( k ) == mylist.at( l ) && k != l )
                        {
                            bol = 0;
                            break;
                        }
                    }
                }

                for ( uint k = 0; k < mylist.size(); k++ )
                {
                    for ( uint l = 0; l < 7; l++ )
                    {
                        if ( mylist.at( k ) == sibkey[l] )
                        {
                            bol = 0;
                            break;
                        }
                    }
                }

                count++;
                // cout<<"bol="<<bol<<endl;
            }

            if ( bol == 1 )
            {
                cout << "bol =" << bol << endl;
                mylist.push_back( key );
            }
        }
        cout << "derefine list size" << mylist.size() << endl;

        // derefine those elements frist

        for ( uint j = 0; j < mylist.size(); j++ )
        {
            Tree<N, value>::derefine( mylist.at( j ) );
        }

        mylist.clear();
        mylist.shrink_to_fit();
        cout << "level\t" << i << endl;
    }
    cout << "count=" << count << endl;
}

#if ( 1 )

template <size_t N, typename value>
bool Voxel<N, value>::IsInsideSegment( morton<N> key, real *xyz )
{
    morton<N> kt;
    bool e = 1;
    uint l = maxlevel - 1;
    // cout<<"\n"<<key<<endl;

    while ( e )
    {
        kt = lookup.at( l ) & key;

        // cout<<"\n"<<kt<<endl;

        // cout<<"\n"<<count(kt)<<endl;

        if ( mesh.count( kt ) == 1 )
        {
            e = 0;
        }
        else
        {
            l--;
        }

        if ( l == -1 )
        {
            cout << RED "fatal error: in finding parent elment in voxel IsInsideSegment" RESET << endl;
            exit( 0 );
        }
    }

    typename bitmap<N, value>::iterator it;
    // real xyz[6];

    // enclosing box is a problem, it calculates the box for voxel not Mesh

    // enclosingBox(key,xyz);
    bool bol = 0, a, b, c;
    it = Tree<N, value>::find( kt );

    if ( it->second != nullptr )
    {
        for ( uint j = 0; j < it->second[0]; j++ )
        {
            a = it->second[3 * j + 1] >= xyz[0] && it->second[3 * j + 1] <= xyz[1];
            b = it->second[3 * j + 2] >= xyz[2] && it->second[3 * j + 2] <= xyz[3];
            c = it->second[3 * j + 3] >= xyz[4] && it->second[3 * j + 3] <= xyz[5];

            if ( a && b && c )
            {
                bol = 1;
                break;
            }
        }
    }
    return ( bol );
}

#endif

template class Voxel<32, real>;

template class Voxel<16, real>;
