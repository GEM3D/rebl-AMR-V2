#include "templateForest.h"
#include "definitions.h"
#include <string.h>


// Methods 0 nd 1 are obsolete, not portable due to byte transfer
// METHOD 0 >>  uses bool for symmetric Comm
// METHOD 1 >>  uses character  for Symmetric Comm
// METHOD 2 >>  uses a cleaned-up version of character for Symmetric Comm
// METHOD 3 >>  specialized for weak scaling
// METHOD 4 >>  uses MPI 3.0 neighbor collectives
// -- The neighbor collective is not easy to implemenet unless the graph
// is updated and this is alot expensive,
// A better approach, instead of trying to eliminate broadcast, is to
// take advantage of communication/computation overlap

// template <size_t N, typename Nvalue, size_t M, typename Mvalue,  template< size_t L, typename Lvalue> class T >
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
TemplateForest<N, Nvalue, M, Mvalue, T>::TemplateForest(int argcs,char *pArgs[], T &proc, real *length, real *coords, uint nx, uint ny, uint nz )
{

/*!<part I, initialize the ancestor coords and length, just like we did for class tree*/

    for ( uint i = 0; i < 3; i++ )
    {
        ancestorlength[i] = length[i];
        ancestorcoords[i] = coords[i];
    }

    real X[6], XC[3], len[3];

    int initFlag=0;

    if ( MPI_Initialized( &initFlag ) != MPI_SUCCESS  && Com.myrank==0 )
    {
        cout << " Exit Code : " << ReblAmrGetErrorEnum( MPI_INIT_CHECK_FAIL ) << endl;
        exit(1 );
    }

    if ( initFlag == 0 )
    {
        if ( MPI_Init( &argcs, &pArgs ) != MPI_SUCCESS  && Com.myrank==0  )
        {
            cout << " Exit Code : " << ReblAmrGetErrorEnum( MPI_INIT_FAIL ) << endl;
            exit( 1 );
        }
    }

    MPIStartUp();
    
    checkInputParams(argcs,pArgs);

    npx = nx;
    npy = ny;
    npz = nz;

    if ( npx == 0 )
    {
        npx = 2;
    }
    if ( npy == 0 )
    {
        npy = 2;
    }
    if ( npz == 0 )
    {
        npz = 2;
    }

    assignSeeds( length, proc );

    geom.construct( ancestorlength, ancestorcoords, 2, 2, 2 );

//    geom.setToZero(  );

    // this error is important as to Remove singualt bits I tag that bit using one bit from far right hand side, i.e. key[0]=1

    if ( (M + N) % 3 == 0 || M > N )
    {
      cout <<RED<<ReblAmrGetErrorEnum(  COMBINED_SIZE ) << RESET<<endl;
      exit( COMBINED_SIZE); 
 //       throw std::runtime_error( RED "(M+N) can not be a multiply of 3 since  one bit is needed to recover singularity in parallel" RESET );
    }
    if ( ZOLTAN_ON )
    {
        zz = Zoltan_Create( Com.mpicom );
    }
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::construct(int argcs,char *pArgs[], T &proc, real *length, real *coords, uint nx, uint ny, uint nz )
{

/*!<part I, initialize the ancestor coords and length, just like we did for class tree*/

    for ( uint i = 0; i < 3; i++ )
    {
        ancestorlength[i] = length[i];
        ancestorcoords[i] = coords[i];
    }

    real X[6], XC[3], len[3];

    int initFlag=0;

    if ( MPI_Initialized( &initFlag ) != MPI_SUCCESS  && Com.myrank==0 )
    {
        cout << " Exit Code : " << ReblAmrGetErrorEnum( MPI_INIT_CHECK_FAIL ) << endl;
        exit(1 );
    }

    if ( initFlag == 0 )
    {
        if ( MPI_Init( &argcs, &pArgs ) != MPI_SUCCESS  && Com.myrank==0  )
        {
            cout << " Exit Code : " << ReblAmrGetErrorEnum( MPI_INIT_FAIL ) << endl;
            exit( 1 );
        }
    }

    MPIStartUp();
    
    checkInputParams(argcs,pArgs);

    npx = nx;
    npy = ny;
    npz = nz;

    if ( npx == 0 )
    {
        npx = 2;
    }
    if ( npy == 0 )
    {
        npy = 2;
    }
    if ( npz == 0 )
    {
        npz = 2;
    }

    assignSeeds( length, proc );

    cout<<"done"<<endl;
    geom.construct( ancestorlength, ancestorcoords, 2, 2, 2 );

    // this error is important as to Remove singualt bits I tag that bit using one bit from far right hand side, i.e. key[0]=1

    if ( (M + N) % 3 == 0 || M > N )
    {
      cout <<RED<<ReblAmrGetErrorEnum(  COMBINED_SIZE ) << RESET<<endl;
      exit( COMBINED_SIZE); 
 //       throw std::runtime_error( RED "(M+N) can not be a multiply of 3 since  one bit is needed to recover singularity in parallel" RESET );
    }
    if ( ZOLTAN_ON )
    {
        zz = Zoltan_Create( Com.mpicom );
    }

}


template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::checkInputParams(int argcs,char *pArgs[])
{
// check input parameters
//

if ( (argcs < 4) && (Com.myrank==0) )
{
  cout <<RED<<ReblAmrGetErrorEnum( NUM_INPUT_ARGS )  << RESET<<endl;
  exit( NUM_INPUT_ARGS );
}

 int proclevel = atoi(pArgs[2]);
 int meshlevel = atoi(pArgs[3]);

  if ( proclevel > ( PROCSIZE / 3 ) &&  Com.myrank==0   )
   {
      cout <<RED<<ReblAmrGetErrorEnum( PROC_LEVEL ) << RESET<<endl;
      exit(PROC_LEVEL); 
  }
  if ( meshlevel > ( TREESIZE / 3 ) &&  Com.myrank==0  )
   {
      cout <<RED<<ReblAmrGetErrorEnum( MESH_LEVEL ) << RESET<<endl;
      exit(MESH_LEVEL); 
  }

}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
TemplateForest<N, Nvalue, M, Mvalue, T>::~TemplateForest()
{

    if ( message != nullptr )
    {
        delete[] message;
    }

    if ( request != nullptr )
    {
        delete[] request;
    }

    if ( request1 != nullptr )
    {
        delete[] request1;
    }
    if ( ZOLTAN_ON )
    {
        Zoltan_Destroy( &zz );
    }
   
    int initFlag;

// free the duplicated MPI_COMM

if(duplicated==1)
{
    MPI_Comm_free(&Com.mpicom);
}
    if ( MPI_Finalized( &initFlag ) != MPI_SUCCESS  && Com.myrank==0 )
    {
        cout << "failure in checking if MPI has already initialize" << endl;
        exit( 1 );
    }

    // cout<<" init flag "<<initFlag<<endl;
    if ( initFlag == 0 )
    {
        if ( MPI_Finalize() != MPI_SUCCESS  && Com.myrank== 0  )
        {
            cout << " Exit Code : " << ReblAmrGetErrorEnum( MPI_FINALIZE_FAIL ) << endl;
            exit( 1 );
        }
    }


}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::assignSeeds( real *length, T &proc )
{

    uint count = 0;
    real X[6], XC[3], len[3];

    for ( auto it = proc.begin(); it != proc.end(); it++ )
    {
        //    cout<<" my_rank "<<Com.myrank <<"  "<<it->first<<" "<<it->second[0]<<endl;
        if ( it->second[0] == Com.myrank )
        {
            proc.enclosingBox( it->first, X );

            XC[0] = ( X[0] + X[1] ) * 0.5;
            XC[1] = ( X[2] + X[3] ) * 0.5;
            XC[2] = ( X[4] + X[5] ) * 0.5;

            len[0] = fabs( X[1] - X[0] );
            len[1] = fabs( X[3] - X[2] );
            len[2] = fabs( X[5] - X[4] );

            // regrdless of the topology, we need a tree in tree list

            trees.push_back( Tree<N, Nvalue>( len, XC ) );

            seeds.push_back( it->first );

            //            cout <<BLUE<< it->first << " XC " << XC[0] << " " << XC[1] << " " << XC[2] <<RESET<< endl;
            count++;
        }
    }

    // cout<<"seedsize "<<seeds.size()<<"count  "<<count<<endl;
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
uint TemplateForest<N, Nvalue, M, Mvalue, T>::getTotalSize()
{
    
    uint size = 0;

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        size += ( *it ).size();
    }

    return ( size );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::encodeGeometry()
{
    morton<N> key;
    morton<M> seedkey;
    auto it3 = seeds.begin();

    uint index;
    bool bol;
    real xyz[6];

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).mortonSTLclear();
    }

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        seedkey = ( *it3 );

        auto it4 = geom.find( seedkey );
        if ( it4->second != nullptr )
        {
            for ( uint j = 0; j < it4->second[0]; j++ )
            {
                xyz[0] = it4->second[3 * j + 1];
                xyz[1] = it4->second[3 * j + 2];
                xyz[2] = it4->second[3 * j + 3];

                ( *it ).convertCoordToMorton( xyz, key );
                  //            cout<<"key "<<key<<endl;
            }
        }
        it3 = std::next( it3, 1 );
    }
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::refineEachTree( uint nlevel )
{
    encodeGeometry();

    for ( uint i = 0; i < nlevel; i++ )
    {
        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            ( *it ).pushToRefinelist( i );

            ( *it ).refineRefineList();
        }
    }
}
//
// move function
//

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::updateGeom( T &proc, real *geom_xyz, uint n)
{
 resetGeom();
  auto it2 = trees.begin();

    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        // empty the geom Trees coords ...
        //
        // realloc with size zero is not equivalent to free, as realloc return a pointer
        //    if ( it->second != nullptr )
        {
            //  free(it->second);
            it->second = (real *)realloc( it->second, 0 );
            // set this to nullptr just to make sure
            it->second = nullptr;
        }
        // empty the mortonSTL
        ( *it2 ).clearMortonSTL();
        it2 = std::next( it2, 1 );
    }


    // cout<<"moveGeom "<<n<<endl;
    // move geom calls assign geom
   assignGeom( proc, geom_xyz, n );

}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::moveGeom( T &proc, real *geom_xyz, uint n, real x[3] )
{
    // need to get rid of the encoded geometry first

    /* commented out for multi geometry 
 
    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).mortonSTLclear();
    }

    */
    auto it2 = trees.begin();

    resetGeom();

    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        // empty the geom Trees coords ...
        //
        // realloc with size zero is not equivalent to free, as realloc return a pointer
        //    if ( it->second != nullptr )
        {
            //  free(it->second);
            it->second = (real *)realloc( it->second, 0 );
            // set this to nullptr just to make sure
            it->second = nullptr;
        }
        // empty the mortonSTL
        ( *it2 ).clearMortonSTL();
        it2 = std::next( it2, 1 );
    }

    for ( uint i = 0; i < n; i++ )
    {
        geom_xyz[3 * i] += x[0];
        geom_xyz[3 * i + 1] += x[1];
        geom_xyz[3 * i + 2] += x[2];
//   cout<<" geom  "<<geom_xyz[3*i+1]<<endl;
    }

     cout<<"moveGeom "<<n<<endl;
    // move geom calls assign geom
   assignGeom( proc, geom_xyz, n );
}


template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::moveGeom( T &proc, GeomSTL *GM, uint nSTL, real *x , int * activeList )
{
    // need to get rid of the encoded geometry first

    /* commented out for multi geometry 
 
    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).mortonSTLclear();
    }

    */
    auto it2 = trees.begin();

    resetGeom();

    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        // empty the geom Trees coords ...
        //
        // realloc with size zero is not equivalent to free, as realloc return a pointer
        //    if ( it->second != nullptr )
        {
            //  free(it->second);
            it->second = (real *)realloc( it->second, 0 );
            // set this to nullptr just to make sure
            it->second = nullptr;
        }
        // empty the mortonSTL
        ( *it2 ).clearMortonSTL();
        it2 = std::next( it2, 1 );
    }

for(int j=0;j<nSTL && activeList!=0;j++)
{
   for ( uint i = 0; i < GM[j].geom_nn; i++ )
    {
        GM[j].geom_xyz[3 * i] += x[j*3+0];
        GM[j].geom_xyz[3 * i + 1] += x[j*3+1];
        GM[j].geom_xyz[3 * i + 2] += x[j*3+2];
//   cout<<" geom  "<<geom_xyz[3*i+1]<<endl;
    }
}
//     cout<<"moveGeom "<<n<<endl;
    // move geom calls assign geom
   assignGeom( proc, GM,nSTL,x,activeList  );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::updateGeom( T &proc, GeomSTL *GM, uint nSTL, real *x ,int *activeList)
{
 resetGeom();
  auto it2 = trees.begin();

    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        // empty the geom Trees coords ...
        //
        // realloc with size zero is not equivalent to free, as realloc return a pointer
        //    if ( it->second != nullptr )
        {
            //  free(it->second);
            it->second = (real *)realloc( it->second, 0 );
            // set this to nullptr just to make sure
            it->second = nullptr;
        }
        // empty the mortonSTL
        ( *it2 ).clearMortonSTL();
        it2 = std::next( it2, 1 );
    }


    // cout<<"moveGeom "<<n<<endl;
    // move geom calls assign geom
 

   assignGeom( proc, GM,nSTL,x,activeList  );
}





#if(0)
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::moveGeomDebug( T &proc, real *geom_xyz, uint n, real x[3] )
{
    // need to get rid of the encoded geometry first

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).mortonSTLclear();
    }

    auto it2 = trees.begin();

    resetGeom();

   // cout<<" geom size  "<<geom.size()<<endl;

    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        // empty the geom Trees coords ...
        //
        // realloc with size zero is not equivalent to free, as realloc return a pointer
        //    if ( it->second != nullptr )
        {

            //  free(it->second);
            it->second = (real *)realloc( it->second, 0 );
            // set this to nullptr just to make sure
            it->second = nullptr;

        }
        // empty the mortonSTL
        ( *it2 ).clearMortonSTL();
        it2 = std::next( it2, 1 );
    }

    for ( uint i = 0; i < n; i++ )
    {
        geom_xyz[3 * i] += x[0];
        geom_xyz[3 * i + 1] += x[1];
        geom_xyz[3 * i + 2] += x[2];
//   cout<<" geom  "<<geom_xyz[3*i+1]<<endl;
    }

    // cout<<"moveGeom "<<n<<endl;
    // move geom calls assign geom

   assignGeom( proc, geom_xyz, n );

}
#endif




#if(0)
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::moveGeomDebug( T &proc, real *geom_xyz, uint n, real x[3] )
{
    // need to get rid of the encoded geometry first

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).mortonSTLclear();
    }

    auto it2 = trees.begin();

    resetGeom();

   // cout<<" geom size  "<<geom.size()<<endl;

    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        // empty the geom Trees coords ...
        //
        // realloc with size zero is not equivalent to free, as realloc return a pointer
        //    if ( it->second != nullptr )
        {

            //  free(it->second);
            it->second = (real *)realloc( it->second, 0 );
            // set this to nullptr just to make sure
            it->second = nullptr;

        }
        // empty the mortonSTL
        ( *it2 ).clearMortonSTL();
        it2 = std::next( it2, 1 );
    }

    for ( uint i = 0; i < n; i++ )
    {
        geom_xyz[3 * i] += x[0];
        geom_xyz[3 * i + 1] += x[1];
        geom_xyz[3 * i + 2] += x[2];
//   cout<<" geom  "<<geom_xyz[3*i+1]<<endl;
    }

    // cout<<"moveGeom "<<n<<endl;
    // move geom calls assign geom

   assignGeom( proc, geom_xyz, n );

}
#endif

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::resetGeom(  )
{
    // cout<<"size of geom tree"<<geom.size()<<endl;

    morton<M> key;

    // it is cheaper to just clear this rather than destroy and create teh object
    geom.clearMesh();

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        key = *it;
        geom.insertKey( key );
         //  cout<<RED<<Com.myrank<<" "<<key<<RESET<<endl;
    }

}

// this function along with getListEachTree can be improved by using encoding
// assigns and encodes the geometry
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::assignGeom( T &proc, real *geom_xyz, uint n )
{
    // cout<<"size of geom tree"<<geom.size()<<endl;

    morton<M> key;
    real xyz[6];
    uint count = 0;
/*
    // it is cheaper to just clear this rather than destroy and create teh object
    geom.clearMesh();

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        key = *it;
        geom.insertKey( key );
           cout<<RED<<Com.myrank<<" "<<key<<RESET<<endl;
    }
*/
  cout<<"geom size before "<<geom.size()<<endl;
#if ( 1 )

    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        key = it->first;
        /*
                 if(fixedlevel==0)
                {
                proc.enclosingBox( key, xyz );
                }
                else
                {

                proc.enclosingBoxFixedLevel( key,fixedlevel ,xyz );
                }
        */

        proc.enclosingBox( key, xyz );

          //cout<<" rank "<<Com.myrank<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" "<<xyz[3]<<" "<<xyz[4]<<" "<<xyz[5]<<endl;

        //     cout<<RED<<" "<<Com.myrank<<" "<<key<<RESET<<endl;


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

        //cout<<RED<<" number of points of the geometry inside each seed = "<<count<<RESET<<endl;

        if ( count > 0 )
        {
            it->second = (real *)realloc( it->second, ( 3 * count + 1 ) * sizeof( real ) );
            // cout<<(it->second)<<endl;

            it->second[0] = count;

            // cout<<count<<endl;
             //cout<<"count "<<(it->second)[0]<<endl;

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

//cout<< "start of encoding "<<endl; 


#endif
    double t1 = MPI_Wtime();

    encodeGeometry();

    double t2 = MPI_Wtime();

  //cout<<"Only encodimg "<<t2-t1<<endl;

#if ( 0 )

    uint sum = 0;
    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        if ( it->second != nullptr )
        {
            sum = sum + it->second[0];
        }
    }

  cout << "number of points " << sum << endl;
#endif

    cout<<"geom_size"<<geom.size()<<endl;
    //cout<<" end of encoding " <<endl;
}



template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::assignGeom( T &proc,GeomSTL *GM,  uint nSTL, real *x , int * activeList )
{
    // cout<<"size of geom tree"<<geom.size()<<endl;

    morton<M> key;
    real xyz[6];
    uint count = 0;
/*
    // it is cheaper to just clear this rather than destroy and create teh object
    geom.clearMesh();

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        key = *it;
        geom.insertKey( key );
           cout<<RED<<Com.myrank<<" "<<key<<RESET<<endl;
    }
*/
  cout<<"geom size before "<<geom.size()<<endl;
#if ( 1 )

    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        key = it->first;
        /*
                 if(fixedlevel==0)
                {
                proc.enclosingBox( key, xyz );
                }
                else
                {

                proc.enclosingBoxFixedLevel( key,fixedlevel ,xyz );
                }
        */

        proc.enclosingBox( key, xyz );

          //cout<<" rank "<<Com.myrank<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" "<<xyz[3]<<" "<<xyz[4]<<" "<<xyz[5]<<endl;

        //     cout<<RED<<" "<<Com.myrank<<" "<<key<<RESET<<endl;


        count = 0;

        for(int k=0;k<nSTL && activeList!=0 ;k++)
        { 
        for ( uint j = 0; j < GM[k].geom_nn; j++ )
        {
            if ( GM[k].geom_xyz[3 * j + 0] >= xyz[0] && GM[k].geom_xyz[3 * j + 0] <= xyz[1] )
            {
                if ( GM[k].geom_xyz[3 * j + 1] >= xyz[2] && GM[k].geom_xyz[3 * j + 1] <= xyz[3] )
                {
                    if ( GM[k].geom_xyz[3 * j + 2] >= xyz[4] && GM[k].geom_xyz[3 * j + 2] <= xyz[5] )

                    {
                        count++;
                    }
                }
            }
        }
        }
//        cout<<RED<<" number of points of the geometry inside each seed = "<<count<<RESET<<endl;

        if ( count > 0 )
        {
            it->second = (real *)realloc( it->second, ( 3 * count + 1 ) * sizeof( real ) );
            // cout<<(it->second)<<endl;

            it->second[0] = count;

            // cout<<count<<endl;
             //cout<<"count "<<(it->second)[0]<<endl;

            count = 0;

        for(int k=0;k<nSTL && activeList!=0 ;k++)
        {
            for ( uint j = 0; j < GM[k].geom_nn; j++ )
            {
                if ( GM[k].geom_xyz[3 * j + 0] >= xyz[0] && GM[k].geom_xyz[3 * j + 0] <= xyz[1] )
                {
                    if ( GM[k].geom_xyz[3 * j + 1] >= xyz[2] && GM[k].geom_xyz[3 * j + 1] <= xyz[3] )
                    {
                        if ( GM[k].geom_xyz[3 * j + 2] >= xyz[4] && GM[k].geom_xyz[3 * j + 2] <= xyz[5] )
                        {
                            it->second[3 * count + 1] = GM[k].geom_xyz[3 * j + 0];
                            it->second[3 * count + 2] = GM[k].geom_xyz[3 * j + 1];
                            it->second[3 * count + 3] = GM[k].geom_xyz[3 * j + 2];

                            count++;
                        }
                    }
                }
            }
         }
        }
    }

//cout<< "start of encoding "<<endl; 


#endif
    double t1 = MPI_Wtime();

    encodeGeometry();

    double t2 = MPI_Wtime();

  //cout<<"Only encodimg "<<t2-t1<<endl;

#if (1 )

    uint sum = 0;
    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        if ( it->second != nullptr )
        {
            sum = sum + it->second[0];
        }
    }

  cout << "number of points " << sum << endl;
#endif

    cout<<"geom_size"<<geom.size()<<endl;
    //cout<<" end of encoding " <<endl;
}



#if(0)
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::assignGeom( T &proc, GeomSTL *GM,int nSTL )
{
    // cout<<"size of geom tree"<<geom.size()<<endl;

    morton<M> key;
    real xyz[6];
    uint count = 0;
/*
    // it is cheaper to just clear this rather than destroy and create teh object
    geom.clearMesh();

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        key = *it;
        geom.insertKey( key );
           cout<<RED<<Com.myrank<<" "<<key<<RESET<<endl;
    }
*/
//  cout<<"geom size "<<geom.size()<<endl;
#if ( 1 )

    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        key = it->first;
        proc.enclosingBox( key, xyz );

          //cout<<" rank "<<Com.myrank<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" "<<xyz[3]<<" "<<xyz[4]<<" "<<xyz[5]<<endl;

        //     cout<<RED<<" "<<Com.myrank<<" "<<key<<RESET<<endl;

        count = 0;
      
        for(int k=0;k<nSTL;k++)
       {
        for ( uint j = 0; j < GM[k].geom_nn; j++ )
        {
            if ( GM[k].geom_xyz[3 * j + 0] >= xyz[0] && GM[k].geom_xyz[3 * j + 0] <= xyz[1] )
            {
                if ( GM[k].geom_xyz[3 * j + 1] >= xyz[2] && GM[k].geom_xyz[3 * j + 1] <= xyz[3] )
                {
                    if ( GM[k].geom_xyz[3 * j + 2] >= xyz[4] && GM[k].geom_xyz[3 * j + 2] <= xyz[5] )
                    {
                        count++;
                    }
                }
            }
        }
       }
        //cout<<RED<<" number of points of the geometry inside each seed = "<<count<<RESET<<endl;

        if ( count > 0 )
        {
            it->second = (real *)realloc( it->second, ( 3 * count + 1 ) * sizeof( real ) );
            // cout<<(it->second)<<endl;

            it->second[0] = count;

            // cout<<count<<endl;
             //cout<<"count "<<(it->second)[0]<<endl;

         count = 0;

        for(int k=0;k<nSTL;k++)
        {
            for ( uint j = 0; j < GM[k].geom_nn; j++ )
            {
                if ( GM[k].geom_xyz[3 * j + 0] >= xyz[0] && GM[k].geom_xyz[3 * j + 0] <= xyz[1] )
                {
                    if ( GM[k].geom_xyz[3 * j + 1] >= xyz[2] && GM[k].geom_xyz[3 * j + 1] <= xyz[3] )
                    {
                        if ( GM[k].geom_xyz[3 * j + 2] >= xyz[4] && GM[k].geom_xyz[3 * j + 2] <= xyz[5] )
                        {
                            it->second[3 * count + 1] = GM[k].geom_xyz[3 * j + 0];
                            it->second[3 * count + 2] = GM[k].geom_xyz[3 * j + 1];
                            it->second[3 * count + 3] = GM[k].geom_xyz[3 * j + 2];

                            count++;
                        }
                    }
                }
            }
        } 
       }
    }

//cout<< "start of encoding "<<endl; 


#endif
    double t1 = MPI_Wtime();

    encodeGeometry();

    double t2 = MPI_Wtime();

  //cout<<"Only encodimg "<<t2-t1<<endl;

#if ( 1 )

    uint sum = 0;
    for ( auto it = geom.begin(); it != geom.end(); it++ )
    {
        if ( it->second != nullptr )
        {
            sum = sum + it->second[0];
        }
    }

  //cout << "number of points " << sum << endl;
#endif

    //cout<<"geom_size"<<geom.size()<<endl;
    //cout<<" end of encoding " <<endl;
}

#endif

//===========================================================
//
// Generating the list that will lead to 4:1 balance
//
//===========================================================

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::getListEachTree()
{
    morton<N> key;

    morton<M> seedkey;
    auto it3 = seeds.begin();

    uint index;
    bool bol;
    real xyz[6];

    it3 = seeds.begin();
    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        seedkey = ( *it3 );

        for ( auto it2 = ( *it ).begin(); it2 != ( *it ).end(); it2++ )
        {
            key = it2->first;
            bol = 0;

            auto it4 = geom.find( seedkey );
            ( *it ).enclosingBox( key, xyz );

            if ( it4->second != nullptr )
            {
                for ( uint j = 0; j < it4->second[0]; j++ )
                {
                    if ( it4->second[3 * j + 1] >= xyz[0] && it4->second[3 * j + 1] <= xyz[1] )
                    {
                        if ( it4->second[3 * j + 2] >= xyz[2] && it4->second[3 * j + 2] <= xyz[3] )
                        {
                            if ( it4->second[3 * j + 3] >= xyz[4] && it4->second[3 * j + 3] <= xyz[5] )
                            {
                                bol = 1;
                                break;
                            }
                        }
                    }
                }
            }

            if ( bol )
            {
                ( *it ).addToList( key );
            }
        }

        // increment the pointer
        it3 = std::next( it3, 1 );
    }
}

//
//
//              Extend the tagged list to account for the local 4:1 balance
//
//

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
bool TemplateForest<N, Nvalue, M, Mvalue, T>::isInSeed( morton<M> &key, uint *counter )
{
    bool bol = false;
    uint count = 0;
    *counter = 0;

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        if ( ( *it ) == key )
        {
            bol = true;
            break;
        }
        else
        {
            count++;
        }
    }

    *counter = count;

    return ( bol );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::flipAll( morton<N> &key, uint *mylevel, uint *direction )
{
    for ( uint i = 0; i < ( *mylevel ); i++ )
    {
        key.flip( N - 3 * i - ( *direction ) - 1 );
    }
}

//
//              extends the list for interior elements
//

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::getDirections( morton<N + M> &key, uint combinedlevel, vector<uint> &directions )
{
    directions.clear();
    uint mylevel;
    bool bol1, bol2;
    morton<N + M> kt1, kt2;
    // level(key, &mylevel );
    bol2 = false;
    // excluding when the level is zero
    // since this is an exception as the other elements are siblings not nonlocal neighbor
    // cout<<"key= "<<key<<"level "<<mylevel<<endl;
    // need to chec for all three directions for every element
    if ( combinedlevel != 1 )
    {
        for ( uint j = 0; j < 3; j++ )
        {
            bol1 = key[( N + M ) - 1 - j];
            // cout<<"j direction"<<bol1<<endl;
            kt1 = 0;
            // if bol1=false basically the generated auxillary code will be zero
            if ( bol1 == true )
            {
                for ( uint k = 0; k < combinedlevel; k++ )
                {
                    kt1[( N + M ) - 1 - 3 * k - j] = bol1;
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

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::recoverAllZeroSingularity( morton<N + M> &key, const uint &combinedlevel )
{
    if ( key[0] == true && combinedlevel != 0 )
    {
        key.flip( 0 );
        key.flip( N + M - 3 * ( combinedlevel - 1 ) - 1 );
    }
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::combinedLevel( const morton<N + M> &key, uint *level )
{
    *level = ( N + M ) / 3;
    uint rem = ( N + M ) % 3;
    uint iend = ( N + M ) / 3;

    /*! to prevent unnecesary bit operation, the morton code is placed from starting from  left hand side*/

    morton<N + M> kt;

    kt = key;
    if ( kt[0] == 1 )
    {
        kt.flip( 0 );
    }

    for ( uint i = 0; i < iend; i++ )
    {
        if ( kt[3 * i + rem] == false && kt[3 * i + rem + 1] == false && kt[3 * i + rem + 2] == false )
        {
            *level = *level - 1;
        }
        else
        {
            break;
        }
    }
}

// check this one

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::findSeedLevelForRcvdMessage( const morton<N + M> &key, uint *mylevel, Tree<M, Mvalue> &proc )
{
    morton<M> kt = 0;

    //    bool bol = false;
    // my level is 1 at the very least, 0 is invalid
    *mylevel = 0;

    // maxseedlevel=1;
    // cout << BLUE << "max seed level " << maxseedlevel << "rank " << Com.myrank << RESET << endl;
    //  cout << key << endl;
    for ( uint i = 0; i < 3 * maxseedlevel; i++ )
    {
        kt[M - i - 1] = key[M + N - 1 - i];
        // cout<<" "<<kt<<endl;
    }

    //    cout << kt << endl;
    uint checklevel;

    for ( uint i = maxseedlevel; i > 0; i-- )
    {
        proc.level( kt, &checklevel );

        if ( std::find( seeds.begin(), seeds.end(), kt ) != seeds.end() && i == checklevel )
        {
            ( *mylevel ) = i;
            //            cout<<RED<<"found"<<RESET<<endl;
            break;
        }

        else
        {
            kt[M - 3 * ( i - 1 ) - 1] = 0;
            kt[M - 3 * ( i - 1 ) - 2] = 0;
            kt[M - 3 * ( i - 1 ) - 3] = 0;
        }
    }

    // cout << GREEN << "seedlevel " << *mylevel << RESET << endl;

    /*
            if ( std::find( seeds.begin(), seeds.end(), kt ) == seeds.end() )
           {
               throw std::runtime_error( RED "Unable to find the seed corresponding to the recieved global key " RESET );
          }
    */
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest
<N, Nvalue, M, Mvalue, T>::findSeedLevelForRcvdMessage( const morton<N + M> &key, uint *mylevel, FullTree<M, Mvalue> &proc )
{

    *mylevel = proc.getLevel();
    // cout<<" ==============  "<<(*mylevel)<<endl;
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest
<N, Nvalue, M, Mvalue, T>::constructSeedKeyForRcvdMessage( const morton<N + M> &key, const uint &seedlevel, morton<M> &seedkey )
{
    seedkey = 0;
    for ( uint i = 0; i < 3 * seedlevel; i++ )
    {
        seedkey[M - i - 1] = key[M + N - i - 1];
    }
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest
<N, Nvalue, M, Mvalue, T>::constructElementKeyForRcvdMessage( const morton<N + M> &key, const uint &seedlevel, morton<N> &elementkey )
{
    // cout << "seedlevel inside" << seedlevel << endl;
    for ( uint i = 0; i < N; i++ )
    {
        elementkey[N - i - 1] = key[M + N - 3 * ( seedlevel ) - i - 1];
    }
    // do not forget the tag for the singularity
    elementkey[0] = key[0];
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::removeAllZeroSingularity( morton<N + M> &key, const uint &combinedlevel )
{
    // will not be able to tell what the combined level is if it ends in 000 at the finest level after flipping the key
    // to prevent this and make this data available and yet not add to communication, the x-direction
    // bit related to the finest level is flipped and then the far right bit is also flipped to tag this change
    // on the reciever side if that bit is true we flip back this bit before using it

    uint NM = M + N;

    if ( key[NM - 3 * ( combinedlevel - 1 ) - 1] == false && key[NM - 3 * ( combinedlevel - 1 ) - 2] == false
         && key[NM - 3 * ( combinedlevel - 1 ) - 3] == false )
    {
        key.flip( N + M - 3 * ( combinedlevel - 1 ) - 1 );
        key[0] = true;
    }
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::getMaxSeedsLevel( T &proc )
{
    uint seedlevel = 0;
    // initializa the maxseed level to 1
    maxseedlevel = 0;

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        proc.level( *it, &seedlevel );
        if ( seedlevel > maxseedlevel )
        {
            maxseedlevel = seedlevel;
        }
    }
    // cout<<GREEN "MAXSEEDLEVEL "<<maxseedlevel<<RESET<<endl;
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
uint TemplateForest<N, Nvalue, M, Mvalue, T>::getMaxSeedLevel() 
{
 
return maxseedlevel;
}


template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::setMaxProcLevel( const uint refinelevel )
{
    maxProcLevel = refinelevel;
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::findFlipLevel( morton<N + M> key, uint *mylevel, uint *changedirectionlevel, uint *direction )
{
    bool bol;

    const uint NM = N + M;

    bol = key[NM - 3 * ( ( *mylevel ) - 1 ) - ( *direction ) - 1];

    // cout<<"index is = "<<bit-3*((*mylevel)-1)-(*direction)-1<<endl;

    // cout<<"bol = "<<bol<<endl;

    *changedirectionlevel = 0;

    for ( uint i = ( *mylevel ) - 1; i > 0; i-- )
    {
        // cout<<"i="<<i <<key[bit-3*(i-1)-(*direction)-1]<<endl;

        if ( key[NM - 3 * ( i - 1 ) - ( *direction ) - 1] != bol )
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

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::flipForNbr( morton<N + M> &key, uint *mylevel, uint *changedirectionlevel, uint *direction )
{
    //    cout << "combinedlevel, changelevel and direction = " << *mylevel << "\t" << *changedirectionlevel << "\t" << *direction << endl;

    uint NM = N + M;
    // cout<<"NM= "<<NM<<endl;
    // if changedirectionlevel==0 that node is a boundary node
    if ( *changedirectionlevel > 0 )
    {
        for ( uint i = ( *changedirectionlevel ); i <= ( *mylevel ); i++ )
        {
            //      cout<<"index "<<(NM-3*(i-1)-(*direction)-1)<<endl;

            // cout<<"direction= "<<*direction<<" "<<"mylevel "<<*mylevel<<"changedir " <<*changedirectionlevel<<endl;
            //         cout<<"index "<<NM-3*i+3-(*direction)-1<<endl;
            // cout<<"i ="<<i<<endl;

            ( key ).flip( NM - 3 * ( i - 1 ) - ( *direction ) - 1 );

            // cout<<"in here"<<(*key)<<endl;
        }
    }

    // if node is a boundary element. nonlocal nbr does not exist
    /*
    else if(*changedirectionlevel==0 )
    {
    for(uint i=0;i<(*mylevel);i++)
    {
    key.flip(NM-3*i-(*direction)-1);
    }

    //         ( key ).flip( NM -  ( *direction ) - 1 );
    }
    */
    // cout<<"exiting flip NBR ..."<<endl;
}

// watch out, this might oveflow for a huge mesh

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::getTotalMeshSize() /*!< creates distributed (acalable) graph for communication */
{

//    unsigned long long meshSize;
    unsigned long long forestsize = getTotalSize();
    int comsize;

    //  cout<<" "<<forestsize<<endl;

    MPI_Reduce( &forestsize, &meshSize, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );

    MPI_Comm_size( MPI_COMM_WORLD, &comsize );

    if ( Com.myrank == 0 )
    {

        std::string filename = "meshSize";
        filename.append( to_string( Com.myrank ) );
        ofstream myfile;
        myfile.open( filename );

        /*
                cout << GREEN "===========================================================\n" << RESET << endl;

                cout << MAGENTA "\ttotal_mesh_across_all_procs_size = " << meshSize << RESET << endl;

                cout << GREEN "\n===========================================================\n" << RESET << endl;
        */

        myfile << "\ttotal_number_of_processors = " << comsize << endl;

        myfile << "\ttotal_mesh_across_all_procs_size = " << std::setprecision( 12 ) << (double)meshSize << endl;

        myfile.close();
    }
}

#if ( 0 )
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest
<N, Nvalue, M, Mvalue, T>::getNbrSeedLevel( morton<N + M> &combinedkey, uint topologylevel, uint *nbrseedlevel, Tree<M, Mvalue> &proc )
{
    morton<M> kt = 0;
    uint l1, mylevel;
    bool bol = false;
    const uint NM = N + M;

    //  cout << RED << topologylevel << RESET << endl;

    // this proc is a balanced tree therefore it either has the flipped key or it has one level lower or one level higher than that key,
    // for higher level the code would exist

    for ( uint j = 0; j < 3 * ( topologylevel ); j++ )
    {
        kt[M - j - 1] = combinedkey[NM - j - 1];
        // proc.level(kt,&mylevel);
    }

    if ( proc.find( kt ) != proc.end() )
    {
        proc.level( kt, &mylevel );
    }
    /*
    if(mylevel>topologylevel)
    {
    throw std::runtime_error(RED "level is higher than current level" RESET);
    }
    */
    else
    {
        kt = 0;
        for ( uint j = 0; j < 3 * ( topologylevel - 1 ); j++ )
        {
            kt[M - j - 1] = combinedkey[NM - j - 1];
            // proc.level(kt,&mylevel);
        }
        proc.level( kt, &mylevel );
    }

    *nbrseedlevel = mylevel;

    // return(bol);
}
#endif

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest
<N, Nvalue, M, Mvalue, T>::getNbrSeedLevel( morton<N + M> &combinedkey, uint topologylevel, uint *nbrseedlevel, Tree<M, Mvalue> &proc )
{
    morton<M> kt = 0;
    uint l1, mylevel;
    bool bol = false;
    const uint NM = N + M;

    //  cout << RED << topologylevel << RESET << endl;

    // this proc is a balanced tree therefore it either has the flipped key or it has one level lower or one level higher than that key,
    // for higher level the code would exist

    for ( uint j = 0; j < 3 * ( topologylevel ); j++ )
    {
        kt[M - j - 1] = combinedkey[NM - j - 1];
        // proc.level(kt,&mylevel);
    }

    uint count = 1;

    while ( proc.find( kt ) == proc.end() )
    {
        kt = 0;
        for ( uint j = 0; j < 3 * ( topologylevel - count ); j++ )
        {
            kt[M - j - 1] = combinedkey[NM - j - 1];
        }

        count++;
    }

    if ( proc.find( kt ) == proc.end() )
    {
    
        cout << " Exit Code : " << ReblAmrGetErrorEnum( NO_SEED ) << endl;
        exit(1);  
//throw std::runtime_error( "seed not found in the proc" );
    }
    else
    {

        proc.level( kt, &mylevel );
    }

    *nbrseedlevel = mylevel;

    // return(bol);
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest
<N, Nvalue, M, Mvalue, T>::getNbrSeedLevel( morton<N + M> &combinedkey, uint topologylevel, uint *nbrseedlevel, FullTree<M, Mvalue> &proc )
{

    *nbrseedlevel = topologylevel;
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
uint TemplateForest<N, Nvalue, M, Mvalue, T>::forestsize()
{
    uint size = 0;
    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        size += ( *it ).size();
    }

    return ( size );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::getElemNbrs( Tree<M, Mvalue> &proc, const morton<M> key, bitvector<M> &nbr )
{
    nbr.clear();
    morton<M> nbrkey = key, sibkey = key, kt, kb;
    uint mylevel, siblevel, nbrlevel, changedirectionlevel, direction;
    // assume that the gut has siblings of same level
    // each direction
    // siblings frist
    // for directon loop

    int forprint = 37;

    vector<uint> directions;

    proc.level( key, &mylevel );
    // cout << mylevel << endl;
    morton<M> kn[4];
    // construct siblings in 3 direction

    if ( !proc.isInMeshList( key ) )
    {
        throw std::runtime_error( "seed does not exist" );
    }

    // cout<<RED<<key<<" "<<RESET<<endl;
    for ( uint direction = 0; direction < 3; direction++ )
    {
        sibkey = key;
        sibkey.flip( M - 3 * ( mylevel - 1 ) - 1 - direction );
        proc.level( sibkey, &siblevel );
        kt = sibkey;

#if ( DEBUG )
        if ( Com.myrank == forprint )
        {
            auto itd = proc.find( kt );

            //  cout << key << " " << direction << " " << BLUE << sibkey << " " << mylevel << " " << siblevel << RESET << "  " <<
            // itd->second[0]
            //      << endl;
            //cout << " sibling " << sibkey << endl;
        }
#endif
        if ( siblevel > mylevel )
        {
            // base one,
            // sibkey.flip(N-3*mylevel-1-direction);
            // now add three more
            // cout<<RED<<"hereeeeeeeeeeeeeeeeeeeeeeeee "<<endl;

            proc.constructHigherLevelNbrs( kt, mylevel, direction, kn );
            for ( uint l = 0; l < 4; l++ )
            {
                nbr.push_back( kn[l] );
#if ( DEBUG )
                if ( Com.myrank == forprint )
                {
                 //   cout << key << " "
                         << "sibling has higher level " << kn[l] << endl;
                }
#endif
            }
        }
        else
        {
            nbr.push_back( kt );
        }
    }

    // need to watch out for the first generation as these
    // only have siblings
    // directions for level one does not mke sense as
    // neighbors are siblings remove this exception
    uint modlevel;
#if ( 1 )

    // look for the possible nonlocal neighbors in 3 directions excluding the boundary element by if(changedirectionlevel!=0)
    // this is for nonlocal and therefore constructing higher order neighbor is different, dont copy paste from above
    for ( uint j = 0; j < 3; j++ )
    {
        // flip combinedkey  to construct global code
        // E
        // uint j=0;
        // may wan to exclude the boundary elements from being checked

        //   direction = directions.at( j );
        direction = j;
        // cout<<"direction "<<direction<<endl;
        proc.findFlipLevel( key, &mylevel, &changedirectionlevel, &direction );
        //             cout<<"CHANGE LEVEL "<<GREEN<<mylevel<<" "<<mylevel<<RESET<<" "<<changedirectionlevel<<endl;
        //        cout << BLUE << key << RESET << endl;
        if ( changedirectionlevel != 0 )
        {
            nbrkey = key;

// cout<<"key "<<key<<endl;

// exclude the boundary elemnts that chanagedirectionlevel=0 for them
// taken care of in tree function

#if ( DEBUG )
            cout << "before " << nbrkey << endl;
#endif

            proc.flipForNbr( &nbrkey, &mylevel, &changedirectionlevel, &direction );

#if ( DEBUG )
            cout << "after  " << nbrkey << endl;
#endif

            proc.level( nbrkey, &nbrlevel );
#if ( DEBUG )

            if ( Com.myrank == forprint )
            {
                auto itd2 = proc.find( nbrkey );

                // cout << GREEN << "key " << key << " nbrkey " << nbrkey <<" mylevel "<<mylevel<<" nbrlevel " <<nbrlevel<<" direction
                // "<<direction<<RESET <<" "<< itd2->second[0] <<endl;
                if ( itd2 != proc.end() )
                {
                   // cout << GREEN << "key " << key << " nbrkey " << nbrkey << " mylevel " << mylevel << " nbrlevel " << nbrlevel
                     //    << " direction " << direction << RESET << " " << itd2->second[0] << endl;
                }

                //  cout<<" this dir "<<directions.at(j)<<endl;
                //   cout<<" this dir "<<j<<endl;
            }
#endif

            kt = nbrkey;

            //      cout << GREEN << "key " << key << "nbrkey " << nbrkey << " mylevel " << mylevel << " nbrlevel " << nbrlevel
            //           << " direction size " << direction << RESET << endl;
            //???????????????????????????????????????????????????????
            // This is very important, is not so important in 4:1 balance
            // this section checks to see if neighbor has a lower level
            //
            auto it5 = proc.find( kt );

            if ( it5 == proc.end() )
            //           if(!proc.isInMeshList(nbrkey))
            {
                //        cout << "look for lower level" << endl;
                kt[M - 3 * ( mylevel - 1 ) - 1] = 0;
                kt[M - 3 * ( mylevel - 1 ) - 2] = 0;
                kt[M - 3 * ( mylevel - 1 ) - 3] = 0;
                nbrlevel = nbrlevel - 1;

#if ( 0 )
                if ( Com.myrank == forprint )
                {
                    auto itd4 = proc.find( kt );

                    // cout << GREEN << "key " << key << " nbrkey " << nbrkey <<" mylevel "<<mylevel<<" nbrlevel " <<nbrlevel<<" direction
                    // "<<direction<<RESET <<" "<< itd2->second[0] <<endl;
                    if ( itd4 != proc.end() )
                    {
                        cout << BLUE << "key " << key << " nbrkey " << kt << " mylevel " << mylevel << " nbrlevel " << nbrlevel
                             << " direction " << direction << RESET << " " << itd4->second[0] << endl;
                    }
                }
#endif
            }
            //  cout << RED << "key " << key << " nbrkey " << kt << RESET << endl;

            if ( nbrlevel > mylevel )
            {
                // base one,
                // sibkey.flip(N-3*mylevel-1-direction);
                // now add three more
                proc.constructNonlocalHigherLevelNbrs( nbrkey, mylevel, direction, kn );

                for ( uint l = 0; l < 4; l++ )
                {
                    nbr.push_back( kn[l] );
                }
#if ( DEBUG )
                if ( Com.myrank == forprint )
                {

                    for ( uint l = 0; l < 4; l++ )
                    {
                        auto itd3 = proc.find( kn[l] );
                        cout << "nonlocal higher order sibling " << kn[l] << endl;
                        // cout << GREEN << "key " << key << " nbrkey " << nbrkey <<" mylevel "<<mylevel<<" nbrlevel " <<nbrlevel<<"
                        // direction
                        // "<<direction<<RESET <<" "<< itd2->second[0] <<endl;
                        if ( itd3 != proc.end() )
                        {
                            cout << GREEN << "key " << key << " nbrkey " << kn[l] << " mylevel " << mylevel << " nbrlevel " << nbrlevel
                                 << " direction " << direction << RESET << " " << itd3->second[0] << endl;
                        }
                    }
                }
#endif
            }

            else
            {
                it5 = proc.find( kt );

                if ( it5 == proc.end() )
                {
                    throw std::runtime_error( "element not in topology" );
                }
                nbr.push_back( kt );
            }
        }
    }

#endif
    /*
        if ( Com.myrank == 2 )
        {
            cout << RED << "myrank " << Com.myrank << "key " << key << "size " << nbr.size() << RESET << endl;
            cout << "====================================================" << endl;
            cout << key << endl;

            for ( uint i = 0; i < nbr.size(); i++ )
            {
                cout << nbr.at( i ) << endl;
            }
        }
    */
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::comPatternConstruct( Tree<M, Mvalue> &proc )
{
    morton<M> seedkey;
    uint mylevel;
    uint myrank = Com.myrank;
     destination.clear();


    bitvector<M> nbr;

    morton<M> key;

    auto it = seeds.begin();

    if ( Com.comsize == 1 )
    {
        return;
    }

    for ( uint i = 0; i < seeds.size(); i++ )
    {
        seedkey = ( *it );
//   cout << YELLOW << seedkey << endl;

#if ( DEBUG )
        if ( Com.myrank == 37 )
        {

            cout << "============== seed==========\n " << seedkey << endl;
            /*
               for ( uint j = 0; j < nbr.size(); j++ )
                    {
                        cout <<  nbr.at( j ) <<  endl;
                    }
            */
            cout << "==========================" << endl;
        }
#endif
        getElemNbrs( proc, seedkey, nbr );
#if ( 1 )
        // cout<<nbr.at(0)<<endl;

        for ( uint j = 0; j < nbr.size(); j++ )
        {
            key = nbr.at( j );

            auto it1 = proc.find( nbr.at( j ) );

            if ( it1 == proc.end() )
            {
                cout << BLUE << key << "not  found" << RESET << endl;
                throw std::runtime_error( RED "Element not found in comPatternConstruct" RESET );
                // std::cout<<RED<<it1->first<<RESET<<endl;
            }

            // cout<<nbr.at(0)<<endl;
            // cout<<key<<end;
            // flip to get the sibling
            if ( it1->second[0] != Com.myrank && std::find( destination.begin(), destination.end(), it1->second[0] ) == destination.end() )
            {
                destination.push_back( it1->second[0] );
            }
#endif
        }

        it = std::next( it, 1 );
    }

    message = new vector<bitset<M + N>>[destination.size()];
    request = new MPI_Request[destination.size()];
    request1 = new MPI_Request[destination.size()];

#if ( 0 )

    if ( Com.myrank == 17 )
    {
        cout << "=========================================" << endl;
        //cout << GREEN << Com.myrank << RESET << endl;

        for ( uint i = 0; i < destination.size(); i++ )
        {
            //cout << RED << "rank " << destination.at( i ) << RESET << endl;
        }
        cout << "=========================================" << endl;
    }
#endif
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::comPatternConstruct( FullTree<M, Mvalue> &proc, vector<uint> &Nbrs )
{

    uint co = 0;

    for ( uint i = 0; i < Nbrs.size(); i++ )
    {
        destination.push_back( Nbrs.at( i ) );
    }

    message = new vector<bitset<M + N>>[destination.size()];
    request = new MPI_Request[destination.size()];
    request1 = new MPI_Request[destination.size()];

    // cout<<"nbrs_size"<<Nbrs.size()<<endl;
    /*
        for ( auto it = proc.begin(); it != proc.end(); it++ )
        {
            if ( it->second[0] == Com.myrank )
            {
                continue;
            }
            else
            {
                destination.push_back( Nbrs.at( co ) );
                // cout<<my_rank<<" "<<Nbrs.at(co)<<endl;
                co++;
            }
        }
    */
    /*
    for(uint i=0; i< destination.size();i++)
    {
    cout<<Com.myrank<<" "<<destination.at(i)<<endl;
    }
    */
}

#if ( METHOD == 0 )
//
// boolean used to communicate tree
//
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::fourToOneBalance( T &proc )
{
    const uint NM = N + M;
    morton<NM> combinedkey = 0, ktcom;
    morton<N> key, kt1, nbrkey;
    morton<M> seedkey, seednbrkey = 0, kt, kt3;
    auto it3 = seeds.begin();
    uint index;
    bool bol;
    real xyz[6];
    bitvector<N> boundaryElem;
    uint mylevel, changedirectionlevel, direction;
    uint topologylevel, nbrseedlevel, nbrlevel, complevel, nbrcomplevel, nbrtopologylevel, combinedlevel, localcombinedlevel;
    uint counter;
    vector<uint> directions;
    it3 = seeds.begin();
    vector<bitset<M + N>> message[destination.size()];
    vector<uint> dest;
    vector<uint> source;
    uint idx;
    int size;
    int flag = 1;
    morton<M + N> rcvkey;
    uint seedlevel;
    morton<N> elementkey;
    uint idex2;
    uint elementlevel, recvmessagecombinedlevel;
    //    morton<M + N> *       sendbuff[destination.size()];
    //    morton<M + N> *       recvbuff[destination.size()];
    vector<uint> start;
    vector<uint> end;
    uint istart, iend;
    uint counter3, counter2, counter5;
    MPI_Request request[destination.size()], request1[destination.size()], request0;
    MPI_Status status;
    bool sb = 0, rb = 0;
    morton<N + M> tempkey;

#if ( 0 )
    std::string filename = "rank";
    filename.append( to_string( Com.myrank ) );
    ofstream myfile;
    myfile.open( filename );

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        ////  cout << RED << "prod id " << Com.myrank << "  seed " << ( *it ) << RESET << endl;
        //   myfile << "prod id " << Com.myrank << "  " << con++ << "  seed " << ( *it ) << endl;
    }
// myfile << " topology size " << proc.size() << endl;

#endif
    bool *sendbuff[destination.size()];
    uint con = 0;
    bool *recvbuff[destination.size()];

    // extract the boundary elements for each tree from the tagged list
    // get the ectent that we need to refine in the initia step
    // we know the number of trees therefore, initialize it here
    // myfile << "Writing this to a file.\n";
    /*
     * previuos version
        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            start.push_back( 0 );
            end.push_back( ( *it ).refineListSize() );
        }
    */
    // the main while loop

    //  cout << "size start=" << start.size() << endl;
    //   cout << "destination ize " << destination.size() << endl;
    //   cout << "size end=" << end.size() << endl;

    bool loopbool = true;
    uint bcstart, bcend;
    bool innerloop = true;

    while ( loopbool )
    {
        // start out by setting loopbool=0
        loopbool = false;
        counter2 = 0;

        // enforce 4: balance for each tree and extract the boundary elements from the list of elements to be refined

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            //            ( *it ).fourToOneP( start.at( counter2 ), end.at( counter2 ) );

            ( *it ).fourToOne();
            ( *it ).refinelistReset();
            counter2++;
        }

        // separates the list between local boundary and non-local boundary
        it3 = seeds.begin();

#if ( 1 )
        counter5 = 0;

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            seedkey = ( *it3 );
            proc.level( seedkey, &topologylevel );

            ( *it ).refinelistReset();

            innerloop = true;
            while ( innerloop )
            {
                innerloop = false;

                // for ( uint i = start.at( counter5 ); i < end.at( counter5 ); i++ )
                for ( auto i = ( *it ).Rbegin(); i != ( *it ).Rend(); i++ )
                {
                    combinedkey = 0;

                    for ( uint j = 0; j < 3 * topologylevel; j++ )
                    {
                        combinedkey[NM - j - 1] = seedkey[M - j - 1];
                    }

                    auto p = ( *it ).readRefineList( i );
                    // if the tag = 1 then we have already investigated this guy, continue
                    if ( p.second == 1 )
                    {
                        continue;
                    }
                    //                    cout<<RED<<p.first<<endl;

                    key = p.first;

                    ( *it ).level( key, &mylevel );

                    combinedlevel = topologylevel + mylevel;

                    for ( uint j = 0; j < 3 * mylevel; j++ )
                    {
                        combinedkey[NM - 3 * ( topologylevel ) - j - 1] = key[N - j - 1];
                    }

                    ktcom = combinedkey;

                    for ( uint direction = 0; direction < 3; direction++ )
                    {
                        if ( !( *it ).isBoundary( key, direction ) )
                        {
                            continue;
                        }
                        //  cout << " direction " << direction << endl;
                        //   cout << " combined key " << combinedkey << endl;
                        findFlipLevel( combinedkey, &combinedlevel, &changedirectionlevel, &direction );

                        if ( changedirectionlevel != 0 )
                        {
                            flipForNbr( combinedkey, &combinedlevel, &changedirectionlevel, &direction );
#if ( 1 )
                            //     cout << " combined key " << combinedkey << " combined level " << combinedlevel << endl;
                            getNbrSeedLevel( combinedkey, maxProclevel, &nbrseedlevel, proc );

                            seednbrkey = 0;
                            for ( uint k = 0; k < 3 * nbrseedlevel; k++ )
                            {
                                seednbrkey[M - 1 - k] = combinedkey[NM - k - 1];
                                // cout<<RED<<combinedkey[NM-k-1]<<RESET<<endl;
                            }

                            //   cout << "seednbrkey " << BLUE << seednbrkey << RESET << endl;
                            // myfile << "  seed key  "<<seednbrkey <<  endl;
                            nbrkey = 0;

                            for ( uint k = 0; k < N; k++ )
                            {
                                nbrkey[N - k - 1] = combinedkey[NM - 3 * ( nbrseedlevel ) - 1 - k];
                            }
                            //  cout << " nbrkey " << nbrkey << endl;

                            if ( isInSeed( seednbrkey, &counter ) )
                            {
                                auto it2 = std::next( trees.begin(), counter );
                                auto it4 = ( *it2 ).find( nbrkey );

                                // seedlevel
                                // ( proc ).level( seednbrkey, &nbrlevel );
                                // initial guess on nbr level as constructed
                                // this is a singularity for level zero, it might be negative for first level
                                // I remove this mannually
                                if ( combinedlevel >= nbrseedlevel )
                                {
                                    nbrlevel = combinedlevel - nbrseedlevel;
                                }
                                else
                                {
                                    nbrlevel = 1;
                                }
                                //      cout<<RED<<nbrlevel<<RESET<<endl;

                                // now find the real levels

                                if ( it4 != ( *it2 ).end() )
                                {
                                    ( *it2 ).level( it4->first, &nbrlevel );
                                    //         cout << YELLOW << nbrlevel << RESET << endl;
                                }
                                else
                                {
                                    // this condition imples level is lower no need for modification and search, I commented it out
                                    //                                nbrlevel                     = mylevel-1;
                                    // ( *it2 ).level( it4->first, &nbrlevel );
                                    nbrlevel = nbrlevel - 1;
                                    nbrkey[N - 3 * nbrlevel - 1] = 0;
                                    nbrkey[N - 3 * nbrlevel - 2] = 0;
                                    nbrkey[N - 3 * nbrlevel - 3] = 0;

                                    //            cout <<RED<< N - 3 * nbrlevel - 1 << " " << mylevel << RESET<<endl;
                                    //                                      auto it4=(*it2).find(nbrkey);
                                    //                                      (*it).level(it4->first,&nbrlevel);
                                }

                                // cout << "nbrkey " << RED << nbrkey << "  count " << counter << RESET << endl;
                                // cout<<"inside"<<endl;
                                if ( ( *it2 ).find( nbrkey ) == ( *it2 ).end() )
                                {
                                    cout << GREEN << nbrkey << "direction " << direction << "level " << nbrlevel << RESET << endl;
                                    throw std::runtime_error( "error in finding key" RESET );
                                }

                                nbrcomplevel = nbrlevel + nbrseedlevel;

                                //                             cout<<RED<<nbrlevel<<" "<<nbrseedlevel<<RESET<<endl;
                                //                             cout<<YELLOW<<nbrcomplevel<<" " <<combinedlevel<<RESET<<endl;
                                // if not in the list add to the list, dont forget to modify here
                                // watch out this element might have already been tagged

                                if ( nbrcomplevel < combinedlevel && ( *it2 ).isInRefineList( nbrkey ) == false )
                                {
                                    ( *it2 ).addToList( nbrkey );
                                    innerloop = true;
                                    // here is where boolwhile is affected
                                    loopbool = true;
                                    //       cout << RED "======================================" << endl;
                                    //       cout << "added to list" << nbrkey << endl;
                                    //       cout << "======================================" RESET << endl;
                                }
                                // cout<<seednbrkey<<RED<<seednbrkey<<RESET<<endl;
                            }

                            else
                            {
                                // sorts for communication

                                // prepare for communication, pack all the elements with the destination, sender is   my_rank

                                // find seednbrkey from global data
                                auto it5 = proc.find( seednbrkey );
                                //                        dest.push_back( it5->second[0] );
                                auto it6 = destination.begin();
                                if ( it5 != proc.end() )
                                {
                                    it6 = find( destination.begin(), destination.end(), it5->second[0] );
                                }

                                // combined key should be avoided
                                idx = it6 - destination.begin();
                                // cout << "---------------------->>> index" << RED << idx << RESET << endl;
                                // before pushing back, need to be careful about the key with all zeros,
                                // to accommodate this, the first bit (far right) is flipped and a sibling is sent
                                removeAllZeroSingularity( combinedkey, combinedlevel );
                                // cout << GREEN << " proc_id " << Com.myrank << "seednbrkey  " << seednbrkey << "combinedkey  " <<
                                // combinedkey
                                //     << "combined level " << combinedlevel << endl;
                                message[idx].push_back( combinedkey );

                                // myfile<< " combined key "<< combinedkey <<endl;
                                //                        cout << message.at( 0 ) << RESET << endl;
                            }
                        }
                        combinedkey = ktcom;
                    }
                    // here change the int value for ith element of the it-th tree

                    ( *it ).flipRefineElemTag( i );
                }
            }
#endif

            counter5++;
            it3 = std::next( it3, 1 );
        }

#endif

// need a communication pattern need to send zero size message if communication is not required  need to add local balance criteria  need to
// loop over and do this again in while loop
// get size of message for each proc send the required information to the neighboring processes

#if ( 1 )
        for ( uint i = 0; i < destination.size(); i++ )
        {
            //    sendbuff[i] = new morton<M + N>[ message[i].size() ];
            /*
                        for ( uint j = 0; j < message[i].size(); j++ )
                        {
                            sendbuff[i][j] = message[i].at( j );
                        }
            */
            // sending with boolean is supposed to be more portable
            sendbuff[i] = new bool[message[i].size() * ( M + N )];
            // counter6=0;

            for ( uint j = 0; j < message[i].size(); j++ )
            {
                tempkey = message[i].at( j );

                //     myfile<<" message to be sent "<<tempkey<<endl;
                for ( uint k = 0; k < M + N; k++ )
                {
                    if ( tempkey[k] == true )
                    {
                        sendbuff[i][j * ( M + N ) + k] = true;
                    }
                    else
                    {
                        sendbuff[i][j * ( M + N ) + k] = false;
                    }
                }
            }

            //            cout << "================================== " << endl;
            //            cout << RESET "messagesize " << message[i].size() << " destination " << destination.at( i ) << RESET << endl;
            //            cout << "================================== " << endl;
            // send messages to destinations and tag each meaasge with self rank (myrank)
            //
            MPI_Isend( sendbuff[i], message[i].size() * sizeof( bool ) * ( M + N ), MPI_BYTE, destination.at( i ), Com.myrank,
                       MPI_COMM_WORLD, &request1[i] );

            //           myfile<<" proc "<<Com.myrank<<" sends "<<message[i].size()*sizeof(bool)*(M+N)<<" to "<<destination.at(i)<<endl;
            //  MPI_Isend( sendbuff[i], message[i].size() * sizeof( bool )*(M+N), MPI_BYTE, destination.at( i ),0  , MPI_COMM_WORLD,
            // &request[i] );
            //  MPI_Isend( sendbuff[i], message[i].size() * sizeof( morton<N + M> ), MPI_BYTE, destination.at( i ), 0, MPI_COMM_WORLD,
            // &request[i] );
        }

#if ( 0 )
        //   enforce 4:1 balance internally to give enough time for the message to be able to be probed
        //   This is to overlap comm and computation

        counter2 = 0;

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            istart = start.at( counter2 );
            iend = end.at( counter2 );
            ( *it ).refineRefineList( istart, iend );
        }
#endif
        // even with overlapping communication and computation it is not garanteed that this part
        // will use enought time for the message to arrive, since this list might be empty
        // however we can use blocking after this section, that we know not all of the time is spent waiting

        for ( uint i = 0; i < destination.size(); i++ )
        {
            //             MPI_Iprobe(destination.at(i), destination.at(i), MPI_COMM_WORLD,&flag, &status);
            // this is a symmetric comm pattern therefore any
            //
            MPI_Probe( destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &status );

            //   MPI_Probe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
            // while(flag!=0)
            {
                // MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&flag, &status);
                //  cout << BLUE << flag << RESET << endl;
                // if ( flag == 1 )
                {
                    MPI_Get_count( &status, MPI_BYTE, &size );

                    // recvbuff[i] = new morton<M + N>[ size ];
                    // notice we revcieved the message by byte so need to specify,  to find the number of elements simly divide it by sizeof
                    // bool
                    recvbuff[i] = new bool[size / sizeof( bool )];
                    //           cout << "******************* " << size << endl;

                    MPI_Irecv( recvbuff[i], size, MPI_BYTE, destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &request[i] );

                    MPI_Wait( &request[i], &status );
                    //                    myfile<<" proc "<<Com.myrank<<" recvs "<<size<<"from"<<destination.at(i)<<endl;

                    /*  MPI_Irecv( recvbuff[i], size, MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &request0 );

                      MPI_Wait( &request0, &status );
                    */

                    // if ( size != 0 )
                    {
                        //                cout << "rank " << Com.myrank << "recvbuff ??????????????? " << recvbuff[0][0] << endl;

                        // if that violates the balance condition add it to the list in the appropriate location

                        for ( uint j = 0; j < size / sizeof( bool ) / ( M + N ); j++ )
                        {
                            rcvkey = 0;

                            for ( uint k = 0; k < N + M; k++ )
                            {
                                if ( recvbuff[i][j * ( N + M ) + k] == true )
                                {
                                    rcvkey.flip( k );
                                }
                            }

                            //                   cout << "p_id " << Com.myrank << " " << rcvkey << endl;

                            //  myfile << " p_id  " << Com.myrank << " " << rcvkey << endl;

                            // restore if modified to handle singularity
                            // restoration needsto be done before this func

                            combinedLevel( rcvkey, &recvmessagecombinedlevel );

                            recoverAllZeroSingularity( rcvkey, recvmessagecombinedlevel );

                            findSeedLevelForRcvdMessage( rcvkey, &seedlevel, proc );

                            //                           cout << rcvkey << endl;
                            //                            cout << "seedlevel " << seedlevel << endl;

                            constructSeedKeyForRcvdMessage( rcvkey, seedlevel, seedkey );
                            //                           cout << "seedkey " << seedkey << endl;

                            /*      	    uint check;
                                               proc.level(seedkey,&check);
                                               cout<<RED<<check<<" "<<seedlevel<<RESET<<endl;
                                               if(seedlevel!=check)
                                              {
                                               throw std::runtime_error("levels calculated for seed inconsistent");
                                               }
                              */
                            constructElementKeyForRcvdMessage( rcvkey, seedlevel, elementkey );
                            //            cout << "element key " << elementkey << endl;

                            // cout << RED "rcvkey level " << recvmessagecombinedlevel << RESET << endl;

                            auto it7 = std::find( seeds.begin(), seeds.end(), seedkey );

                            if ( it7 == seeds.end() )
                            {
                                cout << RED << "myrank " << Com.myrank << " " << seedkey << RESET << endl;

                                throw std::runtime_error( RED "seed not found in Balance Comm" RESET );
                            }

                            idex2 = std::distance( seeds.begin(), it7 );

                            auto it9 = std::next( trees.begin(), idex2 );

                            if ( recvmessagecombinedlevel >= seedlevel )
                            {
                                elementlevel = recvmessagecombinedlevel - seedlevel;
                            }
                            else
                            {
                                elementlevel = 1;
                            }
                            // cout << BLUE << elementlevel << RESET << endl;

                            // myfile <<" recv message combined level " <<  recvmessagecombinedlevel <<" seed level "<< seedlevel << endl;
                            // myfile <<" element level " << elementlevel<< " elementkey  "<< elementkey  << endl;

                            // recover the singularity
                            // error here, I flip the code such that we can get the level
                            // myfile << rcvkey << " " << recvmessagecombinedlevel << endl;
                            // myfile << elementkey << endl;
                            //  myfile << seedkey << " " << seedlevel << endl;
                            //  myfile << elementlevel << endl;
                            ( *it9 ).level( elementkey, &elementlevel );
                            //  myfile << elementlevel << endl;

                            //                            recoverAllZeroSingularity( elementkey, elementlevel );

                            // get the real level
                            if ( ( *it9 ).isInMeshList( elementkey ) == false )
                            {
                                //   cout << "inside mesh" << elementkey << endl;
                                elementkey[N - 3 * ( elementlevel - 1 ) - 1] = 0;

                                elementkey[N - 3 * ( elementlevel - 1 ) - 2] = 0;

                                elementkey[N - 3 * ( elementlevel - 1 ) - 3] = 0;
                            }
                            //     myfile<<elementkey<<endl;
                            ( *it9 ).level( elementkey, &elementlevel );

                            localcombinedlevel = elementlevel + seedlevel;

                            // add if the level is lower

                            if ( localcombinedlevel < recvmessagecombinedlevel && ( *it9 ).isInRefineList( elementkey ) == false )
                            {
                                //                cout << BLUE << localcombinedlevel << "            " << recvmessagecombinedlevel << RESET
                                // << endl;

                                //               cout << "elementkey" << elementkey << endl;
                                ( *it9 ).addToList( elementkey );
                                loopbool = true;
                                /*                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;

                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                */
                            }

                            // see if this elemnt exits, if yes, we are ok, if not add to tree referenced by pointer *it9
                        }
                    }
                }
            }
        }

#endif

        if ( loopbool == 1 )
        {
            sb = 1;
        }
        else
        {
            sb = 0;
        }

        // this value is boolan,use MPI_BYTE, since bool size is implementation dependent, use sizeof(bool)

        MPI_Iallreduce( &sb, &rb, sizeof( bool ), MPI_BYTE, MPI_MAX, MPI_COMM_WORLD, &request0 );

        // do some computation here
        //
        counter3 = 0;

        // for debug for now do the refinement separately
        /*
                for ( auto it = trees.begin(); it != trees.end(); it++ )
                {
                    start.at( counter3 ) = end.at( counter3 );
                    end.at( counter3 ) = ( *it ).refineListSize();
                    counter3++;
                }
        */
        MPI_Wait( &request0, MPI_STATUS_IGNORE );

        for ( uint j = 0; j < destination.size(); j++ )
        {
            MPI_Wait( &request1[j], &status );
            delete[] sendbuff[j];
            message[j].clear();
            message[j].shrink_to_fit();
            delete[] recvbuff[j];
        }
        //
        // need to use the recv buffer
        //
        if ( rb == 1 )
        {
            loopbool = true;
            for ( auto it = trees.begin(); it != trees.end(); it++ )
            {
                ( *it ).refinelistReset();
            }
        }
        if ( loopbool == 1 )
        {
            //           cout << RED << " loopbool " << loopbool << " rb " << rb << " sb " << sb << RESET << endl;
        }
    }

    //    myfile.close();

    //   start.clear();
    //   end.clear();
    // assigne recieved elements to corresponding lists if needed
}

//
//
//  use character to communicate tree
//
//

#elif( METHOD == 1 )
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::fourToOneBalance( T &proc )
{
    // cout << "character solving" << endl;
    const uint NM = N + M;
    morton<NM> combinedkey = 0, ktcom;
    morton<N> key, kt1, nbrkey;
    morton<M> seedkey, seednbrkey = 0, kt, kt3;
    auto it3 = seeds.begin();
    uint index;
    bool bol;
    real xyz[6];
    bitvector<N> boundaryElem;
    uint mylevel, changedirectionlevel, direction;
    uint topologylevel, nbrseedlevel, nbrlevel, complevel, nbrcomplevel, nbrtopologylevel, combinedlevel, localcombinedlevel;
    uint counter;
    vector<uint> directions;
    it3 = seeds.begin();
    vector<bitset<M + N>> message[destination.size()];
    vector<uint> dest;
    vector<uint> source;
    uint idx;
    int size;
    int flag = 1;
    morton<M + N> rcvkey;
    uint seedlevel;
    morton<N> elementkey;
    uint idex2;
    uint elementlevel, recvmessagecombinedlevel;
    //    morton<M + N> *       sendbuff[destination.size()];
    //    morton<M + N> *       recvbuff[destination.size()];
    vector<uint> start;
    vector<uint> end;
    uint istart, iend;
    uint counter3, counter2, counter5;
    MPI_Request request[destination.size()], request1[destination.size()], request0;
    MPI_Status status;
    bool sb = 0, rb = 0;
    morton<N + M> tempkey;
#if ( 0 )
    std::string filename = "rank";
    filename.append( to_string( Com.myrank ) );
    ofstream myfile;
    myfile.open( filename );

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        ////  cout << RED << "prod id " << Com.myrank << "  seed " << ( *it ) << RESET << endl;
        myfile << "prod id " << Com.myrank << "  " << con++ << "  seed " << ( *it ) << endl;
    }
    myfile << " topology size " << proc.size() << endl;
#endif

    uint con = 0;
    char *sendbuff[destination.size()];
    char *recvbuff[destination.size()];

    // extract the boundary elements for each tree from the tagged list
    // get the ectent that we need to refine in the initia step
    // we know the number of trees therefore, initialize it here
    // myfile << "Writing this to a file.\n";
    /*
     * previuos version
        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            start.push_back( 0 );
            end.push_back( ( *it ).refineListSize() );
        }
    */
    // the main while loop

    //  cout << "size start=" << start.size() << endl;
    //   cout << "destination ize " << destination.size() << endl;
    //   cout << "size end=" << end.size() << endl;

    bool loopbool = true;
    uint bcstart, bcend;
    bool innerloop = true;

    while ( loopbool )
    {
        // start out by setting loopbool=0
        loopbool = false;
        counter2 = 0;

        // enforce 4: balance for each tree and extract the boundary elements from the list of elements to be refined

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            //            ( *it ).fourToOneP( start.at( counter2 ), end.at( counter2 ) );

            ( *it ).fourToOne();
            ( *it ).refinelistReset();
            counter2++;
        }

        // separates the list between local boundary and non-local boundary
        it3 = seeds.begin();

#if ( 1 )
        counter5 = 0;

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            seedkey = ( *it3 );
            proc.level( seedkey, &topologylevel );

            ( *it ).refinelistReset();

            innerloop = true;
            while ( innerloop )
            {
                innerloop = false;

                // for ( uint i = start.at( counter5 ); i < end.at( counter5 ); i++ )
                for ( auto i = ( *it ).Rbegin(); i != ( *it ).Rend(); i++ )
                {
                    combinedkey = 0;

                    for ( uint j = 0; j < 3 * topologylevel; j++ )
                    {
                        combinedkey[NM - j - 1] = seedkey[M - j - 1];
                    }

                    auto p = ( *it ).readRefineList( i );
                    // if the tag = 1 then we have already investigated this guy, continue
                    if ( p.second == 1 )
                    {
                        continue;
                    }
                    //                    cout<<RED<<p.first<<endl;

                    key = p.first;

                    ( *it ).level( key, &mylevel );

                    combinedlevel = topologylevel + mylevel;

                    for ( uint j = 0; j < 3 * mylevel; j++ )
                    {
                        combinedkey[NM - 3 * ( topologylevel ) - j - 1] = key[N - j - 1];
                    }

                    ktcom = combinedkey;

                    for ( uint direction = 0; direction < 3; direction++ )
                    {
                        if ( !( *it ).isBoundary( key, direction ) )
                        {
                            continue;
                        }
                        //  cout << " direction " << direction << endl;
                        //   cout << " combined key " << combinedkey << endl;
                        findFlipLevel( combinedkey, &combinedlevel, &changedirectionlevel, &direction );

                        if ( changedirectionlevel != 0 )
                        {
                            flipForNbr( combinedkey, &combinedlevel, &changedirectionlevel, &direction );
#if ( 1 )
                            //     cout << " combined key " << combinedkey << " combined level " << combinedlevel << endl;
                            getNbrSeedLevel( combinedkey, maxProcLevel, &nbrseedlevel, proc );

                            seednbrkey = 0;
                            for ( uint k = 0; k < 3 * nbrseedlevel; k++ )
                            {
                                seednbrkey[M - 1 - k] = combinedkey[NM - k - 1];
                                // cout<<RED<<combinedkey[NM-k-1]<<RESET<<endl;
                            }

                            //   cout << "seednbrkey " << BLUE << seednbrkey << RESET << endl;
                            // myfile << "  seed key  "<<seednbrkey <<  endl;
                            nbrkey = 0;

                            for ( uint k = 0; k < N; k++ )
                            {
                                nbrkey[N - k - 1] = combinedkey[NM - 3 * ( nbrseedlevel ) - 1 - k];
                            }
                            //  cout << " nbrkey " << nbrkey << endl;

                            if ( isInSeed( seednbrkey, &counter ) )
                            {
                                auto it2 = std::next( trees.begin(), counter );
                                auto it4 = ( *it2 ).find( nbrkey );

                                // seedlevel
                                // ( proc ).level( seednbrkey, &nbrlevel );
                                // initial guess on nbr level as constructed
                                // this is a singularity for level zero, it might be negative for first level
                                // I remove this mannually
                                if ( combinedlevel >= nbrseedlevel )
                                {
                                    nbrlevel = combinedlevel - nbrseedlevel;
                                }
                                else
                                {
                                    nbrlevel = 1;
                                }
                                //      cout<<RED<<nbrlevel<<RESET<<endl;

                                // now find the real levels

                                if ( it4 != ( *it2 ).end() )
                                {
                                    ( *it2 ).level( it4->first, &nbrlevel );
                                    //         cout << YELLOW << nbrlevel << RESET << endl;
                                }
                                else
                                {
                                    // this condition imples level is lower no need for modification and search, I commented it out
                                    //                                nbrlevel                     = mylevel-1;
                                    // ( *it2 ).level( it4->first, &nbrlevel );
                                    nbrlevel = nbrlevel - 1;
                                    nbrkey[N - 3 * nbrlevel - 1] = 0;
                                    nbrkey[N - 3 * nbrlevel - 2] = 0;
                                    nbrkey[N - 3 * nbrlevel - 3] = 0;

                                    //            cout <<RED<< N - 3 * nbrlevel - 1 << " " << mylevel << RESET<<endl;
                                    //                                      auto it4=(*it2).find(nbrkey);
                                    //                                      (*it).level(it4->first,&nbrlevel);
                                }

                                // cout << "nbrkey " << RED << nbrkey << "  count " << counter << RESET << endl;
                                // cout<<"inside"<<endl;
                                if ( ( *it2 ).find( nbrkey ) == ( *it2 ).end() )
                                {
                                    cout << GREEN << nbrkey << "direction " << direction << "level " << nbrlevel << RESET << endl;
                                    throw std::runtime_error( "error in finding key" RESET );
                                }

                                nbrcomplevel = nbrlevel + nbrseedlevel;

                                //                             cout<<RED<<nbrlevel<<" "<<nbrseedlevel<<RESET<<endl;
                                //                             cout<<YELLOW<<nbrcomplevel<<" " <<combinedlevel<<RESET<<endl;
                                // if not in the list add to the list, dont forget to modify here
                                // watch out this element might have already been tagged

                                if ( nbrcomplevel < combinedlevel && ( *it2 ).isInRefineList( nbrkey ) == false )
                                {
                                    ( *it2 ).addToList( nbrkey );
                                    innerloop = true;
                                    // here is where boolwhile is affected
                                    loopbool = true;
                                    //       cout << RED "======================================" << endl;
                                    //       cout << "added to list" << nbrkey << endl;
                                    //       cout << "======================================" RESET << endl;
                                }
                                // cout<<seednbrkey<<RED<<seednbrkey<<RESET<<endl;
                            }

                            else
                            {
                                // sorts for communication

                                // prepare for communication, pack all the elements with the destination, sender is   my_rank

                                // find seednbrkey from global data
                                auto it5 = proc.find( seednbrkey );
                                //                        dest.push_back( it5->second[0] );
                                auto it6 = destination.begin();
                                if ( it5 != proc.end() )
                                {
                                    it6 = find( destination.begin(), destination.end(), it5->second[0] );
                                }

                                // combined key should be avoided
                                idx = it6 - destination.begin();
                                // cout << "---------------------->>> index" << RED << idx << RESET << endl;
                                // before pushing back, need to be careful about the key with all zeros,
                                // to accommodate this, the first bit (far right) is flipped and a sibling is sent
                                removeAllZeroSingularity( combinedkey, combinedlevel );
                                // cout << GREEN << " proc_id " << Com.myrank << "seednbrkey  " << seednbrkey << "combinedkey  " <<
                                // combinedkey
                                //     << "combined level " << combinedlevel << endl;
                                message[idx].push_back( combinedkey );

                                // myfile<< " combined key "<< combinedkey <<endl;
                                //                        cout << message.at( 0 ) << RESET << endl;
                            }
                        }
                        combinedkey = ktcom;
                    }
                    // here change the int value for ith element of the it-th tree

                    ( *it ).flipRefineElemTag( i );
                }
            }
#endif

            counter5++;
            it3 = std::next( it3, 1 );
        }

#endif

// need a communication pattern need to send zero size message if communication is not required  need to add local balance criteria  need to
// loop over and do this again in while loop
// get size of message for each proc send the required information to the neighboring processes

#if ( 1 )
        for ( uint i = 0; i < destination.size(); i++ )
        {
            //    sendbuff[i] = new morton<M + N>[ message[i].size() ];
            /*
                        for ( uint j = 0; j < message[i].size(); j++ )
                        {
                            sendbuff[i][j] = message[i].at( j );
                        }
            */
            // sending with boolean is supposed to be more portable
            sendbuff[i] = new char[message[i].size() * ( M + N )];
            // counter6=0;

            for ( uint j = 0; j < message[i].size(); j++ )
            {
                tempkey = message[i].at( j );

                //     myfile<<" message to be sent "<<tempkey<<endl;
                for ( uint k = 0; k < M + N; k++ )
                {
                    if ( tempkey[k] == true )
                    {
                        sendbuff[i][j * ( M + N ) + k] = '1';
                    }
                    else
                    {
                        sendbuff[i][j * ( M + N ) + k] = '0';
                    }
                }
            }

            //            cout << "================================== " << endl;
            //            cout << RESET "messagesize " << message[i].size() << " destination " << destination.at( i ) << RESET << endl;
            //            cout << "================================== " << endl;
            // send messages to destinations and tag each meaasge with self rank (myrank)
            //
            MPI_Isend( sendbuff[i], message[i].size() * sizeof( char ) * ( M + N ), MPI_BYTE, destination.at( i ), Com.myrank,
                       MPI_COMM_WORLD, &request1[i] );

            //           myfile<<" proc "<<Com.myrank<<" sends "<<message[i].size()*sizeof(bool)*(M+N)<<" to "<<destination.at(i)<<endl;
            //  MPI_Isend( sendbuff[i], message[i].size() * sizeof( bool )*(M+N), MPI_BYTE, destination.at( i ),0  , MPI_COMM_WORLD,
            // &request[i] );
            //  MPI_Isend( sendbuff[i], message[i].size() * sizeof( morton<N + M> ), MPI_BYTE, destination.at( i ), 0, MPI_COMM_WORLD,
            // &request[i] );
        }

#if ( 0 )
        //   enforce 4:1 balance internally to give enough time for the message to be able to be probed
        //   This is to overlap comm and computation

        counter2 = 0;

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            istart = start.at( counter2 );
            iend = end.at( counter2 );
            ( *it ).refineRefineList( istart, iend );
        }
#endif
        // even with overlapping communication and computation it is not garanteed that this part
        // will use enought time for the message to arrive, since this list might be empty
        // however we can use blocking after this section, that we know not all of the time is spent waiting

        for ( uint i = 0; i < destination.size(); i++ )
        {
            //             MPI_Iprobe(destination.at(i), destination.at(i), MPI_COMM_WORLD,&flag, &status);
            // this is a symmetric comm pattern therefore any
            //
            MPI_Probe( destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &status );

            //   MPI_Probe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
            // while(flag!=0)
            {
                // MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&flag, &status);
                //  cout << BLUE << flag << RESET << endl;
                // if ( flag == 1 )
                {
                    MPI_Get_count( &status, MPI_BYTE, &size );

                    // recvbuff[i] = new morton<M + N>[ size ];
                    // notice we revcieved the message by byte so need to specify,  to find the number of elements simly divide it by sizeof
                    // bool
                    recvbuff[i] = new char[size / sizeof( char )];
                    //           cout << "******************* " << size << endl;

                    MPI_Irecv( recvbuff[i], size, MPI_BYTE, destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &request[i] );

                    MPI_Wait( &request[i], &status );
                    //                    myfile<<" proc "<<Com.myrank<<" recvs "<<size<<"from"<<destination.at(i)<<endl;

                    /*  MPI_Irecv( recvbuff[i], size, MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &request0 );

                      MPI_Wait( &request0, &status );
                    */

                    // if ( size != 0 )
                    {
                        //                cout << "rank " << Com.myrank << "recvbuff ??????????????? " << recvbuff[0][0] << endl;

                        // if that violates the balance condition add it to the list in the appropriate location

                        for ( uint j = 0; j < size / sizeof( char ) / ( M + N ); j++ )
                        {
                            rcvkey = 0;

                            for ( uint k = 0; k < N + M; k++ )
                            {
                                if ( recvbuff[i][j * ( N + M ) + k] == '1' )
                                {
                                    rcvkey.flip( k );
                                }
                            }

                            //                   cout << "p_id " << Com.myrank << " " << rcvkey << endl;

                            //  myfile << " p_id  " << Com.myrank << " " << rcvkey << endl;

                            // restore if modified to handle singularity
                            // restoration needsto be done before this func

                            combinedLevel( rcvkey, &recvmessagecombinedlevel );

                            recoverAllZeroSingularity( rcvkey, recvmessagecombinedlevel );

                            findSeedLevelForRcvdMessage( rcvkey, &seedlevel, proc );

                            //                           cout << rcvkey << endl;
                            //                            cout << "seedlevel " << seedlevel << endl;

                            constructSeedKeyForRcvdMessage( rcvkey, seedlevel, seedkey );
                            //                           cout << "seedkey " << seedkey << endl;

                            /*      	    uint check;
                                               proc.level(seedkey,&check);
                                               cout<<RED<<check<<" "<<seedlevel<<RESET<<endl;
                                               if(seedlevel!=check)
                                              {
                                               throw std::runtime_error("levels calculated for seed inconsistent");
                                               }
                              */
                            constructElementKeyForRcvdMessage( rcvkey, seedlevel, elementkey );
                            //            cout << "element key " << elementkey << endl;

                            // cout << RED "rcvkey level " << recvmessagecombinedlevel << RESET << endl;

                            auto it7 = std::find( seeds.begin(), seeds.end(), seedkey );

                            if ( it7 == seeds.end() )
                            {
                                //cout << RED << "myrank " << Com.myrank << " " << seedkey << RESET << endl;

                                throw std::runtime_error( RED "seed not found in Balance Comm" RESET );
                            }

                            idex2 = std::distance( seeds.begin(), it7 );

                            auto it9 = std::next( trees.begin(), idex2 );

                            if ( recvmessagecombinedlevel >= seedlevel )
                            {
                                elementlevel = recvmessagecombinedlevel - seedlevel;
                            }
                            else
                            {
                                elementlevel = 1;
                            }
                            // cout << BLUE << elementlevel << RESET << endl;

                            // myfile <<" recv message combined level " <<  recvmessagecombinedlevel <<" seed level "<< seedlevel << endl;
                            // myfile <<" element level " << elementlevel<< " elementkey  "<< elementkey  << endl;

                            // recover the singularity
                            // error here, I flip the code such that we can get the level
                            // myfile << rcvkey << " " << recvmessagecombinedlevel << endl;
                            // myfile << elementkey << endl;
                            // myfile << seedkey << " " << seedlevel << endl;
                            // myfile << elementlevel << endl;
                            ( *it9 ).level( elementkey, &elementlevel );
                            // myfile << elementlevel << endl;

                            //                            recoverAllZeroSingularity( elementkey, elementlevel );

                            // get the real level
                            if ( ( *it9 ).isInMeshList( elementkey ) == false )
                            {
                                //   cout << "inside mesh" << elementkey << endl;
                                elementkey[N - 3 * ( elementlevel - 1 ) - 1] = 0;

                                elementkey[N - 3 * ( elementlevel - 1 ) - 2] = 0;

                                elementkey[N - 3 * ( elementlevel - 1 ) - 3] = 0;
                            }
                            //     myfile<<elementkey<<endl;
                            ( *it9 ).level( elementkey, &elementlevel );

                            localcombinedlevel = elementlevel + seedlevel;

                            // add if the level is lower

                            if ( localcombinedlevel < recvmessagecombinedlevel && ( *it9 ).isInRefineList( elementkey ) == false )
                            {
                                //                cout << BLUE << localcombinedlevel << "            " << recvmessagecombinedlevel << RESET
                                // << endl;

                                //               cout << "elementkey" << elementkey << endl;
                                ( *it9 ).addToList( elementkey );
                                loopbool = true;
                                /*                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;

                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                */
                            }

                            // see if this elemnt exits, if yes, we are ok, if not add to tree referenced by pointer *it9
                        }
                    }
                }
            }
        }

#endif

        if ( loopbool == 1 )
        {
            sb = 1;
        }
        else
        {
            sb = 0;
        }

        // this value is boolan,use MPI_BYTE, since bool size is implementation dependent, use sizeof(bool)

        MPI_Iallreduce( &sb, &rb, sizeof( bool ), MPI_BYTE, MPI_MAX, MPI_COMM_WORLD, &request0 );

        // do some computation here
        //

        counter3 = 0;

        // for debug for now do the refinement separately
        /*
                for ( auto it = trees.begin(); it != trees.end(); it++ )
                {
                    start.at( counter3 ) = end.at( counter3 );
                    end.at( counter3 ) = ( *it ).refineListSize();
                    counter3++;
                }
        */

        MPI_Wait( &request0, MPI_STATUS_IGNORE );

        for ( uint j = 0; j < destination.size(); j++ )
        {
            MPI_Wait( &request1[j], &status );
            delete[] sendbuff[j];
            message[j].clear();
            message[j].shrink_to_fit();
            delete[] recvbuff[j];
        }
        //
        // need to use the recv buffer
        //
        if ( rb == 1 )
        {
            loopbool = true;
            for ( auto it = trees.begin(); it != trees.end(); it++ )
            {
                ( *it ).refinelistReset();
            }
        }
        if ( loopbool == 1 )
        {
            cout << RED << " loopbool " << loopbool << " rb " << rb << " sb " << sb << RESET << endl;
        }
    }

    //    myfile.close();

    //   start.clear();
    //   end.clear();
    // assigne recieved elements to corresponding lists if needed
}
#elif( METHOD == 2 )
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::fourToOneBalance( T &proc )
{
    // cout << "character solving" << endl;
    const uint NM = N + M;
    morton<NM> combinedkey = 0, ktcom;
    morton<N> key, kt1, nbrkey;
    morton<M> seedkey, seednbrkey = 0, kt, kt3;
    auto it3 = seeds.begin();
    uint index;
    bool bol;
    real xyz[6];
    bitvector<N> boundaryElem;
    uint mylevel, changedirectionlevel, direction;
    uint topologylevel, nbrseedlevel, nbrlevel, complevel, nbrcomplevel, nbrtopologylevel, combinedlevel, localcombinedlevel;
    uint counter;
    vector<uint> directions;
    it3 = seeds.begin();
    vector<bitset<M + N>> message[destination.size()];
    vector<uint> dest;
    vector<uint> source;
    uint idx;
    int size;
    int flag = 1;
    morton<M + N> rcvkey;
    uint seedlevel;
    morton<N> elementkey;
    uint idex2;
    uint elementlevel, recvmessagecombinedlevel;
    //    morton<M + N> *       sendbuff[destination.size()];
    //    morton<M + N> *       recvbuff[destination.size()];
    vector<uint> start;
    vector<uint> end;
    uint istart, iend;
    uint counter3, counter2, counter5;
    MPI_Request request[destination.size()], request1[destination.size()], request0;
    MPI_Status status;
    int sb = 0, rb = 0;
    morton<N + M> tempkey;
#if ( 0 )
    std::string filename = "rank";
    filename.append( to_string( Com.myrank ) );
    ofstream myfile;
    myfile.open( filename );

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        ////  cout << RED << "prod id " << Com.myrank << "  seed " << ( *it ) << RESET << endl;
        myfile << "prod id " << Com.myrank << "  " << con++ << "  seed " << ( *it ) << endl;
    }
    myfile << " topology size " << proc.size() << endl;
#endif

    uint con = 0;
    char *sendbuff[destination.size()];
    char *recvbuff[destination.size()];

    // extract the boundary elements for each tree from the tagged list
    // get the ectent that we need to refine in the initia step
    // we know the number of trees therefore, initialize it here
    // myfile << "Writing this to a file.\n";
    /*
     * previuos version
        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            start.push_back( 0 );
            end.push_back( ( *it ).refineListSize() );
        }
    */
    // the main while loop

    //  cout << "size start=" << start.size() << endl;
    //   cout << "destination ize " << destination.size() << endl;
    //   cout << "size end=" << end.size() << endl;

    bool loopbool = true;
    uint bcstart, bcend;
    bool innerloop = true;

    while ( loopbool )
    {
        // start out by setting loopbool=0
        loopbool = false;
        counter2 = 0;

        // enforce 4: balance for each tree and extract the boundary elements from the list of elements to be refined

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            //            ( *it ).fourToOneP( start.at( counter2 ), end.at( counter2 ) );

            ( *it ).fourToOne();
            ( *it ).refinelistReset();
            counter2++;
        }

        // separates the list between local boundary and non-local boundary
        it3 = seeds.begin();

#if ( 1 )
        counter5 = 0;

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            seedkey = ( *it3 );
            proc.level( seedkey, &topologylevel );
            // cout<<GREEN<<"Top lovel "<<topologylevel<<"seedkey "<<seedkey<<RESET<<endl;

            ( *it ).refinelistReset();

            innerloop = true;
            while ( innerloop )
            {
                innerloop = false;

                // for ( uint i = start.at( counter5 ); i < end.at( counter5 ); i++ )
                for ( auto i = ( *it ).Rbegin(); i != ( *it ).Rend(); i++ )
                {
                    combinedkey = 0;

                    for ( uint j = 0; j < 3 * topologylevel; j++ )
                    {
                        combinedkey[NM - j - 1] = seedkey[M - j - 1];
                    }

                    auto p = ( *it ).readRefineList( i );
                    // if the tag = 1 then we have already investigated this guy, continue
                    if ( p.second == 1 )
                    {
                        continue;
                    }
                    //                    cout<<RED<<p.first<<endl;

                    key = p.first;

                    ( *it ).level( key, &mylevel );

                    combinedlevel = topologylevel + mylevel;

                    for ( uint j = 0; j < 3 * mylevel; j++ )
                    {
                        combinedkey[NM - 3 * ( topologylevel ) - j - 1] = key[N - j - 1];
                    }

                    ktcom = combinedkey;

                    for ( uint direction = 0; direction < 3; direction++ )
                    {
                        if ( !( *it ).isBoundary( key, direction ) )
                        {
                            continue;
                        }
                        //  cout << " direction " << direction << endl;
                        //   cout << " combined key " << combinedkey << endl;
                        findFlipLevel( combinedkey, &combinedlevel, &changedirectionlevel, &direction );

                        if ( changedirectionlevel != 0 )
                        {
                            flipForNbr( combinedkey, &combinedlevel, &changedirectionlevel, &direction );
#if ( 1 )
                            //     cout << " combined key " << combinedkey << " combined level " << combinedlevel << endl;
                            getNbrSeedLevel( combinedkey, maxProcLevel, &nbrseedlevel, proc );

                            seednbrkey = 0;
                            for ( uint k = 0; k < 3 * nbrseedlevel; k++ )
                            {
                                seednbrkey[M - 1 - k] = combinedkey[NM - k - 1];
                                // cout<<RED<<combinedkey[NM-k-1]<<RESET<<endl;
                            }

                            //   cout << "seednbrkey " << BLUE << seednbrkey << RESET << endl;
                            // myfile << "  seed key  "<<seednbrkey <<  endl;
                            nbrkey = 0;

                            for ( uint k = 0; k < N; k++ )
                            {
                                nbrkey[N - k - 1] = combinedkey[NM - 3 * ( nbrseedlevel ) - 1 - k];
                            }
                            //  cout << " nbrkey " << nbrkey << endl;

                            if ( isInSeed( seednbrkey, &counter ) )
                            {
                                auto it2 = std::next( trees.begin(), counter );
                                auto it4 = ( *it2 ).find( nbrkey );

                                // seedlevel
                                // ( proc ).level( seednbrkey, &nbrlevel );
                                // initial guess on nbr level as constructed
                                // this is a singularity for level zero, it might be negative for first level
                                // I remove this mannually
                                if ( combinedlevel >= nbrseedlevel )
                                {
                                    nbrlevel = combinedlevel - nbrseedlevel;
                                }
                                else
                                {
                                    nbrlevel = 1;
                                }
                                //      cout<<RED<<nbrlevel<<RESET<<endl;

                                // now find the real levels

                                if ( it4 != ( *it2 ).end() )
                                {
                                    ( *it2 ).level( it4->first, &nbrlevel );
                                    //         cout << YELLOW << nbrlevel << RESET << endl;
                                }
                                else
                                {
                                    // this condition imples level is lower no need for modification and search, I commented it out
                                    //                                nbrlevel                     = mylevel-1;
                                    // ( *it2 ).level( it4->first, &nbrlevel );
                                    nbrlevel = nbrlevel - 1;
                                    nbrkey[N - 3 * nbrlevel - 1] = 0;
                                    nbrkey[N - 3 * nbrlevel - 2] = 0;
                                    nbrkey[N - 3 * nbrlevel - 3] = 0;

                                    //            cout <<RED<< N - 3 * nbrlevel - 1 << " " << mylevel << RESET<<endl;
                                    //                                      auto it4=(*it2).find(nbrkey);
                                    //                                      (*it).level(it4->first,&nbrlevel);
                                }

                                // cout << "nbrkey " << RED << nbrkey << "  count " << counter << RESET << endl;
                                // cout<<"inside"<<endl;
                                if ( ( *it2 ).find( nbrkey ) == ( *it2 ).end() )
                                {
                                    cout << GREEN << nbrkey << "direction " << direction << "level " << nbrlevel << RESET << endl;
                                    throw std::runtime_error( "error in finding key" RESET );
                                }

                                nbrcomplevel = nbrlevel + nbrseedlevel;

                                //                             cout<<RED<<nbrlevel<<" "<<nbrseedlevel<<RESET<<endl;
                                //                             cout<<YELLOW<<nbrcomplevel<<" " <<combinedlevel<<RESET<<endl;
                                // if not in the list add to the list, dont forget to modify here
                                // watch out this element might have already been tagged

                                if ( nbrcomplevel < combinedlevel && ( *it2 ).isInRefineList( nbrkey ) == false )
                                {
                                    ( *it2 ).addToList( nbrkey );
                                    innerloop = true;
                                    // here is where boolwhile is affected
                                    loopbool = true;
                                    //       cout << RED "======================================" << endl;
                                    //       cout << "added to list" << nbrkey << endl;
                                    //       cout << "======================================" RESET << endl;
                                }
                                // cout<<seednbrkey<<RED<<seednbrkey<<RESET<<endl;
                            }

                            else
                            {
                                // sorts for communication

                                // prepare for communication, pack all the elements with the destination, sender is   my_rank

                                // find seednbrkey from global data
                                auto it5 = proc.find( seednbrkey );
                                //                        dest.push_back( it5->second[0] );
                                auto it6 = destination.begin();
                                if ( it5 != proc.end() )
                                {
                                    it6 = find( destination.begin(), destination.end(), it5->second[0] );
                                }
                                else
                                {
                                    std::runtime_error( "Not found" );
                                }
                                // combined key should be avoided
                                idx = it6 - destination.begin();
                                // cout << "---------------------->>> index" << RED << idx << RESET << endl;
                                // before pushing back, need to be careful about the key with all zeros,
                                // to accommodate this, the first bit (far right) is flipped and a sibling is sent
                                removeAllZeroSingularity( combinedkey, combinedlevel );
                                // cout << GREEN << " proc_id " << Com.myrank << "seednbrkey  " << seednbrkey << "combinedkey  " <<
                                // combinedkey
                                //     << "combined level " << combinedlevel << endl;
                                message[idx].push_back( combinedkey );

                                // myfile<< " combined key "<< combinedkey <<endl;
                                //                       cout <<"message "<< message[i].at(0) << RESET << endl;
                            }
                        }
                        combinedkey = ktcom;
                    }
                    // here change the int value for ith element of the it-th tree

                    ( *it ).flipRefineElemTag( i );
                }
            }
#endif

            counter5++;
            it3 = std::next( it3, 1 );
        }

#endif

// do a neighborhood comm to tell the nbrs the size of the message to be recieved

#if ( 1 )
        for ( uint i = 0; i < destination.size(); i++ )
        {
            //    sendbuff[i] = new morton<M + N>[ message[i].size() ];
            /*
                        for ( uint j = 0; j < message[i].size(); j++ )
                        {
                            sendbuff[i][j] = message[i].at( j );
                        }
            */
            // sending with boolean is supposed to be more portable
            sendbuff[i] = new char[message[i].size() * ( M + N )];
            // counter6=0;

            for ( uint j = 0; j < message[i].size(); j++ )
            {
                tempkey = message[i].at( j );

                //     myfile<<" message to be sent "<<tempkey<<endl;
                for ( uint k = 0; k < M + N; k++ )
                {
                    if ( tempkey[k] == true )
                    {
                        sendbuff[i][j * ( M + N ) + k] = '1';
                    }
                    else
                    {
                        sendbuff[i][j * ( M + N ) + k] = '0';
                    }
                }
            }

            //            cout << "================================== " << endl;
            //            cout << RESET "messagesize " << message[i].size() << " destination " << destination.at( i ) << RESET << endl;
            //            cout << "================================== " << endl;
            // send messages to destinations and tag each meaasge with self rank (myrank)
            //
            MPI_Isend( sendbuff[i], message[i].size() * ( M + N ), MPI_CHAR, destination.at( i ), Com.myrank, MPI_COMM_WORLD,
                       &request1[i] );
        }

        for ( uint i = 0; i < destination.size(); i++ )
        {
            //             MPI_Iprobe(destination.at(i), destination.at(i), MPI_COMM_WORLD,&flag, &status);
            // this is a symmetric comm pattern therefore any
            //
            MPI_Probe( destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &status );

            //   MPI_Probe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
            // while(flag!=0)
            {
                // MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&flag, &status);
                //  cout << BLUE << flag << RESET << endl;
                // if ( flag == 1 )
                {
                    MPI_Get_count( &status, MPI_CHAR, &size );

                    // recvbuff[i] = new morton<M + N>[ size ];
                    // notice we revcieved the message by byte so need to specify,  to find the number of elements simly divide it by sizeof
                    // bool
                    recvbuff[i] = new char[size];
                    //           cout << "******************* " << size << endl;

                    MPI_Irecv( recvbuff[i], size, MPI_CHAR, destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &request[i] );

                    MPI_Wait( &request[i], &status );
                    // if ( size != 0 )
                    {
                        //                cout << "rank " << Com.myrank << "recvbuff ??????????????? " << recvbuff[0][0] << endl;

                        // if that violates the balance condition add it to the list in the appropriate location

                        for ( uint j = 0; j < size / ( M + N ); j++ )
                        {
                            rcvkey = 0;

                            for ( uint k = 0; k < N + M; k++ )
                            {
                                if ( recvbuff[i][j * ( N + M ) + k] == '1' )
                                {
                                    rcvkey.flip( k );
                                }
                            }

                            //                   cout << "p_id " << Com.myrank << " " << rcvkey << endl;

                            //  myfile << " p_id  " << Com.myrank << " " << rcvkey << endl;

                            // restore if modified to handle singularity
                            // restoration needsto be done before this func

                            combinedLevel( rcvkey, &recvmessagecombinedlevel );

                            recoverAllZeroSingularity( rcvkey, recvmessagecombinedlevel );

                            findSeedLevelForRcvdMessage( rcvkey, &seedlevel, proc );

                            //                           cout << rcvkey << endl;
                            //                            cout << "seedlevel " << seedlevel << endl;

                            constructSeedKeyForRcvdMessage( rcvkey, seedlevel, seedkey );
                            //                                                      cout << "seedlevel " << seedlevel << endl;

                            /*      	    uint check;
                                               proc.level(seedkey,&check);
                                               cout<<RED<<check<<" "<<seedlevel<<RESET<<endl;
                                               if(seedlevel!=check)
                                              {
                                               throw std::runtime_error("levels calculated for seed inconsistent");
                                               }
                              */
                            constructElementKeyForRcvdMessage( rcvkey, seedlevel, elementkey );
                            //            cout << "element key " << elementkey << endl;

                            // cout << RED "rcvkey level " << recvmessagecombinedlevel << RESET << endl;

                            auto it7 = std::find( seeds.begin(), seeds.end(), seedkey );

                            if ( it7 == seeds.end() )
                            {
                               // cout << RED << "myrank " << Com.myrank << " " << seedkey << RESET << endl;

                                throw std::runtime_error( RED "seed not found in Balance Comm" RESET );
                            }

                            idex2 = std::distance( seeds.begin(), it7 );

                            auto it9 = std::next( trees.begin(), idex2 );

                            if ( recvmessagecombinedlevel >= seedlevel )
                            {
                                elementlevel = recvmessagecombinedlevel - seedlevel;
                            }
                            else
                            {
                                elementlevel = 1;
                            }
                            // cout << BLUE << elementlevel << RESET << endl;

                            // myfile <<" recv message combined level " <<  recvmessagecombinedlevel <<" seed level "<< seedlevel << endl;
                            // myfile <<" element level " << elementlevel<< " elementkey  "<< elementkey  << endl;

                            // recover the singularity
                            // error here, I flip the code such that we can get the level
                            // myfile << rcvkey << " " << recvmessagecombinedlevel << endl;
                            // myfile << elementkey << endl;
                            // myfile << seedkey << " " << seedlevel << endl;
                            // myfile << elementlevel << endl;
                            ( *it9 ).level( elementkey, &elementlevel );
                            // myfile << elementlevel << endl;

                            //                            recoverAllZeroSingularity( elementkey, elementlevel );

                            // get the real level
                            if ( ( *it9 ).isInMeshList( elementkey ) == false )
                            {
                                //   cout << "inside mesh" << elementkey << endl;
                                elementkey[N - 3 * ( elementlevel - 1 ) - 1] = 0;

                                elementkey[N - 3 * ( elementlevel - 1 ) - 2] = 0;

                                elementkey[N - 3 * ( elementlevel - 1 ) - 3] = 0;
                            }
                            //     myfile<<elementkey<<endl;
                            ( *it9 ).level( elementkey, &elementlevel );

                            localcombinedlevel = elementlevel + seedlevel;

                            // add if the level is lower

                            if ( localcombinedlevel < recvmessagecombinedlevel && ( *it9 ).isInRefineList( elementkey ) == false )
                            {
                                //                cout << BLUE << localcombinedlevel << "            " << recvmessagecombinedlevel << RESET
                                // << endl;

                                //               cout << "elementkey" << elementkey << endl;
                                ( *it9 ).addToList( elementkey );
                                loopbool = true;
                                /*                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;

                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                */
                            }

                            // see if this element exits, if yes, we are ok, if not add to tree referenced by pointer *it9
                        }
                    }
                }
            }
        }

#endif

        if ( loopbool == 1 )
        {
            sb = 1;
        }
        else
        {
            sb = 0;
        }

        // this value is boolan,use MPI_BYTE, since bool size is implementation dependent, use sizeof(bool)
        // nor portable switched to short

        //        MPI_Iallreduce( &sb, &rb, sizeof( bool ), MPI_BYTE, MPI_MAX, MPI_COMM_WORLD, &request0 );

        MPI_Iallreduce( &sb, &rb, 1, MPI_SHORT, MPI_MAX, MPI_COMM_WORLD, &request0 );
        // do some computation here
        //

        counter3 = 0;

        // for debug for now do the refinement separately
        /*
                for ( auto it = trees.begin(); it != trees.end(); it++ )
                {
                    start.at( counter3 ) = end.at( counter3 );
                    end.at( counter3 ) = ( *it ).refineListSize();
                    counter3++;
                }
        */

        MPI_Wait( &request0, MPI_STATUS_IGNORE );

        for ( uint j = 0; j < destination.size(); j++ )
        {
            MPI_Wait( &request1[j], &status );
            delete[] sendbuff[j];
            message[j].clear();
            message[j].shrink_to_fit();
            delete[] recvbuff[j];
        }
        //
        // need to use the recv buffer
        //
        if ( rb == 1 )
        {
            loopbool = true;
            for ( auto it = trees.begin(); it != trees.end(); it++ )
            {
                ( *it ).refinelistReset();
            }
        }
        if ( loopbool == 1 )
        {
            //                       cout << RED << " loopbool " << loopbool << " rb " << rb << " sb " << sb << RESET << endl;
        }
    }

    //    myfile.close();

    //   start.clear();
    //   end.clear();
    // assigne recieved elements to corresponding lists if needed
}

#elif( METHOD == 3 )
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::fourToOneBalance( T &proc )
{
    // cout << "character solving" << endl;
    const uint NM = N + M;
    morton<NM> combinedkey = 0, ktcom;
    morton<N> key, kt1, nbrkey;
    morton<M> seedkey, seednbrkey = 0, kt, kt3;
    auto it3 = seeds.begin();
    uint index;
    bool bol;
    real xyz[6];
    bitvector<N> boundaryElem;
    uint mylevel, changedirectionlevel, direction;
    uint topologylevel, nbrseedlevel, nbrlevel, complevel, nbrcomplevel, nbrtopologylevel, combinedlevel, localcombinedlevel;
    uint counter;
    vector<uint> directions;
    it3 = seeds.begin();
    //    vector<bitset<M + N>> message[destination.size()];
    /*vector<vector<bitset<M + N>>> message(destination.size());

    for(int i=0;i<destination.size();i++)
     {
      message[i].push_back();
     }
*/
    //  vector<bitset<M + N>> *message = new vector<bitset<M + N>>[destination.size()];
    //  message = new vector<bitset<M + N>>[destination.size()];

    //  message[0].push_back(0);
    //
    vector<uint> dest;
    vector<uint> source;
    uint idx;
    int size;
    int flag = 1;
    morton<M + N> rcvkey;
    uint seedlevel;
    morton<N> elementkey;
    uint idex2;
    uint elementlevel, recvmessagecombinedlevel;
    //    morton<M + N> *       sendbuff[destination.size()];
    //    morton<M + N> *       recvbuff[destination.size()];
    vector<uint> start;
    vector<uint> end;
    uint istart, iend;
    uint counter3, counter2, counter5;
    MPI_Request request0;

    //    MPI_Request *request=new MPI_Request[destination.size()];
    //    MPI_Request *request1=new MPI_Request[destination.size()];
    //    request = malloc( sizeof(MPI_REQUEST)*destination.size());
    //    request1 = malloc( sizeof(MPI_REQUEST)*destination.size());

    MPI_Status status;
    short sb = 0, rb = 0;
    bool mybool;
    morton<N + M> tempkey;
#if ( DEBUG )
    std::string filename = "rank";
    filename.append( to_string( Com.myrank ) );
    ofstream myfile;
    myfile.open( filename );

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        ////  cout << RED << "prod id " << Com.myrank << "  seed " << ( *it ) << RESET << endl;
        //  myfile << "prod id " << Com.myrank << "  " << con++ << "  seed " << ( *it ) << endl;
    }
    myfile << " topology size " << proc.size() << endl;
    myfile << " destination size " << destination.size() << endl;

    for ( uint i = 0; i < destination.size(); i++ )
    {
        myfile << " destinations " << destination.at( i ) << endl;
    }

#endif

    uint con = 0;
    char *sendbuff[destination.size()];
    char *recvbuff[destination.size()];

    for ( uint i = 0; i < destination.size(); i++ )
    {
        sendbuff[i] = (char *)malloc( 1 * sizeof( char ) );
        recvbuff[i] = (char *)malloc( 1 * sizeof( char ) );
    }

    // the main while loop

    //  cout << "size start=" << start.size() << endl;
    //   cout << "destination ize " << destination.size() << endl;
    //   cout << "size end=" << end.size() << endl;

    bool loopbool = true;
    uint bcstart, bcend;
    bool innerloop = true;

    int intracount = 0;

    bitvector<N> initialList[trees.size()];

    while ( loopbool )
    {
        // start out by setting loopbool=0
        mybool = false;
        loopbool = false;
        counter2 = 0;

        // save the originel list to overlap comp and com

        uint counter4 = 0;
        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {

            for ( auto i = ( *it ).Rbegin(); i != ( *it ).Rend(); i++ )
            {
                auto p = ( *it ).readRefineList( i );
                initialList[counter4].push_back( p.first );
            }
            counter4++;
        }

        // enforce 4: balance for each tree and extract the boundary elements from the list of elements to be refined

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            //            ( *it ).fourToOneP( start.at( counter2 ), end.at( counter2 ) );

            ( *it ).fourToOne();
            ( *it ).refinelistReset();
            counter2++;
        }

        // separates the list between local boundary and non-local boundary
        it3 = seeds.begin();

#if ( 1 )
        counter5 = 0;

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            seedkey = ( *it3 );
            proc.level( seedkey, &topologylevel );

#if ( DEBUG )
            myfile << " ============================== " << endl;
            myfile << " seedkey " << seedkey << "level " << topologylevel << endl;
// cout<<GREEN<<"Top lovel "<<topologylevel<<"seedkey "<<seedkey<<RESET<<endl;
#endif
            ( *it ).refinelistReset();

#if ( 1 )
            innerloop = true;
            while ( innerloop )
            {
                innerloop = false;

                for ( auto i = ( *it ).Rbegin(); i != ( *it ).Rend(); i++ )
                {
                    combinedkey = 0;

                    for ( uint j = 0; j < 3 * topologylevel; j++ )
                    {
                        combinedkey[NM - j - 1] = seedkey[M - j - 1];
                    }

                    auto p = ( *it ).readRefineList( i );
                    // if the tag = 1 then we have already investigated this guy, continue
                    if ( p.second == 1 )
                    {
                        continue;
                    }
                    //                    cout<<RED<<p.first<<endl;

                    key = p.first;

                    ( *it ).level( key, &mylevel );

#if ( DEBUG )
                    myfile << " elem key " << key << "level " << mylevel << endl;
#endif
                    combinedlevel = topologylevel + mylevel;

                    for ( uint j = 0; j < 3 * mylevel; j++ )
                    {
                        combinedkey[NM - 3 * ( topologylevel ) - j - 1] = key[N - j - 1];
                    }

                    ktcom = combinedkey;

#if ( 1 )
                    for ( uint direction = 0; direction < 3; direction++ )
                    {
                        if ( !( *it ).isBoundary( key, direction ) )
                        {
                            continue;
                        }
                        //  cout << " direction " << direction << endl;
                        //   cout << " combined key " << combinedkey << endl;
                        findFlipLevel( combinedkey, &combinedlevel, &changedirectionlevel, &direction );

#if ( DEBUG )
                        myfile << " assembled key " << combinedkey << endl;
#endif

                        if ( changedirectionlevel != 0 )
                        {
                            flipForNbr( combinedkey, &combinedlevel, &changedirectionlevel, &direction );

#if ( DEBUG )
                            myfile << " flipped key " << combinedkey << " change level " << changedirectionlevel << " direction "
                                   << direction << endl;
#endif
#if ( 1 )
                            //     cout << " combined key " << combinedkey << " combined level " << combinedlevel << endl;
                            //                     getNbrSeedLevel( combinedkey, topologylevel, &nbrseedlevel, proc ); //toplogylevel is
                            // wrong
                            getNbrSeedLevel( combinedkey, maxProcLevel, &nbrseedlevel, proc );

#if ( DEBUG )
                            myfile << " nbrseedlevel  " << nbrseedlevel << endl;
#endif

                            seednbrkey = 0;
                            for ( uint k = 0; k < 3 * nbrseedlevel; k++ )
                            {
                                seednbrkey[M - 1 - k] = combinedkey[NM - k - 1];
                                // cout<<RED<<combinedkey[NM-k-1]<<RESET<<endl;
                            }

//   cout << "seednbrkey " << BLUE << seednbrkey << RESET << endl;
#if ( DEBUG )
                            myfile << "  nbr seed key  " << seednbrkey << " nbrseedlevel  " << nbrseedlevel << endl;
#endif

                            nbrkey = 0;

                            for ( uint k = 0; k < N; k++ )
                            {
                                nbrkey[N - k - 1] = combinedkey[NM - 3 * ( nbrseedlevel ) - 1 - k];
                            }
//  cout << " nbrkey " << nbrkey << endl;
#if ( DEBUG )
                            myfile << "  nbr key  " << nbrkey << endl;
                            myfile << "Do I own It  " << isInSeed( seednbrkey, &counter ) << endl;
#endif

                            if ( isInSeed( seednbrkey, &counter ) )
                            {
                                auto it2 = std::next( trees.begin(), counter );
                                auto it4 = ( *it2 ).find( nbrkey );

                                // seedlevel
                                // ( proc ).level( seednbrkey, &nbrlevel );
                                // initial guess on nbr level as constructed
                                // this is a singularity for level zero, it might be negative for first level
                                // I remove this mannually
                                if ( combinedlevel >= nbrseedlevel )
                                {
                                    nbrlevel = combinedlevel - nbrseedlevel;
                                }
                                else
                                {
                                    nbrlevel = 1;
                                }
                                //      cout<<RED<<nbrlevel<<RESET<<endl;

                                // now find the real levels

                                if ( it4 != ( *it2 ).end() )
                                {
                                    ( *it2 ).level( it4->first, &nbrlevel );
                                    //         cout << YELLOW << nbrlevel << RESET << endl;
                                }
                                else
                                {
                                    // this condition imples level is lower no need for modification and search, I commented it out
                                    //                                nbrlevel                     = mylevel-1;
                                    // ( *it2 ).level( it4->first, &nbrlevel );
                                    nbrlevel = nbrlevel - 1;
                                    nbrkey[N - 3 * nbrlevel - 1] = 0;
                                    nbrkey[N - 3 * nbrlevel - 2] = 0;
                                    nbrkey[N - 3 * nbrlevel - 3] = 0;

                                    //            cout <<RED<< N - 3 * nbrlevel - 1 << " " << mylevel << RESET<<endl;
                                    //                                      auto it4=(*it2).find(nbrkey);
                                    //                                      (*it).level(it4->first,&nbrlevel);
                                }

                                // cout << "nbrkey " << RED << nbrkey << "  count " << counter << RESET << endl;
                                // cout<<"inside"<<endl;
                                if ( ( *it2 ).find( nbrkey ) == ( *it2 ).end() )
                                {
                                    //cout << Com.myrank << " " << nbrkey << " direction " << direction << " level " << nbrlevel << endl;
                                    throw std::runtime_error( "error in finding key refinement" RESET );
                                }

                                nbrcomplevel = nbrlevel + nbrseedlevel;

                                //                             cout<<RED<<nbrlevel<<" "<<nbrseedlevel<<RESET<<endl;
                                //                             cout<<YELLOW<<nbrcomplevel<<" " <<combinedlevel<<RESET<<endl;
                                // if not in the list add to the list, dont forget to modify here
                                // watch out this element might have already been tagged

                                if ( nbrcomplevel < combinedlevel && ( *it2 ).isInRefineList( nbrkey ) == false )
                                {
                                    ( *it2 ).addToList( nbrkey );
                                    innerloop = true;
                                    // here is where boolwhile is affected
                                    loopbool = true;
                                    //       cout << RED "======================================" << endl;
                                    //       cout << "added to list" << nbrkey << endl;
                                    //       cout << "======================================" RESET << endl;
                                }
                                // cout<<seednbrkey<<RED<<seednbrkey<<RESET<<endl;
                            }

                            else
                            {
// sorts for communication

#if ( 1 )
// prepare for communication, pack all the elements with the destination, sender is   my_rank
#if ( DEBUG )

                                myfile << "*********************************  " << endl;
                                myfile << "inside else  nbr seed key  " << seednbrkey << " nbrseedlevel  " << nbrseedlevel << endl;
#endif
                                // find seednbrkey from global data
                                appendToMessage( proc, seednbrkey, combinedkey, combinedlevel );

#endif
                            }
                        }
                        combinedkey = ktcom;
                    }
// here change the int value for ith element of the it-th tree
#endif
                    ( *it ).flipRefineElemTag( i );
                }
            }
#endif

#endif
            counter5++;
            it3 = std::next( it3, 1 );
        }

#endif

// do a neighborhood comm to tell the nbrs the size of the message to be recieved

#if ( 1 )
        for ( uint i = 0; i < destination.size(); i++ )
        {
           // sending with boolean is not portable, use characters

            //    sendbuff[i] = new char[message[i].size() * ( M + N )];

            sendbuff[i] = (char *)realloc( sendbuff[i], ( message[i].size() * ( M + N ) + 1 ) * sizeof( char ) );
            // counter6=0;
            //
            // if(sendbuff[i]==NULL && message[i].size()!=0)
            if ( sendbuff[i] == NULL )
            {
               // cout << " myrank " << Com.myrank << " " << message[i].size() << endl;
                throw std::runtime_error( "not able to allocate sendbuf" );
            }
            for ( uint j = 0; j < message[i].size(); j++ )
            {
                tempkey = message[i].at( j );

#if ( DEBUG )
                myfile << " message to be sent " << tempkey << " destination " << destination.at( i ) << endl;
#endif
                for ( uint k = 0; k < M + N; k++ )
                {
                    if ( tempkey[k] == true )
                    {
                        sendbuff[i][j * ( M + N ) + k] = '1';
                    }
                    else
                    {
                        sendbuff[i][j * ( M + N ) + k] = '0';
                    }
                }
            }

            //            cout << "================================== " << endl;
            //            cout << RESET "messagesize " << message[i].size() << " destination " << destination.at( i ) << RESET << endl;
            //            cout << "================================== " << endl;
            // send messages to destinations and tag each meaasge with self rank (myrank)
            //
            MPI_Isend( sendbuff[i], message[i].size() * ( M + N ), MPI_CHAR, destination.at( i ), Com.myrank, MPI_COMM_WORLD,
                       &request1[i] );
        }

        for ( uint i = 0; i < destination.size(); i++ )
        {
            //             MPI_Iprobe(destination.at(i), destination.at(i), MPI_COMM_WORLD,&flag, &status);
            // this is a symmetric comm pattern therefore any
            //
            MPI_Probe( destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &status );

            MPI_Get_count( &status, MPI_CHAR, &size );

            // recvbuff[i] = new morton<M + N>[ size ];
            // notice we revcieved the message by byte so need to specify,  to find the number of elements simly divide it by sizeof
            // bool
            //   recvbuff[i] = new char[size];
            //   size+1 is to ensure deallocation
            recvbuff[i] = (char *)realloc( recvbuff[i], ( size + 1 ) * sizeof( char ) );

            // if(recvbuff[i]==NULL && size!=0)
            if ( recvbuff[i] == NULL )
            {
              //  cout << " myrank " << Com.myrank << " " << size << endl;
                throw std::runtime_error( "bad alloc in recv buffer" );
            }

            //           cout << "******************* " << size << endl;

            MPI_Irecv( recvbuff[i], size, MPI_CHAR, destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &request[i] );

            MPI_Wait( &request[i], &status );
            // if ( size != 0 )
            //                cout << "rank " << Com.myrank << "recvbuff ??????????????? " << recvbuff[0][0] << endl;

            // if that violates the balance condition add it to the list in the appropriate location

            for ( uint j = 0; j < size / ( M + N ); j++ )
            {
                rcvkey = 0;

                for ( uint k = 0; k < N + M; k++ )
                {
                    if ( recvbuff[i][j * ( N + M ) + k] == '1' )
                    {
                        rcvkey.flip( k );
                    }
                }

//             cout << "p_id " << Com.myrank << " " << rcvkey << endl;
#if ( DEBUG )
                myfile << " p_id  " << Com.myrank << " " << rcvkey << endl;
#endif
                // restore if modified to handle singularity
                // restoration needsto be done before this func

                combinedLevel( rcvkey, &recvmessagecombinedlevel );

                recoverAllZeroSingularity( rcvkey, recvmessagecombinedlevel );

                findSeedLevelForRcvdMessage( rcvkey, &seedlevel, proc );

                //                           cout << rcvkey << endl;
                //                            cout << " seedlevel " << seedlevel << endl;

                constructSeedKeyForRcvdMessage( rcvkey, seedlevel, seedkey );
                // cout <<" myrank "<<Com.myrank <<" seedlevel " << seedlevel << endl;

                /*      	    uint check;
                                   proc.level(seedkey,&check);
                                   cout<<RED<<check<<" "<<seedlevel<<RESET<<endl;
                                   if(seedlevel!=check)
                                  {
                                   throw std::runtime_error("levels calculated for seed inconsistent");
                                   }
                  */
                constructElementKeyForRcvdMessage( rcvkey, seedlevel, elementkey );
                //            cout << "element key " << elementkey << endl;

                // cout << RED "rcvkey level " << recvmessagecombinedlevel << RESET << endl;

                auto it7 = std::find( seeds.begin(), seeds.end(), seedkey );

                if ( it7 == seeds.end() )
                {
                  //  cout << " p_id  " << Com.myrank << " " << rcvkey << "  " << recvmessagecombinedlevel << endl;
                  //  cout << RED << "myrank " << Com.myrank << " " << seedkey << RESET << endl;
                    //  cout << "myrank " << Com.myrank << " " << seedkey << endl;

                    throw std::runtime_error( RED "seed not found in Balance Comm" RESET );
                }

                idex2 = std::distance( seeds.begin(), it7 );

                auto it9 = std::next( trees.begin(), idex2 );

                if ( recvmessagecombinedlevel >= seedlevel )
                {
                    elementlevel = recvmessagecombinedlevel - seedlevel;
                }
                else
                {
                    elementlevel = 1;
                }

#if ( DEBUG )
                myfile << " recv message combined level " << recvmessagecombinedlevel << " seed level " << seedlevel << endl;
                myfile << " element level " << elementlevel << " elementkey  " << elementkey << endl;

                recover the singularity
                // error here, I flip the code such that we can get the level
                myfile
                << rcvkey << " " << recvmessagecombinedlevel << endl;
                myfile << elementkey << endl;
                myfile << seedkey << " " << seedlevel << endl;
                myfile << elementlevel << endl;
#endif
                ( *it9 ).level( elementkey, &elementlevel );
                // myfile << elementlevel << endl;

                //                            recoverAllZeroSingularity( elementkey, elementlevel );

                // get the real level
                if ( ( *it9 ).isInMeshList( elementkey ) == false )
                {
                    //   cout << "inside mesh" << elementkey << endl;
                    elementkey[N - 3 * ( elementlevel - 1 ) - 1] = 0;

                    elementkey[N - 3 * ( elementlevel - 1 ) - 2] = 0;

                    elementkey[N - 3 * ( elementlevel - 1 ) - 3] = 0;
                }
                //     myfile<<elementkey<<endl;
                ( *it9 ).level( elementkey, &elementlevel );

                localcombinedlevel = elementlevel + seedlevel;

                // add if the level is lower

                if ( localcombinedlevel < recvmessagecombinedlevel && ( *it9 ).isInRefineList( elementkey ) == false )
                {
                    //                cout << BLUE << localcombinedlevel << "            " << recvmessagecombinedlevel << RESET
                    // << endl;

                    //               cout << "elementkey" << elementkey << endl;
                    ( *it9 ).addToList( elementkey );
                    loopbool = true;

                    mybool = true;
                }

                // see if this element exits, if yes, we are ok, if not add to tree referenced by pointer *it9
            }
        }

#endif

#if ( 1 )

        if ( mybool == true )
        {
            intracount++;
        }

        if ( loopbool == 1 )
        {
            sb = 1;
        }
        else
        {
            sb = 0;
        }

        // use MPI_SHORT as boolean is implementation dependent as transfering MPI_BYTE is not robust due to big endian little endian issue

        MPI_Iallreduce( &sb, &rb, 1, MPI_SHORT, MPI_MAX, MPI_COMM_WORLD, &request0 );

        // do some computation here
        //
        counter3 = 0;

        // switch the tag for initial list
if(OVERLAP==1)
{
        auto it = trees.begin();

        for ( uint i1 = 0; i1 < trees.size(); i1++ )
        {
            ( *it ).refineRefineList( initialList[i1] );
            it = std::next( it, 1 );
        }

        for ( uint i = 0; i < trees.size(); i++ )
        {
            initialList[i].clear();
        }
}
        // for debug for now do the refinem
        MPI_Wait( &request0, MPI_STATUS_IGNORE );

        for ( uint j = 0; j < destination.size(); j++ )
        {
            MPI_Wait( &request1[j], &status );

            while ( message[j].size() > 0 )
            {
                //   MPI_Wait( &request1[j], &status );
                message[j].clear();

                // delete[] sendbuff[j];
            }

            //     message[j].shrink_to_fit();
            //    delete[] recvbuff[j];
        }

#endif
        //
        // need to use the recv buffer
        //
        if ( rb == 1 )
        {
            loopbool = true;
            for ( auto it = trees.begin(); it != trees.end(); it++ )
            {
                ( *it ).refinelistReset();
            }
        }
        if ( loopbool == 1 )
        {
#if ( 1 )
            cout << RED << " loopbool " << loopbool << " rb " << rb << " sb " << sb << RESET << endl;
            cout<<" count " <<intracount<<endl;
#endif
        }
    }

            cout <<  " -------------------------------------- " << endl;
#if ( DEBUG )
    myfile.close();
#endif

    for ( int i = 0; i < destination.size(); i++ )
    {
        free( sendbuff[i] );
        free( recvbuff[i] );
    }

    //   start.clear();
    //   end.clear();
    // assigne recieved elements to corresponding lists if needed

    // cout<<"============================"<<RED<<Com.myrank<<" "<<intracount<<endl;
}

#elif( METHOD == 4 )
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::fourToOneBalance( T &proc )
{
    // cout << "character solving" << endl;
    const uint NM = N + M;
    morton<NM> combinedkey = 0, ktcom;
    morton<N> key, kt1, nbrkey;
    morton<M> seedkey, seednbrkey = 0, kt, kt3;
    auto it3 = seeds.begin();
    uint index;
    bool bol;
    real xyz[6];
    bitvector<N> boundaryElem;
    uint mylevel, changedirectionlevel, direction;
    uint topologylevel, nbrseedlevel, nbrlevel, complevel, nbrcomplevel, nbrtopologylevel, combinedlevel, localcombinedlevel;
    uint counter;
    vector<uint> directions;
    it3 = seeds.begin();
    vector<uint> dest;
    vector<uint> source;
    uint idx;
    int size;
    int flag = 1;
    morton<M + N> rcvkey;
    uint seedlevel;
    morton<N> elementkey;
    uint idex2;
    uint elementlevel, recvmessagecombinedlevel;
    vector<uint> start;
    vector<uint> end;
    uint istart, iend;
    uint counter3, counter2, counter5;
    MPI_Request request0;
    MPI_Status status;
    short sb = 0, rb = 0;
    bool mybool;
    morton<N + M> tempkey;
#if ( DEBUG )
    std::string filename = "rank";
    filename.append( to_string( Com.myrank ) );
    ofstream myfile;
    myfile.open( filename );

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        ////  cout << RED << "prod id " << Com.myrank << "  seed " << ( *it ) << RESET << endl;
        //        myfile << "prod id " << Com.myrank << "  " << con++ << "  seed " << ( *it ) << endl;
    }
    myfile << " topology size " << proc.size() << endl;
#endif

    int *dispS = nullptr;
    int *dispR = nullptr;
    uint con = 0;
    //    char *sendbuff[destination.size()];
    //    char *recvbuff[destination.size()];
    char *sendbuff = nullptr;
    char *recvbuff = nullptr;

    sendbuff = new char[1];
    recvbuff = new char[1];
    dispS = new int[1];
    dispR = new int[1];

    int sendSize[destination.size()];
    // for debug
    int recvSize[destination.size()];

    // the main while loop

    //  cout << "size start=" << start.size() << endl;
    //   cout << "destination ize " << destination.size() << endl;
    //   cout << "size end=" << end.size() << endl;

    bool loopbool = true;
    uint bcstart, bcend;
    bool innerloop = true;

    int intracount = 0;

    bitvector<N> initialList[trees.size()];

    while ( loopbool )
    {
        // start out by setting loopbool=0
        mybool = false;
        loopbool = false;
        counter2 = 0;

        // save the originel list to overlap comp and com

        uint counter4 = 0;
        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {

            for ( auto i = ( *it ).Rbegin(); i != ( *it ).Rend(); i++ )
            {
                auto p = ( *it ).readRefineList( i );
                initialList[counter4].push_back( p.first );
            }
            counter4++;
        }

        // enforce 4: balance for each tree and extract the boundary elements from the list of elements to be refined

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            //            ( *it ).fourToOneP( start.at( counter2 ), end.at( counter2 ) );

            ( *it ).fourToOne();
            ( *it ).refinelistReset();
            counter2++;
        }

        // separates the list between local boundary and non-local boundary
        it3 = seeds.begin();

#if ( 1 )
        counter5 = 0;

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            seedkey = ( *it3 );
            proc.level( seedkey, &topologylevel );

#if ( DEBUG )
            myfile << " ============================== " << endl;
            myfile << " seedkey " << seedkey << "level " << topologylevel << endl;
// cout<<GREEN<<"Top lovel "<<topologylevel<<"seedkey "<<seedkey<<RESET<<endl;
#endif
            ( *it ).refinelistReset();

            innerloop = true;
            while ( innerloop )
            {
                innerloop = false;

                for ( auto i = ( *it ).Rbegin(); i != ( *it ).Rend(); i++ )
                {
                    combinedkey = 0;

                    for ( uint j = 0; j < 3 * topologylevel; j++ )
                    {
                        combinedkey[NM - j - 1] = seedkey[M - j - 1];
                    }

                    auto p = ( *it ).readRefineList( i );
                    // if the tag = 1 then we have already investigated this guy, continue
                    if ( p.second == 1 )
                    {
                        continue;
                    }
                    //                    cout<<RED<<p.first<<endl;

                    key = p.first;

                    ( *it ).level( key, &mylevel );

#if ( DEBUG )
                    myfile << " elem key " << key << "level " << mylevel << endl;
#endif
                    combinedlevel = topologylevel + mylevel;

                    for ( uint j = 0; j < 3 * mylevel; j++ )
                    {
                        combinedkey[NM - 3 * ( topologylevel ) - j - 1] = key[N - j - 1];
                    }

                    ktcom = combinedkey;

                    for ( uint direction = 0; direction < 3; direction++ )
                    {
                        if ( !( *it ).isBoundary( key, direction ) )
                        {
                            continue;
                        }
                        //  cout << " direction " << direction << endl;
                        //   cout << " combined key " << combinedkey << endl;
                        findFlipLevel( combinedkey, &combinedlevel, &changedirectionlevel, &direction );

#if ( DEBUG )
                        myfile << " assembled key " << combinedkey << endl;
#endif

                        if ( changedirectionlevel != 0 )
                        {
                            flipForNbr( combinedkey, &combinedlevel, &changedirectionlevel, &direction );

#if ( DEBUG )
                            myfile << " flipped key " << combinedkey << " change level " << changedirectionlevel << " direction "
                                   << direction << endl;
#endif
#if ( 1 )

                            //     cout << " combined key " << combinedkey << " combined level " << combinedlevel << endl;
                            // getNbrSeedLevel( combinedkey, topologylevel, &nbrseedlevel, proc ); //toplogylevel is  wrong
                            getNbrSeedLevel( combinedkey, maxProcLevel, &nbrseedlevel, proc );

#if ( DEBUG )
                            myfile << " nbrseedlevel  " << nbrseedlevel << endl;
#endif

                            seednbrkey = 0;
                            for ( uint k = 0; k < 3 * nbrseedlevel; k++ )
                            {
                                seednbrkey[M - 1 - k] = combinedkey[NM - k - 1];
                                // cout<<RED<<combinedkey[NM-k-1]<<RESET<<endl;
                            }

//   cout << "seednbrkey " << BLUE << seednbrkey << RESET << endl;
#if ( DEBUG )
                            myfile << "  nbr seed key  " << seednbrkey << " nbrseedlevel  " << nbrseedlevel << endl;
#endif

                            nbrkey = 0;

                            for ( uint k = 0; k < N; k++ )
                            {
                                nbrkey[N - k - 1] = combinedkey[NM - 3 * ( nbrseedlevel ) - 1 - k];
                            }
//  cout << " nbrkey " << nbrkey << endl;
#if ( DEBUG )
                            myfile << "  nbr key  " << nbrkey << endl;
                            myfile << "Do I own It  " << isInSeed( seednbrkey, &counter ) << endl;
#endif

                            if ( isInSeed( seednbrkey, &counter ) )
                            {
                                auto it2 = std::next( trees.begin(), counter );
                                auto it4 = ( *it2 ).find( nbrkey );

                                // seedlevel
                                // ( proc ).level( seednbrkey, &nbrlevel );
                                // initial guess on nbr level as constructed
                                // this is a singularity for level zero, it might be negative for first level
                                // I remove this mannually
                                if ( combinedlevel >= nbrseedlevel )
                                {
                                    nbrlevel = combinedlevel - nbrseedlevel;
                                }
                                else
                                {
                                    nbrlevel = 1;
                                }
                                //      cout<<RED<<nbrlevel<<RESET<<endl;

                                // now find the real levels

                                if ( it4 != ( *it2 ).end() )
                                {
                                    ( *it2 ).level( it4->first, &nbrlevel );
                                    //         cout << YELLOW << nbrlevel << RESET << endl;
                                }
                                else
                                {
                                    // this condition imples level is lower no need for modification and search, I commented it out
                                    //                                nbrlevel                     = mylevel-1;
                                    // ( *it2 ).level( it4->first, &nbrlevel );
                                    nbrlevel = nbrlevel - 1;
                                    nbrkey[N - 3 * nbrlevel - 1] = 0;
                                    nbrkey[N - 3 * nbrlevel - 2] = 0;
                                    nbrkey[N - 3 * nbrlevel - 3] = 0;

                                    //            cout <<RED<< N - 3 * nbrlevel - 1 << " " << mylevel << RESET<<endl;
                                    //                                      auto it4=(*it2).find(nbrkey);
                                    //                                      (*it).level(it4->first,&nbrlevel);
                                }
#if ( DEBUG )
                                myfile << nbrkey << "direction " << direction << "level " << nbrlevel << RESET << "my fault" << endl;
#endif

                                // cout << "nbrkey " << RED << nbrkey << "  count " << counter << RESET << endl;
                                // cout<<"inside"<<endl;
                                if ( ( *it2 ).find( nbrkey ) == ( *it2 ).end() )
                                {
                                    cout << GREEN << nbrkey << "direction " << direction << "level " << nbrlevel << RESET << endl;
                                    throw std::runtime_error( "error in finding key refinement" RESET );
                                }

                                nbrcomplevel = nbrlevel + nbrseedlevel;

                                //                             cout<<RED<<nbrlevel<<" "<<nbrseedlevel<<RESET<<endl;
                                //                             cout<<YELLOW<<nbrcomplevel<<" " <<combinedlevel<<RESET<<endl;
                                // if not in the list add to the list, dont forget to modify here
                                // watch out this element might have already been tagged

                                if ( nbrcomplevel < combinedlevel && ( *it2 ).isInRefineList( nbrkey ) == false )
                                {
                                    ( *it2 ).addToList( nbrkey );
                                    innerloop = true;
                                    // here is where boolwhile is affected
                                    loopbool = true;
                                    //       cout << RED "======================================" << endl;
                                    //       cout << "added to list" << nbrkey << endl;
                                    //       cout << "======================================" RESET << endl;
                                }
                                // cout<<seednbrkey<<RED<<seednbrkey<<RESET<<endl;
                            }

                            else
                            {
                                // sorts for communication

                                appendToMessage( proc, seednbrkey, combinedkey, combinedlevel );
// prepare for communication, pack all the elements with the destination, sender is   my_rank
#if ( DEBUG )

                                myfile << "*********************************  " << endl;
                                myfile << "inside else  nbr seed key  " << seednbrkey << " nbrseedlevel  " << nbrseedlevel << endl;
                                // find seednbrkey from global data
                                auto it5 = proc.find( seednbrkey );
                                //                        dest.push_back( it5->second[0] );
                                auto it6 = destination.begin();
                                if ( it5 != proc.end() )
                                {
                                    it6 = find( destination.begin(), destination.end(), it5->second[0] );
                                }
                                else
                                {
                                    std::runtime_error( "Not found" );
                                }
                                // combined key should be avoided
                                idx = it6 - destination.begin();
                                // cout << "---------------------->>> index" << RED << idx << RESET << endl;
                                // before pushing back, need to be careful about the key with all zeros,
                                // to accommodate this, the first bit (far right) is flipped and a sibling is sent
                                //  myfile<< " combined key before "<< combinedkey << "level "<<combinedlevel <<endl;
                                removeAllZeroSingularity( combinedkey, combinedlevel );
                                // cout << GREEN << " proc_id " << Com.myrank << "seednbrkey  " << seednbrkey << "combinedkey  " <<
                                // combinedkey
                                //     << "combined level " << combinedlevel << endl;
                                message[idx].push_back( combinedkey );

#if ( DEBUG )
                                myfile << " combined key pushed back  " << combinedkey << "level " << combinedlevel << endl;
                                myfile << "*********************************  " << endl;
#endif
//                       cout <<"message "<< message[i].at(0) << RESET << endl;

#endif
                            }
                        }
                        combinedkey = ktcom;
                    }
                    // here change the int value for ith element of the it-th tree

                    ( *it ).flipRefineElemTag( i );
                }
            }
#endif

            counter5++;
            it3 = std::next( it3, 1 );
        }

#endif

        // do a neighborhood comm to tell the nbrs the size of the message to be recieved
        // using neighborhood collectives

        for ( uint i = 0; i < destination.size(); i++ )
        {
            sendSize[i] = message[i].size() * ( M + N );
            recvSize[i] = 0;
        }

        MPI_Ineighbor_alltoall( sendSize, 1, MPI_INT, recvSize, 1, MPI_INT, graphComm, &request[0] );

        // MPI_Neighbor_alltoall(sendSize, 1, MPI_INT, recvSize, 1 ,MPI_INT, graphComm);

        int totalSendSize = 0;
        int totalRecvSize = 0;

        for ( uint i = 0; i < destination.size(); i++ )
        {
            totalSendSize += sendSize[i];
        }

        // dispS=new int [destination.size()];

        dispS = (int *)realloc( dispS, destination.size() * sizeof( int ) );
        dispS[0] = 0;

        for ( int i = 1; i < destination.size(); i++ )
        {
            dispS[i] = dispS[i - 1] + sendSize[i - 1];
        }

        // sendbuff=new char[totalSendSize];

        sendbuff = (char *)realloc( sendbuff, totalSendSize * sizeof( char ) );

        if ( sendbuff == NULL && totalSendSize != 0 )
        {
            throw std::runtime_error( "bad allocation in sendbuff in METHOD 4" );
        }

        // do a neighborhood comm to tell the nbrs the size of the message to be recieved
        for ( uint i = 0; i < destination.size(); i++ )
        {

            for ( uint j = 0; j < message[i].size(); j++ )
            {
                tempkey = message[i].at( j );

                //    myfile<<" message to be sent "<<tempkey<<" destination "<<destination.at(i)<<endl;
                for ( uint k = 0; k < M + N; k++ )
                {
                    if ( tempkey[k] == true )
                    {
                        sendbuff[dispS[i] + j * ( M + N ) + k] = '1';
                    }
                    else
                    {
                        sendbuff[dispS[i] + j * ( M + N ) + k] = '0';
                    }
                }

#if ( DEBUG )
                myfile << " myrank " << Com.myrank << " message converted " << tempkey << " destination " << destination.at( i ) << endl;
#endif
            }
        }

        MPI_Wait( &request[0], &status );

#if ( DEBUG )
        for ( uint i = 0; i < destination.size(); i++ )
        {
            myfile << " rank " << Com.myrank << " send data " << sendSize[i] << " to " << destination.at( i ) << " and will rcv "
                   << recvSize[i] << endl;
        }
#endif

        for ( uint i = 0; i < destination.size(); i++ )
        {
            totalRecvSize += recvSize[i];
        }

        // myfile<<" rank "<<Com.myrank<<" total send size " <<  totalSendSize<<"total recieve  "<<totalRecvSize<<endl;

        // recvbuff=new char[totalRecvSize];
        recvbuff = (char *)realloc( recvbuff, totalRecvSize * sizeof( char ) );

        if ( recvbuff == NULL && totalRecvSize != 0 )
        {
            throw std::runtime_error( "bad allocation in sendbuff in METHOD 4" );
        }

        dispR = (int *)realloc( dispR, destination.size() * sizeof( int ) );

        // dispR=new int [destination.size()];

        dispR[0] = 0;

        for ( int i = 1; i < destination.size(); i++ )
        {
            dispR[i] = dispR[i - 1] + recvSize[i - 1];
            // myfile<<" rank "<<Com.myrank<<" dispS" <<  dispS[i]<<" dispR  "<<dispR[i]<<endl;
        }

        for ( int i = 0; i < destination.size(); i++ )
        {
#if ( DEBUG )
            myfile << " rank " << Com.myrank << " dispSend " << dispS[i] << " dispR  " << dispR[i] << endl;
#endif
        }

        // one more neighborhood collective to tranfer the data

        MPI_Ineighbor_alltoallv( sendbuff, sendSize, dispS, MPI_CHAR, recvbuff, recvSize, dispR, MPI_CHAR, graphComm, &request[0] );

#if ( DEBUG )
        myfile << " rank " << Com.myrank << " totalRecvSize " << totalRecvSize << endl;
#endif

        MPI_Wait( &request[0], &status );

        for ( uint j = 0; j < totalRecvSize / ( M + N ); j++ )
        {
            rcvkey = 0;

            for ( uint k = 0; k < N + M; k++ )
            {
                if ( recvbuff[j * ( N + M ) + k] == '1' )
                {
                    rcvkey.flip( k );
                }
            }

#if ( DEBUG )
            myfile << " rank " << Com.myrank << " recevdKey " << rcvkey << endl;
#endif
#if ( 1 )
            combinedLevel( rcvkey, &recvmessagecombinedlevel );

            recoverAllZeroSingularity( rcvkey, recvmessagecombinedlevel );

            findSeedLevelForRcvdMessage( rcvkey, &seedlevel, proc );

            constructSeedKeyForRcvdMessage( rcvkey, seedlevel, seedkey );

            constructElementKeyForRcvdMessage( rcvkey, seedlevel, elementkey );

            //            cout << "element key " << elementkey << endl;

            // cout << RED "rcvkey level " << recvmessagecombinedlevel << RESET << endl;

            auto it7 = std::find( seeds.begin(), seeds.end(), seedkey );

            if ( it7 == seeds.end() )
            {
               // cout << " p_id  " << Com.myrank << " " << rcvkey << "  " << recvmessagecombinedlevel << endl;
                // cout << RED << "myrank " << Com.myrank << " " << seedkey << RESET << endl;
                //  cout << "myrank " << Com.myrank << " " << seedkey << endl;

                throw std::runtime_error( RED "seed not found in Balance Comm" RESET );
            }

            idex2 = std::distance( seeds.begin(), it7 );

            auto it9 = std::next( trees.begin(), idex2 );

            if ( recvmessagecombinedlevel >= seedlevel )
            {
                elementlevel = recvmessagecombinedlevel - seedlevel;
            }
            else
            {
                elementlevel = 1;
            }
            // cout << BLUE << elementlevel << RESET << endl;

            // myfile <<" recv message combined level " <<  recvmessagecombinedlevel <<" seed level "<< seedlevel << endl;
            // myfile <<" element level " << elementlevel<< " elementkey  "<< elementkey  << endl;

            // recover the singularity
            // error here, I flip the code such that we can get the level
            // myfile << rcvkey << " " << recvmessagecombinedlevel << endl;
            // myfile << elementkey << endl;
            // myfile << seedkey << " " << seedlevel << endl;
            // myfile << elementlevel << endl;
            ( *it9 ).level( elementkey, &elementlevel );
            // myfile << elementlevel << endl;

            //                            recoverAllZeroSingularity( elementkey, elementlevel );

            // get the real level
            if ( ( *it9 ).isInMeshList( elementkey ) == false )
            {
                //   cout << "inside mesh" << elementkey << endl;
                elementkey[N - 3 * ( elementlevel - 1 ) - 1] = 0;

                elementkey[N - 3 * ( elementlevel - 1 ) - 2] = 0;

                elementkey[N - 3 * ( elementlevel - 1 ) - 3] = 0;
            }
            //     myfile<<elementkey<<endl;
            ( *it9 ).level( elementkey, &elementlevel );

            localcombinedlevel = elementlevel + seedlevel;

            // add if the level is lower

            if ( localcombinedlevel < recvmessagecombinedlevel && ( *it9 ).isInRefineList( elementkey ) == false )
            {
                //               cout << "elementkey" << elementkey << endl;
                ( *it9 ).addToList( elementkey );
                loopbool = true;

                mybool = true;
            }

#endif
            // see if this element exits, if yes, we are ok, if not add to tree referenced by pointer *it9
        }

        if ( mybool == true )
        {
            intracount++;
        }

        if ( loopbool == 1 )
        {
            sb = 1;
        }
        else
        {
            sb = 0;
        }

        // this value is boolan,use MPI_BYTE, since bool size is implementation dependent, use sizeof(bool)

        // MPI_Iallreduce( &sb, &rb, sizeof( bool ), MPI_BYTE, MPI_MAX, MPI_COMM_WORLD, &request0 );

        // this value is boolan,use MPI_BYTE, since bool size is implementation dependent, use sizeof(bool)

        MPI_Iallreduce( &sb, &rb, 1, MPI_SHORT, MPI_MAX, MPI_COMM_WORLD, &request0 );

        // do some computation here
        //

        counter3 = 0;

        // switch the tag for initial list

        auto it = trees.begin();

        for ( uint i1 = 0; i1 < trees.size(); i1++ )
        {
            ( *it ).refineRefineList( initialList[i1] );
            it = std::next( it, 1 );
        }

        for ( uint i = 0; i < trees.size(); i++ )
        {
            initialList[i].clear();
        }

        // for debug for now do the refinem
        MPI_Wait( &request0, MPI_STATUS_IGNORE );

        for ( uint j = 0; j < destination.size(); j++ )
        {
            //      MPI_Wait( &request1[j], &status );
            //         delete[] sendbuff;
            message[j].clear();
            //     message[j].shrink_to_fit();
            //           delete[] recvbuff;
        }
#if ( 1 )
        //
        // need to use the recv buffer
        //
        if ( rb == 1 )
        {
            loopbool = true;
            for ( auto it = trees.begin(); it != trees.end(); it++ )
            {
                ( *it ).refinelistReset();
            }
        }
        if ( loopbool == 1 )
        {
            //                       cout << RED << " loopbool " << loopbool << " rb " << rb << " sb " << sb << RESET << endl;
        }

#endif

#if ( 0 )
#endif
    }

#if ( DEBUG )
    myfile.close();

#endif
    //   start.clear();
    //   end.clear();
    // assigne recieved elements to corresponding lists if needed

    // cout<<"============================"<<RED<<Com.myrank<<" "<<intracount<<endl;
}

#endif

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::refineForestBalanced( uint nlevel, T &proc )
{
    //   encodeGeometry();

    getMaxSeedsLevel( proc );

    //  cout << " max seed level " << maxseedlevel << endl;

    for ( uint i = 0; i < nlevel; i++ )
    {
        // cout<<RED" nlevel "<<nlevel<<RESET<<endl;
        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            ( *it ).pushToRefinelist( i );
        }
        fourToOneBalance( proc );

        //        MPI_Barrier(MPI_COMM_WORLD);
        if ( OVERLAP == 0 )
        {
            for ( auto it = trees.begin(); it != trees.end(); it++ )
            {
                ( *it ).refineRefineList();
                // ( *it ).boundarylist.clear();
            }
        }
    }
}

#if ( 0 )
/*
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void Forest<N, Nvalue, M, Mvalue>::removeZeroCombined( const morton<N + M> &key, const uint &seedlevel, morton<M> &seedkey )
{


}
*/

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void Forest<N, Nvalue, M, Mvalue>::debug( Tree<M, Mvalue> &proc )
{
    morton<M> seedkey = 0;
    uint counter = 0;

    if ( isInSeed( seedkey, &counter ) )
    {
        //   cout << "seedkey found " << seedkey << endl;
        //    cout << "counter is " << counter << endl;
    }
    morton<N> elem = 0;
    int myrank;

    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

    if ( myrank == 0 )
    {
        auto it = std::next( trees.begin(), counter );
        for ( auto it = seeds.begin(); it != seeds.end(); it++ )
        {
            cout << RED << ( *it ) << RESET << endl;
        }

        ( *it ).addToList( elem );
        cout << "elem_new " << elem << endl;

        for ( auto it2 = trees.begin(); it2 != trees.end(); it2++ )
        {
            ( *it2 ).refineRefineList();
        }

        elem.flip( N - 1 );
        elem.flip( N - 2 );
        elem.flip( N - 3 );

        ( *it ).addToList( elem );
        cout << "elem =" << elem << endl;
    }

    fourToOneBalance( proc );

#if ( 0 )
    {
        for ( auto it2 = trees.begin(); it2 != trees.end(); it2++ )
        {
            ( *it2 ).refineRefineList();
        }

        elem.flip( N - 4 );
        elem.flip( N - 5 );
        elem.flip( N - 6 );

        ( *it ).addToList( elem );

        extendListLocalBalance( proc );

        cout << "elem " << elem << endl;
        for ( auto it2 = trees.begin(); it2 != trees.end(); it2++ )
        {
            ( *it2 ).refineRefineList();
        }
    }
#endif
}

//====================================================================
//
//  redistribution of the mesh after generation, this is to enforce load
//  balance at the end of the preocess for solving
//
//=====================================================================

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void Forest<N, Nvalue, M, Mvalue>::zoltanGeomrepart( Tree<M, Mvalue> &proc, uint setmethod )
{
    uint myvalue = seeds.size();
    uint q = myvalue;

    CommPoint2Point<uint> com( &q, 1 );

    uint offset = 0;

    com.getOffset( q, &offset );
    // cout << RED << "my value" << myvalue << RESET << endl;
    // later on need to gather and sum up all the proc parts, need to chop the proc tree
    // for better scaling

    //       uint totalvalue=proc.size();
    uint totalvalue;
    CommCollective<uint> comc( nullptr, 1, Com.comsize - 1 );
    comc.IgetTotalNumber( &offset, &myvalue, &totalvalue );
    real *weight = nullptr;
    weight = new real[myvalue];
    Center_coords XYZ;
    auto it1 = seeds.begin();
    //    morton<N> key;
    uint co = 0;
    real xyzc[3];

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        proc.centroid( ( *it1 ), xyzc );

        XYZ.push_back( CenterCoords() );
        XYZ.at( co ).x = xyzc[0];
        XYZ.at( co ).y = xyzc[1];
        XYZ.at( co ).z = xyzc[2];
        weight[co] = ( *it ).size();
        //            cout<<" weights "<<Com.myrank<<" "<<weight[co]<<endl;
        co++;
        it1 = std::next( it1, 1 );
    }

    comc.waitOnRequest();
    //    cout<<RED<<Com.myrank<<" offset "<<offset<<RESET<<endl;

    if ( !zoltanGeometricPartitioner( myvalue, totalvalue, offset, setmethod, zz, XYZ, weight, &zoltan_out ) )
    {
        //    return;
    }

    // after repartition, move them around here
    // need to tranfer seed as well as trees in one message;
    // need to preserve the first element to specify the number of seedsi, sizeof tree and tree elements, sizeof nexte tree, tree elements
    // ,...
    // and so on and so forth
    //

    std::string filename = "zoltan";
    filename.append( to_string( Com.myrank ) );
    ofstream f;
    f.open( filename );

    // sort which procesors are the destination here

    uint countdest;
    vector<uint> dest;
    for ( uint i = 0; i < zoltan_out.numExport; i++ )
    {
        if ( std::find( dest.begin(), dest.end(), zoltan_out.exportProcs[i] ) == dest.end() )
        {
            dest.push_back( zoltan_out.exportProcs[i] );
        }
    }
    f << "number of trees to export" << zoltan_out.numExport << endl;

    for ( uint i = 0; i < dest.size(); i++ )
    {
        f << " dest " << dest.at( i ) << endl;
    }

    // which seeds to send,
    bitvector<M> message[dest.size()];
    auto it2 = seeds.begin();
    for ( uint i = 0; i < dest.size(); i++ )
    {
        for ( uint j = 0; j < zoltan_out.numExport; j++ )
        {
            if ( zoltan_out.exportProcs[i] == dest.at( i ) )
            {
                it2 = std::next( seeds.begin(), zoltan_out.exportLocalGids[j] );
                message[i].push_back( ( *it2 ) );
            }
        }
    }

    for ( uint i = 0; i < dest.size(); i++ )
    {
        for ( uint j = 0; j < message[i].size(); j++ )
        {
            cout << " seeds " << message[i].at( j ) << endl;
        }
    }

    //
    // sort the tree in there
    //

    bool *sendbuff[dest.size()];
    uint msize[dest.size()];
    uint idx;

    for ( uint i = 0; i < dest.size(); i++ )
    {
        msize[i] = 0;
        for ( uint j = 0; j < message[i].size(); j++ )
        {
            auto it7 = std::find( seeds.begin(), seeds.end(), message[i].at( j ) );
            idx = std::distance( seeds.begin(), it7 );
            auto it3 = std::next( trees.begin(), idx );
            msize[i] += ( *it3 ).size();
        }
        // one is added to tell the number of trees for the reciever
        sendbuff[i] = new bool[( msize[i] + 2 * message[i].size() + 1 ) * N];
        memset( sendbuff[i], 0, sizeof( bool ) * ( msize[i] + 2 * message[i].size() + 1 ) * N );
        f << "size of total tree elements message " << msize[i] << endl;
    }

#if ( 1 )
    // assign the trees to the buffer
    uint totalsize;
    morton<N> kt;
    co = 0;
    morton<M> kt2;
    MPI_Request request[dest.size()];

    MPI_Status status;
    for ( uint i = 0; i < dest.size(); i++ )
    {
        // seed morton code and size of the following tree makes it 2*message[i]
        totalsize = msize[i] + 2 * message[i].size() + 1;
        morton<N> key( message[i].size() );
        // get the entire size
        co = 0;

        for ( uint j = 0; j < N; j++ )
        {
            sendbuff[i][co * N + j] = key[j];
        }
        co++;

        for ( uint j = 0; j < message[i].size(); j++ )
        {
            auto it7 = std::find( seeds.begin(), seeds.end(), message[i].at( j ) );

            kt2 = message[i].at( j );
            // shift to fit M bits in N bits N>M

            for ( uint k = 0; k < M; k++ )
            {
                sendbuff[i][N * co + ( N - M ) + k - 1] = kt2[k];
            }
            co++;

            idx = std::distance( seeds.begin(), it7 );

            if ( it7 == seeds.end() )
            {
                throw std::runtime_error( RED "seed not in the list in function geomrepart" RESET );
            }

            auto it3 = std::next( trees.begin(), idx );

            morton<N> key1( ( *it3 ).size() );
            //           f<<"seed code " << message[i].at(j)<<endl;
            //           f<<"tree size sent "<<(*it3).size()<<endl;
            //           f<<" size morton code"<<key1<<endl;
            for ( uint k = 0; k < N; k++ )
            {
                sendbuff[i][N * co + k] = key1[k];
            }
            co++;
            for ( uint k = 0; k < ( *it3 ).size(); k++ )
            {
                ( *it3 ).getKey( k, key );
                for ( uint l = 0; l < N; l++ )
                {
                    sendbuff[i][N * co + l] = key[l];
                }
                co++;
            }
        }
        //   cout<<">>>>>>>>>>>>>>"<<dest.at(i)<<endl;
        //       f<<"bytes sent "<<totalsize* N * sizeof( bool )<<endl;
        //       f<<" number of elements sent "<<totalsize<<endl;
        MPI_Isend( sendbuff[i], totalsize * N * sizeof( bool ), MPI_BYTE, dest.at( i ), Com.myrank, MPI_COMM_WORLD, &request[i] );
    }

#endif

    // delete the seed and the corresponding tree for the sender

    // f<<Com.myrank<<" seedsize "<<seeds.size()<<" treesize "<<trees.size()<<" forest size "<<forestsize()<<endl;

    for ( uint i = 0; i < dest.size(); i++ )
    {
        for ( uint j = 0; j < message[i].size(); j++ )
        {
            auto it7 = std::find( seeds.begin(), seeds.end(), message[i].at( j ) );

            std::find( seeds.begin(), seeds.end(), message[i].at( j ) );

            uint idx = std::distance( seeds.begin(), it7 );

            if ( it7 == seeds.end() )
            {
                throw std::runtime_error( RED "invalid cell" RESET );
            }

            auto it3 = std::next( trees.begin(), idx );

            seeds.remove( message[i].at( j ) );
            trees.erase( ( it3 ) );
        }
    }

    // update proc topology now that it is modified

    // now probe for the message to be recieved
    vector<uint> src;
    for ( uint i = 0; i < zoltan_out.numImport; i++ )
    {
        if ( std::find( src.begin(), src.end(), zoltan_out.importProcs[i] ) == src.end() )
        {
            src.push_back( zoltan_out.importProcs[i] );
        }
    }

    bool *recvbuff[src.size()];
    // probe the message
    // recieve the message, if size bigger than zero then Irecv the message

    // cout<<YELLOW<<Com.myrank<<" src size "<<src.size()<<"dest size"<<dest.size()<<" "<<zoltan_out.numImport<<RESET<<endl;
    // f<<Com.myrank<<" src size "<<src.size()<<" dest size "<<dest.size()<<" import size  "<<zoltan_out.numImport<<endl;
    int size;
    uint recvtotalsize = 0;
    MPI_Request request2[src.size()];

    uint co2 = 0;
    real X[6];
    real len[3], XC[3];
    uint treesize;
    uint nx = 2, ny = 2, nz = 2;

    for ( uint i = 0; i < src.size(); i++ )
    {
        //  cout << Com.myrank << "waiting" << endl;
        MPI_Probe( src.at( i ), src.at( i ), MPI_COMM_WORLD, &status );
        MPI_Get_count( &status, MPI_BYTE, &size );
        recvbuff[i] = new bool[size / sizeof( bool )];
        memset( recvbuff[i], 0, size );
        //  f<<" revd size "<<size<<" nelem "<<size/sizeof(bool)<<endl;
        MPI_Irecv( recvbuff[i], size, MPI_BYTE, src.at( i ), src.at( i ), MPI_COMM_WORLD, &request2[i] );
        MPI_Wait( &request2[i], &status );
//    MPI_Recv( recvbuff[i], size, MPI_BYTE, src.at( i ), src.at( i ), MPI_COMM_WORLD, &status );

// debug to see if recvd message is the same
#if ( 1 )
        if ( size != 0 )
        {
            co2 = 0;
            for ( uint j = 0; j < N; j++ )
            {
                kt[j] = recvbuff[i][co2 * N + j];
            }
            co2++;

            //   f<<"size of seeds"<<kt<<endl;

            // get the number of trees
            recvtotalsize = kt.to_ulong();
            //   f<<" converted recieved size "<<recvtotalsize<<endl;

            for ( uint l = 0; l < recvtotalsize; l++ )
            { // get the seed

                for ( uint k = 0; k < M; k++ )
                {
                    kt2[k] = recvbuff[i][co2 * N + ( N - M ) + k - 1];
                }
                co2++;
                // insert this to the list of seeds
                // f<<"seed code "<<kt2<<endl;

                seeds.push_back( kt2 );
                // insert a  new tree to the list
                proc.enclosingBox( kt2, X );

                XC[0] = ( X[0] + X[1] ) * 0.5;
                XC[1] = ( X[2] + X[3] ) * 0.5;
                XC[2] = ( X[4] + X[5] ) * 0.5;

                len[0] = fabs( X[1] - X[0] );
                len[1] = fabs( X[3] - X[2] );
                len[2] = fabs( X[5] - X[4] );

                trees.push_back( Tree<N, Nvalue>( len, XC, nx, ny, nz ) );

                // get the size
                for ( uint j = 0; j < N; j++ )
                {
                    kt[j] = recvbuff[i][co2 * N + j];
                }
                co2++;
                treesize = kt.to_ulong();
                // f<<kt<<endl;
                // f<<treesize<<endl;
                for ( uint j = 0; j < treesize; j++ )
                {
                    for ( uint k = 0; k < N; k++ )
                    {
                        kt[k] = recvbuff[i][co2 * N + k];
                    }
                    // add to mesh list
                    ( trees.back() ).insertKey( kt );
                    // f<<co2<<endl;
                    co2++;
                }
                //  f<<co2<<endl;
                // insert this tree
                //  loop over trees
            }
        }
#endif
    }

#if ( 1 )
    //    if ( Com.myrank == 0 )
    {
        printf( RED "to be imported by rank(%d) %d\n" RESET, Com.myrank, zoltan_out.numImport );
        for ( int i = 0; i < zoltan_out.numExport; i++ )
        {
            printf( " %d will send %d   %d  max size %d \n", Com.myrank, zoltan_out.exportLocalGids[i], zoltan_out.exportProcs[i],
                    seeds.size() );
            f << Com.myrank << " "
              << "will send " << zoltan_out.exportLocalGids[i] << " " << zoltan_out.exportProcs[i] << endl;
        }
        cout << "===============================================================" << endl;
        for ( int i = 0; i < zoltan_out.numImport; i++ )
        {
            printf( "%d will recieve  %d  %d  max size %d\n", Com.myrank, zoltan_out.importLocalGids[i], zoltan_out.importProcs[i],
                    seeds.size() );
            f << Com.myrank << " "
              << "will recvie" << zoltan_out.importLocalGids[i] << " " << zoltan_out.importProcs[i] << endl;
        }
    }
// f<<"source  size "<<src.size()<<endl;
#endif
    Zoltan_LB_Free_Part( &zoltan_out.importGlobalGids, &zoltan_out.importLocalGids, &zoltan_out.importProcs, &zoltan_out.importToPart );

    Zoltan_LB_Free_Part( &zoltan_out.exportGlobalGids, &zoltan_out.exportLocalGids, &zoltan_out.exportProcs, &zoltan_out.exportToPart );
    delete[] weight;

    // watch out for deleting the pointers before making sure it has arrived

    for ( uint i = 0; i < dest.size(); i++ )
    {
        MPI_Wait( &request[i], &status );
        delete[] sendbuff[i];
    }

    for ( uint i = 0; i < src.size(); i++ )
    {
        delete[] recvbuff[i];
    }
    f << Com.myrank << " seedsize " << seeds.size() << " treesize " << trees.size() << " forest size " << forestsize() << endl;
}
//
//
// \brief
//   Move the geometry
//   4:1 balance needs to be retained for keeping the balance
//   need to ask for permission from the neighboring tree (or processor)
//   to keep a boundary element in the derefinement list
//
//
template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void Forest<N, Nvalue, M, Mvalue>::debugDerefine( Tree<M, Mvalue> &proc )
{
#if ( 0 )
    morton<N> key;
    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).addToList( key );
        ( *it ).refineRefineList();
    }

    morton<N> key1 = 0;
    std::unordered_map<morton<N>, int> refinelist;
    key1.flip( N - 1 );
    key1.flip( N - 2 );
    key1.flip( N - 3 );

    auto it1 = seeds.begin();

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        if ( *it1 == 0 )
        {
            ( *it ).addToList( key1 );
            ( *it ).refineRefineList();
        }
        cout << "seed of tree" << ( *it1 ) << endl;
        it1 = std::next( it1, 1 );
    }

    morton<M> key2 = 0;
    key2.flip( M - 1 );
    cout << "seed " << key2 << endl;

    it1 = seeds.begin();

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        if ( *it1 == key2 )
        {
            cout << "the processor who holds the seed number" << Com.myrank << endl;
            ( *it ).addToDerefineList( 0 );
        }
        it1 = std::next( it1, 1 );
    }

    retainFourToOneBalance( proc );

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).derefineDerefineList();
    }

#endif
}
// check the true false situation for 2:1 balance with neighbors

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
bool Forest
<N, Nvalue, M, Mvalue>::checkWithNbrs( bool *sendbuf, bool *recvbuf ) /*!< creates distributed (acalable) graph for communication */
{
    bool bol = false;

    MPI_Neighbor_allgather( sendbuf, sizeof( bool ), MPI_BYTE, recvbuf, sizeof( bool ), MPI_BYTE, graphComm );

    /*
    for(uint i=0;i<destination.size(); i++)
    {
      cout<<YELLOW<<"myrank "<<Com.myrank<<" rcvbuf "<<recvbuf[i]<<RESET<<endl;
    }
    */

    for ( uint i = 0; i < destination.size(); i++ )
    {
        if ( recvbuf[i] == true )
        {
            cout << recvbuf[i] << endl;
            bol = true;
            break;
        }
    }
    return ( bol );
}

// this is to send the destinations the sizes of the messages to be recieved
// if the size == 0, there will be no send and no recieves
// This will reove the necessity to send 0 size messages like I have done in
// the previous 4:1 message at the expense of one message

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void Forest
<N, Nvalue, M, Mvalue>::rcvrMessageSize( int *sendbuf, int *recvbuf ) /*!< creates distributed (acalable) graph for communication */
{

    MPI_Neighbor_alltoall( sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, graphComm );
    /*
        for(uint i=0;i<destination.size(); i++)
        {
          cout<<MAGENTA<<"myrank "<<Com.myrank<<" rcvbuf "<<"["<<i<<"] "<<destination.at(i)<<" "<<recvbuf[i]<<RESET<<endl;
        }
      */
}
//
// nbrsOfNbrs consistency
//

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void Forest<N, Nvalue, M, Mvalue>::checkNbrsOfNbrsConsistency() /*!< creates distributed (acalable) graph for communication */
{

    int size;
    MPI_Status status, status1;
    MPI_Request request[nbrsOfNbrs.size()], request1;
    uint dest[nbrsOfNbrs.size()];

    for ( uint i = 0; i < nbrsOfNbrs.size(); i++ )
    {
        dest[i] = nbrsOfNbrs.at( i );
    }

    for ( uint i = 0; i < nbrsOfNbrs.size(); i++ )
    {
        MPI_Isend( dest, nbrsOfNbrs.size(), MPI_UNSIGNED, nbrsOfNbrs.at( i ), Com.myrank, MPI_COMM_WORLD, &request1 );
    }

    uint *recvbuff[nbrsOfNbrs.size()];

    uint val;

    for ( uint i = 0; i < nbrsOfNbrs.size(); i++ )
    {
        cout << RED << "myrank " << Com.myrank << " " << nbrsOfNbrs.at( i ) << RESET << endl;
        MPI_Probe( nbrsOfNbrs.at( i ), nbrsOfNbrs.at( i ), MPI_COMM_WORLD, &status );

        MPI_Get_count( &status, MPI_UNSIGNED, &size );

        // recvbuff[i] = new morton<M + N>[ size ];
        // notice we revcieved the message by byte so need to specify,  to find the number of elements simly divide it by sizeof
        // bool
        recvbuff[i] = new uint[size];
        //           cout << "******************* " << size << endl;

        MPI_Irecv( recvbuff[i], size, MPI_UNSIGNED, nbrsOfNbrs.at( i ), nbrsOfNbrs.at( i ), MPI_COMM_WORLD, &request[i] );

        MPI_Wait( &request[i], &status1 );

        if ( std::find( recvbuff[i], recvbuff[i] + size, Com.myrank ) == recvbuff[i] + size )
        {
            cout << "inconsistency occured at rank " << Com.myrank << endl;
        }
    }

    if ( Com.myrank == 0 )
    {
        cout << GREEN << "**************************************************************\n" << RESET << endl;

        cout << MAGENTA << "\t Graph Consistency Check: Successful \n" << RESET << endl;

        cout << GREEN << "**************************************************************\n" << RESET << endl;
    }
    // cout << YELLOW << Com.myrank << " = " << nbrsOfNbrs.size() << RESET << endl;

    // cout << YELLOW << Com.myrank << " = " << destination.size() << RESET << endl;
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue>
void Forest
<N, Nvalue, M, Mvalue>::checkZoltanPartConsistency( Tree<M, Mvalue> &proc ) /*!< creates distributed (acalable) graph for communication */
{

    uint numtrees = 0;

    for ( auto it = proc.begin(); it != proc.end(); it++ )
    {
        if ( Com.myrank == it->second[0] )
        {
            numtrees++;
        }
    }

    //    cout<<"num Tree "<<my_rank<<'\t'<<numtrees<<endl;

    if ( numtrees == 0 )
    {
        throw std::runtime_error( RED "Invalid partitioning: every Processor should at least have one cube" RESET );
    }

    uint total_trees;

    MPI_Reduce( &numtrees, &total_trees, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );

    if ( Com.myrank == 0 )
    {
        cout << RED << "tree size after partitioning =  " << total_trees << RESET << endl;
        if ( proc.size() != total_trees )
        {

            throw std::runtime_error( "inconsistent partitioning by Zoltan" );
        }
    }

    if ( proc.size() < Com.comsize )
    {
        throw std::runtime_error( RED "invalid topology " RESET );
    }
}

#endif

// this basically tells wich procesors are active for the while loop
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::nonCollectiveNbrComm()
{

    uint *sendbuf = new uint[destination.size()];
    uint *recvbuf = new uint[destination.size()];
    ;

    for ( uint i = 0; i < destination.size(); i++ )
    {
        sendbuf[i] = 0;
        recvbuf[i] = 0;
        sendtag.push_back( 1 );
        recvtag.push_back( 1 );
        cout << Com.myrank << "  " << destination.at( i ) << endl;
    }

    MPI_Request request, request0;
    MPI_Status status;

    if ( Com.myrank == 0 )
    {
        sendbuf[Com.myrank + 0] = 1;

        sendbuf[Com.myrank + 1] = 1;
    }

    uint Scount = 0;

    uint Rcount = 0;

    uint tmp;

    for ( uint i = 0; i < destination.size(); i++ )
    {
        if ( sendtag[i] == 1 )
        {
            MPI_Isend( sendbuf + i, 1, MPI_UNSIGNED, destination.at( i ), Com.myrank, Com.mpicom, &request0 );
        }
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {
        if ( recvtag[i] == 1 )
        {
            MPI_Irecv( recvbuf + i, 1, MPI_UNSIGNED, destination.at( i ), destination.at( i ), Com.mpicom, &request );

            MPI_Wait( &request, &status );
        }
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {

        cout << RED << Com.myrank << " " << recvbuf[i] << RESET << endl;
    }

    // shift send and recieve now

    for ( uint i = 0; i < destination.size(); i++ )
    {
        if ( sendbuf[i] == 1 )
        {
            sendbuf[i] = 0;
            recvbuf[i] = 1;
        }
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {
        if ( sendbuf[i] == 0 && recvbuf[i] == 1 )
        {
            recvbuf[i] = 1;
        }
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {

        cout << GREEN << Com.myrank << " " << sendbuf[i] << RESET << endl;
    }

    // one-sidedly communicte this to every neighbor
    /*
     uint max_buff_size=6;
     uint *rma_buff;
    MPI_Win win;

     MPI_Alloc_mem(sizeof(uint)*max_buff_size, MPI_INFO_NULL,&rma_buff);

     MPI_Win_create(rma_buff,max_buff_size*sizeof(uint),sizeof(uint),MPI_INFO_NULL,Com.mpicom,&win);

     MPI_Win_fence(0, win);

   // add my_rank communicate and subtract to get the rank of the sender

   for(i=0;i<destination.size();i++)
   {
   if(sendbuff[i]==1)
   {
     MPI_Put(rma_buff+i,my_buff_size,MPI_UNSIGNED,destination,displacement2,my_buff_size,MPI_UNSIGNED,win);
   }

   }

     MPI_Win_fence(0, win);
   */

    /*
    //
    // recommunicate sen/recv tags
    //
        for ( uint i = 0; i < destination.size(); i++ )
        {
            if ( sendtag[i] == 1 )
            {
                MPI_Isend( sendbuf+i, 1, MPI_UNSIGNED,destination.at( i ), Com.myrank, Com.mpicom, &request0 );
            }
        }

       for ( uint i = 0; i < destination.size(); i++ )
        {
            if ( recvtag[i] == 1 )
            {
               MPI_Irecv( recvbuf + i, 1, MPI_UNSIGNED, destination.at( i ), destination.at( i ), Com.mpicom, &request );

               MPI_Wait(&request,&status);

            }
        }

    //
    */

    //  cout<<Scount<<" "<<Rcount<<endl;

    //
}
//
//
// derefinemnt functions
//
//

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::pushToDerefineEachTree( uint nlevel, Tree<M, uint> &proc )
{
    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).pushToDerefinelist( nlevel );

    }

    retainFourToOneBalance( proc );

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).derefineDerefineList();
    }
}

// this version reuires use of the parallel piped enclosing the geometry 
// and the number of inactive geometries
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::pushToDerefineEachTree( uint nlevel, Tree<M, uint> &proc,int nInactive, real *enClosingBoxForInactiveGeoms )
{
    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).pushToDerefinelist( nlevel );

        //( *it ).ignoreInactive(nInactive, enClosingBoxForInactiveGeoms);       

        ( *it ).ignoreInactiveVertices(nInactive, enClosingBoxForInactiveGeoms);       

    }

    retainFourToOneBalance( proc );

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).derefineDerefineList();
    }
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::retainFourToOneBalance( Tree<M, uint> &proc )
{

    const uint NM = N + M;
    morton<NM> combinedkey = 0, ktcom;
    morton<N> key, kt1, nbrkey;
    morton<M> seedkey, seednbrkey = 0, kt, kt3;
    auto it3 = seeds.begin();
    uint index;
    bool bol;
    real xyz[6];
    bitvector<N> boundaryElem;
    uint mylevel, changedirectionlevel, direction;
    uint topologylevel, nbrseedlevel, nbrlevel, complevel, nbrcomplevel, nbrtopologylevel, combinedlevel, localcombinedlevel;
    uint counter;
    vector<uint> directions;
    it3 = seeds.begin();
    vector<bitset<M + N>> message[destination.size()];
    vector<bitset<M + N>> message3[destination.size()];
    //    vector<bitset<M + N>> message2[destination.size()];
    vector<uint> message2[destination.size()];
    vector<uint> dest;
    vector<uint> source;
    uint idx;
    int size;
    int flag = 1;
    morton<M + N> rcvkey;
    uint seedlevel;
    morton<N> elementkey;
    uint idex2;
    uint elementlevel, recvmessagecombinedlevel;
    vector<uint> start;
    vector<uint> end;
    uint istart, iend;
    uint counter3, counter2, counter5;
    MPI_Request request[destination.size()], request1[destination.size()], request0;
    MPI_Status status;
    morton<N + M> tempkey;

    std::string filename = "rank0";
    filename.append( to_string( Com.myrank ) );
    ofstream myfile;
    myfile.open( filename );

    uint con = 0;
    char *sendbuff[destination.size()];
    char *recvbuff[destination.size()];
    uint *sendbuff2[destination.size()];
    uint *recvbuff2[destination.size()];

    // extract the boundary elements for each tree from the tagged list
    // get the ectent that we need to refine in the initia step
    // we know the number of trees therefore, initialize it here
    // myfile << "Writing this to a file.\n";
    // the main while loop

    cout << "size start=" << start.size() << endl;
    cout << "destination ize " << destination.size() << endl;
    cout << "size end=" << end.size() << endl;

    bool loopbool = true;
    uint bcstart, bcend;
    //
    // This algorithm is a bit different :
    // if the element tagged to be removed destroys the 4:1 balance
    // remove it from the tagged list, no ripple effect in this situation is possible
    // while ( loopbool ), basically each tree askes for permission to remove the element
    //

    // start out by setting loopbool=0
    counter2 = 0;
    //
    //
    //      separates the list between local boundary and non-local boundary
    //
    //

    it3 = seeds.begin();

    counter5 = 0;

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        ( *it ).retainFourToOne();
    }

    for ( auto it = trees.begin(); it != trees.end(); it++ )
    {
        seedkey = ( *it3 );
        proc.level( seedkey, &topologylevel );

        for ( auto i = ( *it ).Dbegin(); i != ( *it ).Dend(); i++ )
        {
            combinedkey = 0;

            for ( uint j = 0; j < 3 * topologylevel; j++ )
            {
                combinedkey[NM - j - 1] = seedkey[M - j - 1];
            }

            // gives the key in the list
            //

            key = ( *it ).readDerefineList( i );

            ( *it ).level( key, &mylevel );

            //     cout << "key " << key << " mylevel " << mylevel << endl;

            combinedlevel = topologylevel + mylevel;

            for ( uint j = 0; j < 3 * mylevel; j++ )
            {
                combinedkey[NM - 3 * ( topologylevel ) - j - 1] = key[N - j - 1];
            }

            ktcom = combinedkey;
//  cout<<" combined key "<<combinedkey<<endl;
#if ( 1 )
            for ( uint direction = 0; direction < 3; direction++ )
            {
                if ( !( *it ).isBoundary( key, direction ) )
                {
                    continue;
                }
                // cout << RED" direction " << direction << endl;
                // cout << RED" combined key "RESET << combinedkey << endl;
                findFlipLevel( combinedkey, &combinedlevel, &changedirectionlevel, &direction );

                if ( changedirectionlevel != 0 )
                {
                    flipForNbr( combinedkey, &combinedlevel, &changedirectionlevel, &direction );
                    //                        cout << " combined key " << combinedkey << " combined level " << combinedlevel << endl;
                    // toplogy level stinks, bug fixed
                    //                   getNbrSeedLevel( combinedkey, topologylevel, &nbrseedlevel, proc );
                    getNbrSeedLevel( combinedkey, maxProcLevel, &nbrseedlevel, proc );

                    seednbrkey = 0;

                    for ( uint k = 0; k < 3 * nbrseedlevel; k++ )
                    {
                        seednbrkey[M - 1 - k] = combinedkey[NM - k - 1];
                        // cout<<RED<<combinedkey[NM-k-1]<<RESET<<endl;
                    }

                    // cout << "seednbrkey " << BLUE << seednbrkey << RESET << endl;
                    // myfile << "  seed key  "<<seednbrkey <<  endl;
                    nbrkey = 0;

                    for ( uint k = 0; k < N; k++ )
                    {
                        nbrkey[N - k - 1] = combinedkey[NM - 3 * ( nbrseedlevel ) - 1 - k];
                    }

                    // cout << " nbrkey " << nbrkey << endl;

                    if ( isInSeed( seednbrkey, &counter ) )
                    {
                        auto it2 = std::next( trees.begin(), counter );
                        auto it4 = ( *it2 ).find( nbrkey );

                        // seedlevel
                        // ( proc ).level( seednbrkey, &nbrlevel );
                        // initial guess on nbr level as constructed
                        // this is a singularity for level zero, it might be negative for first level
                        // I remove this mannually

                        if ( combinedlevel >= nbrseedlevel )
                        {
                            nbrlevel = combinedlevel - nbrseedlevel;
                        }
                        else
                        {
                            nbrlevel = 1;
                        }
                        //      cout<<RED<<nbrlevel<<RESET<<endl;

                        // now find the real levels

                        if ( it4 != ( *it2 ).end() )
                        {
                            ( *it2 ).level( it4->first, &nbrlevel );
                            //    cout << YELLOW << nbrlevel << RESET << endl;
                        }
                        else
                        {
                            // this condition imples level is lower no need for modification and search, I commented it out
                            //                                nbrlevel                     = mylevel-1;
                            // ( *it2 ).level( it4->first, &nbrlevel );
                            nbrlevel = nbrlevel - 1;
                            nbrkey[N - 3 * nbrlevel - 1] = 0;
                            nbrkey[N - 3 * nbrlevel - 2] = 0;
                            nbrkey[N - 3 * nbrlevel - 3] = 0;

                            //            cout <<RED<< N - 3 * nbrlevel - 1 << " " << mylevel << RESET<<endl;
                            //                                      auto it4=(*it2).find(nbrkey);
                            //                                      (*it).level(it4->first,&nbrlevel);
                        }

                        // cout << "nbrkey " << RED << nbrkey << "  count " << counter << RESET << endl;
                        // cout<<"inside"<<endl;
                        if ( ( *it2 ).find( nbrkey ) == ( *it2 ).end() )
                        {
                            cout << GREEN << nbrkey << "direction " << direction << "level " << nbrlevel << RESET << endl;
                            throw std::runtime_error( "error in finding key in derefine" RESET );
                        }

                        nbrcomplevel = nbrlevel + nbrseedlevel;

                        //                         cout<<RED<<nbrlevel<<" "<<nbrseedlevel<<RESET<<endl;
                        cout << YELLOW << nbrcomplevel << " " << combinedlevel << RESET << endl;
                        // if not in the list add to the list, dont forget to modify here
                        // watch out this element might have already been tagged

                        if ( nbrcomplevel > combinedlevel )
                        {
                            cout << " inside here" << key << endl;
                            //                                   remove from list of derefinement actually
                            ( *it ).removeFromDerefineList( i );
                            // removes
                            //                            counter5=counter5-8;
                            //
                            //

                            // here is where boolwhile is affected
                            cout << RED "======================================" << endl;
                            cout << "removed from list" << nbrkey << endl;
                            cout << "======================================" RESET << endl;
                        }
                        // cout<<seednbrkey<<RED<<seednbrkey<<RESET<<endl;
                    }

                    else
                    {
                        // sorts for communication

                        // prepare for communication, pack all the elements with the destination, sender is   my_rank

                        // find seednbrkey from global data
                        auto it5 = proc.find( seednbrkey );
                        //                        dest.push_back( it5->second[0] );
                        auto it6 = destination.begin();

                        if ( it5 != proc.end() )
                        {
                            it6 = find( destination.begin(), destination.end(), it5->second[0] );
                        }

                        // combined key should be avoided
                        idx = it6 - destination.begin();
                        // cout << "---------------------->>> index" << RED << idx << RESET << endl;
                        // before pushing back, need to be careful about the key with all zeros,
                        // to accommodate this, the first bit (far right) is flipped and a sibling is sent

                        removeAllZeroSingularity( combinedkey, combinedlevel );
                        //                             cout << GREEN << " proc_id " << Com.myrank  << "combinedkey  " << combinedkey;
                        //     << "combined level " << combinedlevel << endl;

                        message[idx].push_back( combinedkey );
                        message3[idx].push_back( ktcom );

                        // myfile<< " combined key "<< combinedkey <<endl;
                        //                        cout << message.at( 0 ) << RESET << endl;
                    }
                }
                combinedkey = ktcom;
            }

#endif
        }
        //       counter5++;
        it3 = std::next( it3, 1 );
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {
        // cout << RED << "my Rank " << Com.myrank << "message size" << message[i].size() << RESET << endl;
        myfile << "Rank " << Com.myrank << "sends " << message[i].size() << "size message to " << destination.at( i ) << endl;
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {
        for ( uint j = 0; j < message[i].size(); j++ )
            // cout << RED << "my Rank " << Com.myrank << " " << message[i].at( j ) << " " << destination.at( i ) << RESET << endl;
            myfile << "my Rank " << Com.myrank << " " << message[i].at( j ) << " " << destination.at( i ) << endl;
    }

    // need a communication pattern, need to send zero size message if communication is not required  need to add local balance criteria
    // need to
    // loop over and do this again in while loop
    // get size of message for each proc send the required information to the neighboring processes

    morton<N + M> rcvkey2;

#if ( 1 )

    for ( uint i = 0; i < destination.size(); i++ )
    {
        sendbuff[i] = new char[message[i].size() * ( M + N )];

        for ( uint j = 0; j < message[i].size(); j++ )
        {
            tempkey = message[i].at( j );

            for ( uint k = 0; k < M + N; k++ )
            {
                if ( tempkey[k] == true )
                {
                    sendbuff[i][j * ( M + N ) + k] = '1';
                }
                else
                {
                    sendbuff[i][j * ( M + N ) + k] = '0';
                }
            }
        }

        //  cout << "================================== " << endl;
        //  cout << RESET "messagesize " << message[i].size() << " destination " << destination.at( i ) << RESET << endl;
        //  cout << "================================== " << endl;
        // send messages to destinations and tag each meaasge with self rank (myrank)
        //
        MPI_Isend( sendbuff[i], message[i].size() * ( M + N ), MPI_CHAR, destination.at( i ), Com.myrank, MPI_COMM_WORLD, &request1[i] );
    }

    // even with overlapping communication and computation it is not garanteed that this part
    // will use enought time for the message to arrive, since this list might be empty
    // however we can use blocking after this section, that we know not all of the time is spent waiting
    for ( uint i = 0; i < destination.size(); i++ )
    {
        //             MPI_Iprobe(destination.at(i), destination.at(i), MPI_COMM_WORLD,&flag, &status);
        // this is a symmetric comm pattern therefore any
        //
        MPI_Probe( destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &status );

        // while(flag!=0)
        {
            // MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&flag, &status);
            if ( flag == 1 )
            {
                MPI_Get_count( &status, MPI_CHAR, &size );

                // recvbuff[i] = new morton<M + N>[ size ];
                // notice we revcieved the message by byte so need to specify,  to find the number of elements simly divide it by sizeof
                // bool
                recvbuff[i] = new char[size];
                //           cout << "******************* " << size << endl;

                MPI_Irecv( recvbuff[i], size, MPI_CHAR, destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &request[i] );

                MPI_Wait( &request[i], &status );

                // if ( size != 0 )
                {
                    for ( uint j = 0; j < size / ( M + N ); j++ )
                    {
                        rcvkey = 0;

                        for ( uint k = 0; k < N + M; k++ )
                        {
                            if ( recvbuff[i][j * ( N + M ) + k] == '1' )
                            {
                                rcvkey.flip( k );
                            }
                        }

                        rcvkey2 = rcvkey;
                        combinedLevel( rcvkey, &recvmessagecombinedlevel );

                        recoverAllZeroSingularity( rcvkey, recvmessagecombinedlevel );

                        findSeedLevelForRcvdMessage( rcvkey, &seedlevel, proc );

                        // cout << rcvkey << endl;
                        //  cout << "seedlevel " << seedlevel << endl;

                        constructSeedKeyForRcvdMessage( rcvkey, seedlevel, seedkey );
                        //  cout << "seedkey " << seedkey << endl;

                        /*      	    uint check;
                                           proc.level(seedkey,&check);
                                           cout<<RED<<check<<" "<<seedlevel<<RESET<<endl;
                                           if(seedlevel!=check)
                                          {
                                           throw std::runtime_error("levels calculated for seed inconsistent");
                                           }
                          */
                        constructElementKeyForRcvdMessage( rcvkey, seedlevel, elementkey );
                        //            cout << "element key " << elementkey << endl;

                        // cout << RED "rcvkey level " << recvmessagecombinedlevel << RESET << endl;

                        auto it7 = std::find( seeds.begin(), seeds.end(), seedkey );

                        if ( it7 == seeds.end() )
                        {
                            cout << RED << "myrank " << Com.myrank << " " << seedkey << RESET << endl;

                            throw std::runtime_error( RED "seed not found in Balance Comm" RESET );
                        }

                        idex2 = std::distance( seeds.begin(), it7 );

                        auto it9 = std::next( trees.begin(), idex2 );

                        if ( recvmessagecombinedlevel >= seedlevel )
                        {
                            elementlevel = recvmessagecombinedlevel - seedlevel;
                        }
                        else
                        {
                            elementlevel = 1;
                        }
                        //   myfile << rcvkey << " " << recvmessagecombinedlevel << endl;

                        //  cout << RED << Com.myrank << rcvkey << " " << recvmessagecombinedlevel << RESET << endl;
                        //    myfile << elementkey << endl;
                        //    myfile << seedkey << " " << seedlevel << endl;
                        //    myfile << elementlevel << endl;
                        ( *it9 ).level( elementkey, &elementlevel );
                        //    myfile << elementlevel << endl;

                        if ( ( *it9 ).isInMeshList( elementkey ) == false )
                        {
                            //   cout << "inside mesh" << elementkey << endl;
                            elementkey[N - 3 * ( elementlevel - 1 ) - 1] = 0;

                            elementkey[N - 3 * ( elementlevel - 1 ) - 2] = 0;

                            elementkey[N - 3 * ( elementlevel - 1 ) - 3] = 0;
                        }
                        //     myfile<<elementkey<<endl;
                        ( *it9 ).level( elementkey, &elementlevel );

                        localcombinedlevel = elementlevel + seedlevel;

                        // add if the level is higher put this in list of rejected derefinement list and sent it back

                        //   cout << BLUE << localcombinedlevel << "            " << recvmessagecombinedlevel << RESET << endl;
                        if ( localcombinedlevel > recvmessagecombinedlevel )
                        {
                            //     cout << BLUE << localcombinedlevel << "            " << recvmessagecombinedlevel << RESET << endl;

                            //               cout << "elementkey" << elementkey << endl;

                            // send the recvkey back implying that derefinement is rejected
                            //  puts in destination index
                            // get the index of the source and put it in appropriate location in message2
                            idx = std::find( destination.begin(), destination.end(), destination.at( i ) ) - destination.begin();

                            //   message2[idx].push_back( rcvkey2 );
                            message2[idx].push_back( j );
                        }
                    }
                }
            }
        }
    }

#endif

    // send your response from the recieved message

    for ( uint i = 0; i < destination.size(); i++ )
    {
        //      cout << YELLOW << "my Rank " << Com.myrank << "message size" << message2[i].size() << RESET << endl;
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {
        for ( uint j = 0; j < message2[i].size(); j++ )
        {
            cout << YELLOW << "my Rank " << Com.myrank << " " << message2[i].at( j ) << " " << destination.at( i ) << RESET << endl;
        }
    }

// send it back like a tennis ball
#if ( 1 )
    for ( uint i = 0; i < destination.size(); i++ )
    {
        sendbuff2[i] = new uint[message2[i].size()];

        for ( uint j = 0; j < message2[i].size(); j++ )
        {

            sendbuff2[i][j] = message2[i].at( j );
        }

        // send messages to destinations and tag each meaasge with self rank (myrank)
        //
        MPI_Isend( sendbuff2[i], message2[i].size(), MPI_UNSIGNED, destination.at( i ), Com.myrank, MPI_COMM_WORLD, &request1[i] );
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {
        //
        MPI_Probe( destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &status );
        MPI_Get_count( &status, MPI_UNSIGNED, &size );

        // recvbuff[i] = new morton<M + N>[ size ];
        // notice we revcieved the message by byte so need to specify,  to find the number of elements simly divide it by sizeof
        // bool
        recvbuff2[i] = new uint[size];
        //           cout << "******************* " << size << endl;
        MPI_Irecv( recvbuff2[i], size, MPI_UNSIGNED, destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &request[i] );

        MPI_Wait( &request[i], &status );

        for ( uint j = 0; j < size; j++ )
        {
            rcvkey = message3[i].at( recvbuff2[i][j] );

            myfile << " rcvkey " << rcvkey << endl;

            getNbrSeedLevel( rcvkey, maxProcLevel, &seedlevel, proc );

            seedkey = 0;

            for ( uint k = 0; k < 3 * seedlevel; k++ )
            {
                seedkey[M - 1 - k] = rcvkey[NM - k - 1];
                // cout<<RED<<combinedkey[NM-k-1]<<RESET<<endl;
            }

            for ( uint k = 0; k < N; k++ )

            {
                elementkey[N - k - 1] = rcvkey[NM - 3 * ( seedlevel ) - 1 - k];
            }

            auto it7 = std::find( seeds.begin(), seeds.end(), seedkey );

            idex2 = std::distance( seeds.begin(), it7 );
            auto it9 = std::next( trees.begin(), idex2 );

            if ( it9 == trees.end() )
            {
                cout << " rank " << Com.myrank << " seedkey " << seedkey << endl;
                throw std::runtime_error( "seed not found in derefine operations" );
            }
            else
            {
                myfile << " element " << elementkey << endl;
                ( ( *it9 ).findInDerefine( elementkey ) );

                ( *it9 ).removeFromDerefineList( ( *it9 ).findInDerefine( elementkey ) );
            }
        }
    }
// remove from the list after recieve

#endif

    for ( uint j = 0; j < destination.size(); j++ )
    {
        MPI_Wait( &request[j], &status );
        delete[] sendbuff[j];
        message[j].clear();
        message[j].shrink_to_fit();
        delete[] recvbuff[j];

        MPI_Wait( &request1[j], &status );
        delete[] sendbuff2[j];
        message2[j].clear();
        //                    message2[j].shrink_to_fit();
        message3[j].clear();
        delete[] recvbuff2[j];
    }

    //
    // need to use the recv buffer
    //

    // assigne recieved elements to corresponding lists if needed
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest
<N, Nvalue, M, Mvalue, T>::checkZoltanPartConsistency( Tree
                                                       <M, Mvalue> &proc ) /*!< creates distributed (acalable) graph for communication */
{

    uint numtrees = 0;

    for ( auto it = proc.begin(); it != proc.end(); it++ )
    {
        if ( Com.myrank == it->second[0] )
        {
            numtrees++;
        }
    }

    //    cout<<"num Tree "<<my_rank<<'\t'<<numtrees<<endl;

    if ( numtrees == 0 )
    {
        throw std::runtime_error( RED "Invalid partitioning: every Processor should at least have one cube" RESET );
    }

    uint total_trees;

    MPI_Reduce( &numtrees, &total_trees, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );

    if ( Com.myrank == 0 )
    {
        cout << RED << "tree size after partitioning =  " << total_trees << RESET << endl;
        if ( proc.size() != total_trees )
        {

            throw std::runtime_error( "inconsistent partitioning by Zoltan" );
        }
    }

    if ( proc.size() < Com.comsize )
    {
        throw std::runtime_error( RED "invalid topology " RESET );
    }
}

//****************************************************************
//
//
//       generate graph Communicator (graphComm) comm
//
//
//****************************************************************

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::createCommGraph( uint Lnbr ) /*!< creates distributed (scalable) graph for communication */
{
    //
    // this is a symmetric graph so source and destination ranks are the same and hence,
    // unweighted for now, if Lnbr is set to zero it considers immediate neighbors
    // if set to 1 considers neighbors of neighbors
    //

    uint size;
    uint indegree;
    uint outdegree;
    int *sources = nullptr;

    if ( Lnbr == 0 )
    {
        size = destination.size();
        sources = new int[size];
        indegree = size;
        outdegree = size;

        for ( uint i = 0; i < size; i++ )
        {
            sources[i] = destination[i];

            cout << "myrank " << Com.myrank << " " << destination[i] << endl;
        }
    }

    else if ( Lnbr == 1 )
    {

        size = nbrsOfNbrs.size();
        sources = new int[size];
        indegree = size;
        outdegree = size;

        for ( uint i = 0; i < size; i++ )
        {
            sources[i] = nbrsOfNbrs[i];

            //        cout << "myrank " << Com.myrank << " " << nbrOfNbrs[i] << endl;
        }
    }
    else
    {

        throw std::runtime_error( "Only neighbors and neighbors of neighbors are supported\n" );
    }

    int reorder = REORDER;

    MPI_Dist_graph_create_adjacent( MPI_COMM_WORLD, indegree, sources, MPI_UNWEIGHTED, indegree, sources, MPI_UNWEIGHTED, MPI_INFO_NULL,
                                    reorder, &graphComm );

    delete[] sources;

#if ( DEBUG )
    int status;
    MPI_Topo_test( graphComm, &status );
    if ( status == MPI_DIST_GRAPH )
    {
        // cout<<GREEN<<"status= "<<status<<RESET<<endl;
    }
    else
    {
        throw std::runtime_error( "distributed graph generation failed\n" );
    }
    // int comsize;
    //  MPI_Comm_size( graphComm, &comsize );
    //    cout<<comsize<<endl;
    int incount, outcount;
    // weighted is a flag indicating if the graph is weighted
    int weighted;
    MPI_Dist_graph_neighbors_count( graphComm, &incount, &outcount, &weighted );
    //    cout<<"myrank "<<Com.myrank<<" graph size "<<incount<<endl;
    int maxdegree = size;
    int srcs[maxdegree], dests[maxdegree];
    int sourceweights[maxdegree], destweights[maxdegree];
    MPI_Dist_graph_neighbors( graphComm, maxdegree, srcs, sourceweights, maxdegree, dests, destweights );
    /*
     std::string filename = "rank";
     filename.append( to_string( Com.myrank ) );
     ofstream myfile;
     myfile.open( filename );


    for(uint i=0; i<size ;i++)
    {
        myfile<<"myrank "<<Com.myrank<<" sources "<<srcs[i]<<" destinations "<<dests[i]<<endl;
    }

    myfile.close();
    */

    bool sendbuf = {false};
    bool *recvbuf = nullptr;
    recvbuf = new bool[size];
    // MPI_Neighbor_allgather(&sendbuf, 1 , MPI_UNSIGNED , recvbuf, 1 , MPI_UNSIGNED , graphComm);

    //   checkWithNbrs( &sendbuf, recvbuf );

    for ( uint i = 0; i < size; i++ )
    {
        //  cout<<YELLOW<<"myrank "<<Com.myrank<<" rcvbuf "<<recvbuf[i]<<RESET<<endl;
    }
    delete[] recvbuf;
#endif
}

//=========================================================================================
//
//   need to construct a graph with nbrsOfNbrs as discrete graph due to the ripple effect
//   due to the absorbtion that occure in high levels, this condition is
//   expected to occur at low levels of adaptation, and since we do not perform edge balance
//   precautionarily second degree nbrs are being considered
//
//==========================================================================================

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::createNbrsOfNbrs() /*!< creates distributed (acalable) graph for communication */
{
    /*
        cout << YELLOW << " = " << Com.myrank << '\t' << endl;
        for ( uint i = 0; i < destination.size(); i++ )
        {
            cout << YELLOW << " " << destination.at( i ) << '\t';
        }
        cout << "----------------------- " << RESET << endl;
    */

    uint extendedSize = destination.size() + 1;
    uint extendedDest[extendedSize];

    extendedDest[0] = Com.myrank;
    int size;
    MPI_Status status, status1;
    MPI_Request request[destination.size()], request1;

    for ( uint i = 0; i < destination.size(); i++ )
    {
        nbrsOfNbrs.push_back( destination.at( i ) );
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {
        extendedDest[i + 1] = destination.at( i );
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {
        MPI_Isend( extendedDest, extendedSize, MPI_UNSIGNED, destination.at( i ), Com.myrank, MPI_COMM_WORLD, &request1 );
    }

    uint *recvbuff[destination.size()];

    uint val;

    for ( uint i = 0; i < destination.size(); i++ )
    {
        MPI_Probe( destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &status );

        MPI_Get_count( &status, MPI_UNSIGNED, &size );

        // recvbuff[i] = new morton<M + N>[ size ];
        // notice we revcieved the message by byte so need to specify,  to find the number of elements simly divide it by sizeof
        // bool
        recvbuff[i] = new uint[size];
        //           cout << "******************* " << size << endl;

        MPI_Irecv( recvbuff[i], size, MPI_UNSIGNED, destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &request[i] );

        MPI_Wait( &request[i], &status1 );

        for ( uint j = 0; j < size; j++ )
        {
            val = recvbuff[i][j];
            if ( std::find( nbrsOfNbrs.begin(), nbrsOfNbrs.end(), val ) == nbrsOfNbrs.end() && val != Com.myrank )
            {
                nbrsOfNbrs.push_back( val );
            }
        }
    }

    //    cout << YELLOW << Com.myrank << " nbrsOfNbrs size = " << nbrsOfNbrs.size() << RESET << endl;
    for ( uint i = 0; i < nbrsOfNbrs.size(); i++ )
    {

        cout << YELLOW << Com.myrank << " nbrsOfNbrs size[" << i << "] = " << nbrsOfNbrs.at( i ) << RESET << endl;
    }

    //    cout << YELLOW << Com.myrank << " nbrsSize = " << destination.size() << RESET << endl;
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::checkGraphConsistency() /*!< creates distributed (acalable) graph for communication */
{

    int size;
    MPI_Status status, status1;
    MPI_Request request[destination.size()], request1;
    uint dest[destination.size()];

    for ( uint i = 0; i < destination.size(); i++ )
    {
        dest[i] = destination.at( i );
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {
        MPI_Isend( dest, destination.size(), MPI_UNSIGNED, destination.at( i ), Com.myrank, MPI_COMM_WORLD, &request1 );
    }

    uint *recvbuff[destination.size()];

    uint val;

    for ( uint i = 0; i < destination.size(); i++ )
    {
        cout << RED << "myrank " << Com.myrank << " " << destination.at( i ) << RESET << endl;
        MPI_Probe( destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &status );

        MPI_Get_count( &status, MPI_UNSIGNED, &size );

        // recvbuff[i] = new morton<M + N>[ size ];
        // notice we revcieved the message by byte so need to specify,  to find the number of elements simly divide it by sizeof
        // bool
        recvbuff[i] = new uint[size];
        //           cout << "******************* " << size << endl;

        MPI_Irecv( recvbuff[i], size, MPI_UNSIGNED, destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &request[i] );

        MPI_Wait( &request[i], &status1 );

        if ( std::find( recvbuff[i], recvbuff[i] + size, Com.myrank ) == recvbuff[i] + size )
        {
            cout << "inconsistency occured at rank " << Com.myrank << endl;
        }
    }

    for ( uint i = 0; i < destination.size(); i++ )
    {

        delete[] recvbuff[i];
    }

    if ( Com.myrank == 0 )
    {
        cout << GREEN << "**************************************************************\n" << RESET << endl;

        cout << MAGENTA << "\t Graph Consistency Check: Successful \n" << RESET << endl;

        cout << GREEN << "**************************************************************\n" << RESET << endl;
    }
    // cout << YELLOW << Com.myrank << " = " << nbrsOfNbrs.size() << RESET << endl;
    // cout << YELLOW << Com.myrank << " = " << destination.size() << RESET << endl;
}

// keeping the old version for now, will delete once I am done with it
#if ( 0 )

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N, Nvalue, M, Mvalue, T>::fourToOneBalance( T &proc )
{
    // cout << "character solving" << endl;
    const uint NM = N + M;
    morton<NM> combinedkey = 0, ktcom;
    morton<N> key, kt1, nbrkey;
    morton<M> seedkey, seednbrkey = 0, kt, kt3;
    auto it3 = seeds.begin();
    uint index;
    bool bol;
    real xyz[6];
    bitvector<N> boundaryElem;
    uint mylevel, changedirectionlevel, direction;
    uint topologylevel, nbrseedlevel, nbrlevel, complevel, nbrcomplevel, nbrtopologylevel, combinedlevel, localcombinedlevel;
    uint counter;
    vector<uint> directions;
    it3 = seeds.begin();
    //    vector<bitset<M + N>> message[destination.size()];
    /*vector<vector<bitset<M + N>>> message(destination.size());

    for(int i=0;i<destination.size();i++)
     {
      message[i].push_back();
     }
*/
    vector<bitset<M + N>> *message = new vector<bitset<M + N>>[destination.size()];

    //  message[0].push_back(0);
    //
    vector<uint> dest;
    vector<uint> source;
    uint idx;
    int size;
    int flag = 1;
    morton<M + N> rcvkey;
    uint seedlevel;
    morton<N> elementkey;
    uint idex2;
    uint elementlevel, recvmessagecombinedlevel;
    //    morton<M + N> *       sendbuff[destination.size()];
    //    morton<M + N> *       recvbuff[destination.size()];
    vector<uint> start;
    vector<uint> end;
    uint istart, iend;
    uint counter3, counter2, counter5;
    MPI_Request request[destination.size()], request1[destination.size()], request0;
    MPI_Status status;
    short sb = 0, rb = 0;
    bool mybool;
    morton<N + M> tempkey;
#if ( DEBUG )
    std::string filename = "rank";
    filename.append( to_string( Com.myrank ) );
    ofstream myfile;
    myfile.open( filename );

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
        ////  cout << RED << "prod id " << Com.myrank << "  seed " << ( *it ) << RESET << endl;
        //  myfile << "prod id " << Com.myrank << "  " << con++ << "  seed " << ( *it ) << endl;
    }
    myfile << " topology size " << proc.size() << endl;
    myfile << " destination size " << destination.size() << endl;

    for ( uint i = 0; i < destination.size(); i++ )
    {
        myfile << " destinations " << destination.at( i ) << endl;
    }

#endif

    uint con = 0;
    char *sendbuff[destination.size()];
    char *recvbuff[destination.size()];

    for ( uint i = 0; i < destination.size(); i++ )
    {
        sendbuff[i] = (char *)malloc( 1 * sizeof( char ) );
        recvbuff[i] = (char *)malloc( 1 * sizeof( char ) );
    }

    // the main while loop

    //  cout << "size start=" << start.size() << endl;
    //   cout << "destination ize " << destination.size() << endl;
    //   cout << "size end=" << end.size() << endl;

    bool loopbool = true;
    uint bcstart, bcend;
    bool innerloop = true;

    int intracount = 0;

    bitvector<N> initialList[trees.size()];

    for ( uint i = 0; i < destination.size(); i++ )
    {
        sendtag.push_back( 1 );

        recvtag.push_back( 1 );
    }

    //    while ( loopbool )
    {
        // start out by setting loopbool=0
        mybool = false;
        loopbool = false;
        counter2 = 0;

        /*
            for ( uint j = 0; j < destination.size(); j++ )
                {
                    message[j].clear();
                }
        */
        //    vector<bitset<M + N>> message[destination.size()];
        // save the originel list to overlap comp and com

        uint counter4 = 0;
        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {

            for ( auto i = ( *it ).Rbegin(); i != ( *it ).Rend(); i++ )
            {
                auto p = ( *it ).readRefineList( i );
                initialList[counter4].push_back( p.first );
            }
            counter4++;
        }

        // enforce 4: balance for each tree and extract the boundary elements from the list of elements to be refined

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            //            ( *it ).fourToOneP( start.at( counter2 ), end.at( counter2 ) );

            ( *it ).fourToOne();
            ( *it ).refinelistReset();
            counter2++;
        }

        // separates the list between local boundary and non-local boundary
        it3 = seeds.begin();

#if ( 1 )
        counter5 = 0;

        for ( auto it = trees.begin(); it != trees.end(); it++ )
        {
            seedkey = ( *it3 );
            proc.level( seedkey, &topologylevel );

#if ( DEBUG )
            myfile << " ============================== " << endl;
            myfile << " seedkey " << seedkey << "level " << topologylevel << endl;
// cout<<GREEN<<"Top lovel "<<topologylevel<<"seedkey "<<seedkey<<RESET<<endl;
#endif
            ( *it ).refinelistReset();

#if ( 1 )
            innerloop = true;
            while ( innerloop )
            {
                innerloop = false;

                for ( auto i = ( *it ).Rbegin(); i != ( *it ).Rend(); i++ )
                {
                    combinedkey = 0;

                    for ( uint j = 0; j < 3 * topologylevel; j++ )
                    {
                        combinedkey[NM - j - 1] = seedkey[M - j - 1];
                    }

                    auto p = ( *it ).readRefineList( i );
                    // if the tag = 1 then we have already investigated this guy, continue
                    if ( p.second == 1 )
                    {
                        continue;
                    }
                    //                    cout<<RED<<p.first<<endl;

                    key = p.first;

                    ( *it ).level( key, &mylevel );

#if ( DEBUG )
                    myfile << " elem key " << key << "level " << mylevel << endl;
#endif
                    combinedlevel = topologylevel + mylevel;

                    for ( uint j = 0; j < 3 * mylevel; j++ )
                    {
                        combinedkey[NM - 3 * ( topologylevel ) - j - 1] = key[N - j - 1];
                    }

                    ktcom = combinedkey;

#if ( 1 )
                    for ( uint direction = 0; direction < 3; direction++ )
                    {
                        if ( !( *it ).isBoundary( key, direction ) )
                        {
                            continue;
                        }
                        //  cout << " direction " << direction << endl;
                        //   cout << " combined key " << combinedkey << endl;
                        findFlipLevel( combinedkey, &combinedlevel, &changedirectionlevel, &direction );

#if ( DEBUG )
                        myfile << " assembled key " << combinedkey << endl;
#endif

                        if ( changedirectionlevel != 0 )
                        {
                            flipForNbr( combinedkey, &combinedlevel, &changedirectionlevel, &direction );

#if ( DEBUG )
                            myfile << " flipped key " << combinedkey << " change level " << changedirectionlevel << " direction "
                                   << direction << endl;
#endif
#if ( 1 )
                            //     cout << " combined key " << combinedkey << " combined level " << combinedlevel << endl;
                            //                     getNbrSeedLevel( combinedkey, topologylevel, &nbrseedlevel, proc ); //toplogylevel is
                            // wrong
                            getNbrSeedLevel( combinedkey, maxProcLevel, &nbrseedlevel, proc );

#if ( DEBUG )
                            myfile << " nbrseedlevel  " << nbrseedlevel << endl;
#endif

                            seednbrkey = 0;
                            for ( uint k = 0; k < 3 * nbrseedlevel; k++ )
                            {
                                seednbrkey[M - 1 - k] = combinedkey[NM - k - 1];
                                // cout<<RED<<combinedkey[NM-k-1]<<RESET<<endl;
                            }

//   cout << "seednbrkey " << BLUE << seednbrkey << RESET << endl;
#if ( DEBUG )
                            myfile << "  nbr seed key  " << seednbrkey << " nbrseedlevel  " << nbrseedlevel << endl;
#endif

                            nbrkey = 0;

                            for ( uint k = 0; k < N; k++ )
                            {
                                nbrkey[N - k - 1] = combinedkey[NM - 3 * ( nbrseedlevel ) - 1 - k];
                            }
//  cout << " nbrkey " << nbrkey << endl;
#if ( DEBUG )
                            myfile << "  nbr key  " << nbrkey << endl;
                            myfile << "Do I own It  " << isInSeed( seednbrkey, &counter ) << endl;
#endif

                            if ( isInSeed( seednbrkey, &counter ) )
                            {
                                auto it2 = std::next( trees.begin(), counter );
                                auto it4 = ( *it2 ).find( nbrkey );

                                // seedlevel
                                // ( proc ).level( seednbrkey, &nbrlevel );
                                // initial guess on nbr level as constructed
                                // this is a singularity for level zero, it might be negative for first level
                                // I remove this mannually
                                if ( combinedlevel >= nbrseedlevel )
                                {
                                    nbrlevel = combinedlevel - nbrseedlevel;
                                }
                                else
                                {
                                    nbrlevel = 1;
                                }
                                //      cout<<RED<<nbrlevel<<RESET<<endl;

                                // now find the real levels

                                if ( it4 != ( *it2 ).end() )
                                {
                                    ( *it2 ).level( it4->first, &nbrlevel );
                                    //         cout << YELLOW << nbrlevel << RESET << endl;
                                }
                                else
                                {
                                    // this condition imples level is lower no need for modification and search, I commented it out
                                    //                                nbrlevel                     = mylevel-1;
                                    // ( *it2 ).level( it4->first, &nbrlevel );
                                    nbrlevel = nbrlevel - 1;
                                    nbrkey[N - 3 * nbrlevel - 1] = 0;
                                    nbrkey[N - 3 * nbrlevel - 2] = 0;
                                    nbrkey[N - 3 * nbrlevel - 3] = 0;

                                    //            cout <<RED<< N - 3 * nbrlevel - 1 << " " << mylevel << RESET<<endl;
                                    //                                      auto it4=(*it2).find(nbrkey);
                                    //                                      (*it).level(it4->first,&nbrlevel);
                                }

                                // cout << "nbrkey " << RED << nbrkey << "  count " << counter << RESET << endl;
                                // cout<<"inside"<<endl;
                                if ( ( *it2 ).find( nbrkey ) == ( *it2 ).end() )
                                {
                                    cout << Com.myrank << " " << nbrkey << " direction " << direction << " level " << nbrlevel << endl;
                                    throw std::runtime_error( "error in finding key refinement" RESET );
                                }

                                nbrcomplevel = nbrlevel + nbrseedlevel;

                                //                             cout<<RED<<nbrlevel<<" "<<nbrseedlevel<<RESET<<endl;
                                //                             cout<<YELLOW<<nbrcomplevel<<" " <<combinedlevel<<RESET<<endl;
                                // if not in the list add to the list, dont forget to modify here
                                // watch out this element might have already been tagged

                                if ( nbrcomplevel < combinedlevel && ( *it2 ).isInRefineList( nbrkey ) == false )
                                {
                                    ( *it2 ).addToList( nbrkey );
                                    innerloop = true;
                                    // here is where boolwhile is affected
                                    loopbool = true;
                                    //       cout << RED "======================================" << endl;
                                    //       cout << "added to list" << nbrkey << endl;
                                    //       cout << "======================================" RESET << endl;
                                }
                                // cout<<seednbrkey<<RED<<seednbrkey<<RESET<<endl;
                            }

                            else
                            {
// sorts for communication

#if ( 1 )
// prepare for communication, pack all the elements with the destination, sender is   my_rank
#if ( DEBUG )

                                myfile << "*********************************  " << endl;
                                myfile << "inside else  nbr seed key  " << seednbrkey << " nbrseedlevel  " << nbrseedlevel << endl;
#endif
                                // find seednbrkey from global data
                                auto it5 = proc.find( seednbrkey );
                                //                        dest.push_back( it5->second[0] );

                                if ( it5 == proc.end() )
                                {
                                    throw std::runtime_error( "destination not found to send " );
                                }

                                auto it6 = destination.begin();
                                if ( it5 != proc.end() )
                                {
                                    it6 = find( destination.begin(), destination.end(), it5->second[0] );
                                }
                                else
                                {
                                    std::runtime_error( "Destination Not found" );
                                }
                                // combined key should be avoided
                                idx = it6 - destination.begin();
                                // cout << "---------------------->>> index" << RED << idx << RESET << endl;
                                // before pushing back, need to be careful about the key with all zeros,
                                // to accommodate this, the first bit (far right) is flipped and a sibling is sent
                                //  myfile<< " combined key before "<< combinedkey << "level "<<combinedlevel <<endl;
                                removeAllZeroSingularity( combinedkey, combinedlevel );
// cout << GREEN << " proc_id " << Com.myrank << "seednbrkey  " << seednbrkey << "combinedkey  " <<
// combinedkey
//     << "combined level " << combinedlevel << endl;
#if ( 0 )
                                message[idx].push_back( combinedkey );
// message[2].push_back( combinedkey );
#endif
#if ( DEBUG )
                                myfile << "idx  " << idx << " dest size " << destination.size() << endl;

                                myfile << "diff " << destination.size() - idx << endl;

                                myfile << " combined key pushed back  " << combinedkey << "level " << combinedlevel << endl;
                                myfile << "*********************************  " << endl;
#endif
//                       cout <<"message "<< message[i].at(0) << RESET << endl;
#endif
                            }
                        }
                        combinedkey = ktcom;
                    }
// here change the int value for ith element of the it-th tree
#endif
                    ( *it ).flipRefineElemTag( i );
                }
            }
#endif

#endif

            counter5++;
            it3 = std::next( it3, 1 );
        }

#endif

// do a neighborhood comm to tell the nbrs the size of the message to be recieved

#if ( 0 )
        for ( uint i = 0; i < destination.size(); i++ )
        {
            // prrfomr send/recieve only when the switch is on

            /*
                        if ( sendtag[i] == 0 )
                        {
                            continue;
                        }
            */
            // sending with boolean is supposed to be more portable

            //    sendbuff[i] = new char[message[i].size() * ( M + N )];

            sendbuff[i] = (char *)realloc( sendbuff[i], ( message[i].size() * ( M + N ) + 1 ) * sizeof( char ) );
            // counter6=0;
            //
            // if(sendbuff[i]==NULL && message[i].size()!=0)
            if ( sendbuff[i] == NULL )
            {
                cout << " myrank " << Com.myrank << " " << message[i].size() << endl;
                throw std::runtime_error( "not able to allocate sendbuf" );
            }
            for ( uint j = 0; j < message[i].size(); j++ )
            {
                tempkey = message[i].at( j );

#if ( DEBUG )
                myfile << " message to be sent " << tempkey << " destination " << destination.at( i ) << endl;
#endif
                for ( uint k = 0; k < M + N; k++ )
                {
                    if ( tempkey[k] == true )
                    {
                        sendbuff[i][j * ( M + N ) + k] = '1';
                    }
                    else
                    {
                        sendbuff[i][j * ( M + N ) + k] = '0';
                    }
                }
            }

            //            cout << "================================== " << endl;
            //            cout << RESET "messagesize " << message[i].size() << " destination " << destination.at( i ) << RESET << endl;
            //            cout << "================================== " << endl;
            // send messages to destinations and tag each meaasge with self rank (myrank)
            //
            MPI_Isend( sendbuff[i], message[i].size() * ( M + N ), MPI_CHAR, destination.at( i ), Com.myrank, MPI_COMM_WORLD,
                       &request1[i] );
        }

        for ( uint i = 0; i < destination.size(); i++ )
        {
            //             MPI_Iprobe(destination.at(i), destination.at(i), MPI_COMM_WORLD,&flag, &status);
            // this is a symmetric comm pattern therefore any
            //
            /*
                        if ( recvtag[i] == 0 )
                        {
                            continue;
                        }
            */
            MPI_Probe( destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &status );

            //   MPI_Probe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
            // while(flag!=0)
            {
                // MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,&flag, &status);
                //  cout << BLUE << flag << RESET << endl;
                // if ( flag == 1 )
                {
                    MPI_Get_count( &status, MPI_CHAR, &size );

                    // recvbuff[i] = new morton<M + N>[ size ];
                    // notice we revcieved the message by byte so need to specify,  to find the number of elements simly divide it by sizeof
                    // bool
                    //   recvbuff[i] = new char[size];
                    //   size+1 is to ensure deallocation
                    recvbuff[i] = (char *)realloc( recvbuff[i], ( size + 1 ) * sizeof( char ) );

                    // if(recvbuff[i]==NULL && size!=0)
                    if ( recvbuff[i] == NULL )
                    {
                        cout << " myrank " << Com.myrank << " " << size << endl;
                        throw std::runtime_error( "bad alloc in recv buffer" );
                    }

                    //           cout << "******************* " << size << endl;

                    MPI_Irecv( recvbuff[i], size, MPI_CHAR, destination.at( i ), destination.at( i ), MPI_COMM_WORLD, &request[i] );

                    MPI_Wait( &request[i], &status );
                    // if ( size != 0 )
                    {
                        //                cout << "rank " << Com.myrank << "recvbuff ??????????????? " << recvbuff[0][0] << endl;

                        // if that violates the balance condition add it to the list in the appropriate location

                        for ( uint j = 0; j < size / ( M + N ); j++ )
                        {
                            rcvkey = 0;

                            for ( uint k = 0; k < N + M; k++ )
                            {
                                if ( recvbuff[i][j * ( N + M ) + k] == '1' )
                                {
                                    rcvkey.flip( k );
                                }
                            }

//             cout << "p_id " << Com.myrank << " " << rcvkey << endl;
#if ( DEBUG )
                            myfile << " p_id  " << Com.myrank << " " << rcvkey << endl;
#endif
                            // restore if modified to handle singularity
                            // restoration needsto be done before this func

                            combinedLevel( rcvkey, &recvmessagecombinedlevel );

                            recoverAllZeroSingularity( rcvkey, recvmessagecombinedlevel );

                            findSeedLevelForRcvdMessage( rcvkey, &seedlevel, proc );

                            //                           cout << rcvkey << endl;
                            //                            cout << " seedlevel " << seedlevel << endl;

                            constructSeedKeyForRcvdMessage( rcvkey, seedlevel, seedkey );
                            // cout <<" myrank "<<Com.myrank <<" seedlevel " << seedlevel << endl;

                            /*      	    uint check;
                                               proc.level(seedkey,&check);
                                               cout<<RED<<check<<" "<<seedlevel<<RESET<<endl;
                                               if(seedlevel!=check)
                                              {
                                               throw std::runtime_error("levels calculated for seed inconsistent");
                                               }
                              */
                            constructElementKeyForRcvdMessage( rcvkey, seedlevel, elementkey );
                            //            cout << "element key " << elementkey << endl;

                            // cout << RED "rcvkey level " << recvmessagecombinedlevel << RESET << endl;

                            auto it7 = std::find( seeds.begin(), seeds.end(), seedkey );

                            if ( it7 == seeds.end() )
                            {
                                cout << " p_id  " << Com.myrank << " " << rcvkey << "  " << recvmessagecombinedlevel << endl;
                                cout << RED << "myrank " << Com.myrank << " " << seedkey << RESET << endl;
                                //  cout << "myrank " << Com.myrank << " " << seedkey << endl;

                                throw std::runtime_error( RED "seed not found in Balance Comm" RESET );
                            }

                            idex2 = std::distance( seeds.begin(), it7 );

                            auto it9 = std::next( trees.begin(), idex2 );

                            if ( recvmessagecombinedlevel >= seedlevel )
                            {
                                elementlevel = recvmessagecombinedlevel - seedlevel;
                            }
                            else
                            {
                                elementlevel = 1;
                            }
                            // cout << BLUE << elementlevel << RESET << endl;

                            // myfile <<" recv message combined level " <<  recvmessagecombinedlevel <<" seed level "<< seedlevel << endl;
                            // myfile <<" element level " << elementlevel<< " elementkey  "<< elementkey  << endl;

                            // recover the singularity
                            // error here, I flip the code such that we can get the level
                            // myfile << rcvkey << " " << recvmessagecombinedlevel << endl;
                            // myfile << elementkey << endl;
                            // myfile << seedkey << " " << seedlevel << endl;
                            // myfile << elementlevel << endl;
                            ( *it9 ).level( elementkey, &elementlevel );
                            // myfile << elementlevel << endl;

                            //                            recoverAllZeroSingularity( elementkey, elementlevel );

                            // get the real level
                            if ( ( *it9 ).isInMeshList( elementkey ) == false )
                            {
                                //   cout << "inside mesh" << elementkey << endl;
                                elementkey[N - 3 * ( elementlevel - 1 ) - 1] = 0;

                                elementkey[N - 3 * ( elementlevel - 1 ) - 2] = 0;

                                elementkey[N - 3 * ( elementlevel - 1 ) - 3] = 0;
                            }
                            //     myfile<<elementkey<<endl;
                            ( *it9 ).level( elementkey, &elementlevel );

                            localcombinedlevel = elementlevel + seedlevel;

                            // add if the level is lower

                            if ( localcombinedlevel < recvmessagecombinedlevel && ( *it9 ).isInRefineList( elementkey ) == false )
                            {
                                //                cout << BLUE << localcombinedlevel << "            " << recvmessagecombinedlevel << RESET
                                // << endl;

                                //               cout << "elementkey" << elementkey << endl;
                                ( *it9 ).addToList( elementkey );
                                loopbool = true;
                                /*                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;

                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                                                cout << "===================" << endl;
                */

                                // ????????????????????/
                                // modify the tags
                                mybool = true;
                            }

                            // see if this element exits, if yes, we are ok, if not add to tree referenced by pointer *it9
                        }
                    }
                }
            }
        }

#endif

#if ( 0 )
        /*
        if(mybool==false)
        {
        for(uint i=0;i<destination.size();i++)
        {
        destination.at(i)=0;
        }
        }
        */
        //
        //
        //       define recievers
        //
        //
        //

        if ( mybool == true )
        {
            intracount++;
        }

        if ( loopbool == 1 )
        {
            sb = 1;
        }
        else
        {
            sb = 0;
        }

        // this value is boolan,use MPI_BYTE, since bool size is implementation dependent, use sizeof(bool)

        // MPI_Iallreduce( &sb, &rb, sizeof( bool ), MPI_BYTE, MPI_MAX, MPI_COMM_WORLD, &request0 );

        // this value is boolan,use MPI_BYTE, since bool size is implementation dependent, use sizeof(bool)

        MPI_Iallreduce( &sb, &rb, 1, MPI_SHORT, MPI_MAX, MPI_COMM_WORLD, &request0 );

        // do some computation here
        //
        counter3 = 0;

        // switch the tag for initial list

        auto it = trees.begin();

        for ( uint i1 = 0; i1 < trees.size(); i1++ )
        {
            ( *it ).refineRefineList( initialList[i1] );
            it = std::next( it, 1 );
        }

        for ( uint i = 0; i < trees.size(); i++ )
        {
            initialList[i].clear();
        }

        // for debug for now do the refinem
        MPI_Wait( &request0, MPI_STATUS_IGNORE );

        for ( uint j = 0; j < destination.size(); j++ )
        {
            //   MPI_Wait( &request1[j], &status );

            while ( message[j].size() > 0 )
            {
                MPI_Wait( &request1[j], &status );
                message[j].clear();

                // delete[] sendbuff[j];
            }

            //     message[j].shrink_to_fit();
            //    delete[] recvbuff[j];
        }

#endif
        //
        // need to use the recv buffer
        //
        if ( rb == 1 )
        {
            loopbool = true;
            for ( auto it = trees.begin(); it != trees.end(); it++ )
            {
                ( *it ).refinelistReset();
            }
        }
        if ( loopbool == 1 )
        {
#if ( DEBUG )
            cout << RED << " loopbool " << loopbool << " rb " << rb << " sb " << sb << RESET << endl;
#endif
        }
    }

#if ( DEBUG )
    myfile.close();

#endif
    /*
      for(int i=0;i<destination.size();i++)
      {
       free(sendbuff[i]);
       free(recvbuff[i]);
      }
    */

    //   start.clear();
    //   end.clear();
    // assigne recieved elements to corresponding lists if needed

    // cout<<"============================"<<RED<<Com.myrank<<" "<<intracount<<endl;
}

#endif

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
uint TemplateForest<N, Nvalue, M, Mvalue, T>::findIndexInSeed( T &proc, morton<M> &seedkey )
{
    auto it5 = proc.find( seedkey );
    auto it6 = destination.begin();

    if ( it5 != proc.end() )
    {
        it6 = find( destination.begin(), destination.end(), it5->second[0] );
    }
    else
    {
        std::runtime_error( "Destination Not found" );
    }
    // combined key should be avoided
    uint idx = it6 - destination.begin();
    /*
    #if ( DEBUG )

        if ( idx >= destination.size() )
        {

            std::runtime_error( "Unacceptable Index for destination" );
        }
    #endif
    */
    return ( idx );
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest
<N, Nvalue, M, Mvalue, T>::appendToMessage( T &proc, morton<M> &seednbrkey, morton<M + N> &combinedkey, const uint combinedlevel )
{

    // there are two possibilities,
    // a) the neighborseed has level not greater than seedlevel which we only have one destination

    // morton<M> kn[4];
    uint idx;
    idx = findIndexInSeed( proc, seednbrkey );

    // if an element has a higher order neighbor, nonlocal neighbor algorithm will find the sibling which might belong to another processor,
    // this ponly happens if the neighbor level is higher so we can simply ignore this
    if ( idx < destination.size() )
    {
        removeAllZeroSingularity( combinedkey, combinedlevel );
        // cout << GREEN << " proc_id " << Com.myrank << "seednbrkey  " << seednbrkey << "combinedkey  " <<
        // combinedkey
        //     << "combined level " << combinedlevel << endl;

        message[idx].push_back( combinedkey );
        // message[2].push_back( combinedkey );
    }

#if ( 0 )
#endif
    /*
        myfile << "idx  " << idx << " dest size " << destination.size() << endl;

        myfile << "diff " << destination.size() - idx << endl;

        myfile << " combined key pushed back  " << combinedkey << "level " << combinedlevel << endl;
        myfile << "*********************************  " << endl;
    */
    //                       cout <<"message "<< message[i].at(0) << RESET << endl;
    // b) the neighbor seed is bigger than the seed, therefore we need to construct 4 seeds, send them to different destinations

    // no need to check this as this automatically implies that the neighbor has a higher level
}


// some instantiations
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N,Nvalue,M,Mvalue,T>::MPIStartUp()
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

    if ( MPI_ERROR_DISABLE == 1 )
    {
        if ( MPI_SUCCESS != MPI_Comm_set_errhandler( Com.mpicom, MPI_ERRORS_RETURN ) )
        {
            cout << BLUE << " Rank(" << Com.myrank << ") > Exit Code : " << MPI_ERROR_HANDLE_FAIL << RESET << endl;
            cout << BLUE << ReblAmrGetErrorEnum( MPI_ERROR_HANDLE_FAIL ) << RESET << endl;
            exit( 1 );
        }
    }
  
  duplicated=1;
  
}


const char *ReblAmrGetErrorEnum( ReblAmrResult error )
{
    switch ( error )
    {
        case SUCCESS:
            return "SUCCESS";
        case NUM_INPUT_ARGS:
            return (
            "\n \t INPUT_FAILURE: \n \n  3 inputs are needed: 1. geometry file (.stl), 2. topology level (int) , 3. refinement level (int) \n " );
        case PROC_LEVEL:
            return ( "\n PROCSIZE is not big enough to accomodate desired processor topology level\n" );
         case MPI_INIT_CHECK_FAIL:
            return ( "\n Error in checking to see if MPI is already initialized \n " );
        case MPI_INIT_FAIL:
            return ( "\n Error in MPI initialization\n " );
        case MPI_DUP_FAIL:
            return ( "\n Error in MPI Comm Duplication\n " );
        case COMSIZE_FAIL:
            return ( "\n Comsize sould be i*i, where 'i' is an integer\n " );
        case MPI_GET_RANK_FAIL:
            return ( "\n MPI_Comm_rank in PencilDcmp::MPIStartUp failed  \n " );
        case MPI_COMSIZE_FAIL:
            return ( "\n   MPI_Comm_size in PencilDcmp::MPIStartUp failed   \n " );
        case COMBINED_SIZE:
            return ( "\n   (M+N) can not be a multiply of 3 since  one bit is needed to recover singularity in parallel   \n " );
        case NO_SEED:
            return ( "\n   seed not found in the proc \n " );
        case MPI_ERROR_HANDLE_FAIL:
            return ( "\n   failed to in setting error handler for MPI \n " );
        case THOMAS_FAIL:
            return ( "\n   Thomas algorithm Failed \n " );
        case MESH_LEVEL:
            return ( "\n TREESIZE is not big enough to accomodate desired mesh level\n" );
   
    }

    return "<unknown>";
}


template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N,Nvalue,M,Mvalue,T>::currentDateTime()
{
    time_t    now = time( 0 );
    struct tm tstruct;
    // char       nameAppendix[80];
    tstruct = *localtime( &now );
    strftime( nameAppendix, sizeof( nameAppendix ), "%Y_%m_%d_%X", &tstruct );
    // return nameAppendix;
}


template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N,Nvalue,M,Mvalue,T>::runInfo()
{
    currentDateTime();
 //   double meshSize = nChunk * nChunk * nChunk * nxChunk * nyChunk * nzChunk;

    if ( Com.myrank == 0 )
    {
        ofstream    ReblAmrOut;
        std::string filename = "ReblAmr_";
        filename.append( nameAppendix );
        filename.append( "_out" );

        ReblAmrOut.open( filename );
        ReblAmrOut << "*********************************************************\n" << endl;
        ReblAmrOut << "              ReblAmr: " << nameAppendix << "        \n" << endl;
        ReblAmrOut << "*********************************************************\n" << endl;
        ReblAmrOut << "---------------------------------------------------------\n" << endl;
        ReblAmrOut << "                     Parameters                        \n" << endl;
        ReblAmrOut << "---------------------------------------------------------\n" << endl;
        ReblAmrOut << "# STL geometry mesh quality metrics report = " << CHECK_MESH  <<endl;
        ReblAmrOut << "# PROC container size = " << PROCSIZE  <<endl;
        ReblAmrOut << "# TREE container size = " << TREESIZE  <<endl;
        ReblAmrOut << "# nxChunk = " << npx  <<endl;
        ReblAmrOut << "# nyChunk = " << npy  <<endl;
        ReblAmrOut << "# nzChunk = " << npz  <<endl;
        ReblAmrOut << "# Processor Tree Topology Level : " <<   endl;
        ReblAmrOut << "# Mesh Level (Tree Level)  " <<   endl;
 
        if(ZOLTAN_ON==1)
        {
        ReblAmrOut << "# Zoltan = ON" << endl;
        if(WEIGHT==1) 
        {
        ReblAmrOut << "# Perfoming Weighted Partitioning" << endl;
        }
        } 

        if ( WR == 0 )
        {
            ReblAmrOut << "# I/O = OFF"
                    << "\n"
                    << endl;
        }
        else if(WR==1)
        {
            ReblAmrOut << "# I/O =  " << WR<<" (i.e. writing only the points) " << "\n" << endl; 
        }
       else
       {
            ReblAmrOut << "# I/O writing out AMR blocks  " <<WR<<" (i.e. Writing the AMR Blocks)" << "\n" << endl; 
       } 
       ReblAmrOut << "---------------------------------------------------------\n" << endl;
       ReblAmrOut << "\n         # Global mesh size = "<<  FormatWithCommas(meshSize)<<" ~ "<<setprecision(2)<<std::fixed<< meshSize/(1.e6)<<" (M) " << "\n"<< endl;
       ReblAmrOut << "---------------------------------------------------------\n" << endl;
       ReblAmrOut.close();
      }
}

template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N,Nvalue,M,Mvalue,T>::updateSeeds(T &proc)
{

 seeds.clear();
    
// since seeds are used in here 

for ( auto it = proc.begin(); it != proc.end(); it++ )
 {
          //cout<<" my_rank "<<Com.myrank <<"  "<<it->first<<" "<<it->second[0]<<endl;
          cout<<" my_rank "<<Com.myrank <<"  "<<it->first<<endl;
  //      if ( it->second[0] == Com.myrank )
        {
           seeds.push_back( it->first );

 //            cout <<BLUE<< it->first << " XC " << XC[0] << " " << XC[1] << " " << XC[2] <<RESET<< endl;
        }
  }

    
 cout<<" seedsize  " <<seeds.size()<<" porocsize  "<<proc.size()<<endl;


morton<N> key;
/*
for(int i=0;i<proc.size();i++)
{
key=proc.mesh.at(i);

}
*/

}


template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N,Nvalue,M,Mvalue,T>::updateSeedsAndTrees(T &proc)
{

  bitlist<M> seedsTmp; 
  bitlist<M> seedsDiff; 
/*
// save the old seeds 
  for(auto it=seeds.begin();it!=seeds.end();it++)
  {
    seedsTmp.push_back((*it));

    cout<<" temp  "<<(*it)<<endl;
  }

    cout<<" temp  "<<seedsTmp.size()<<endl;
*/
 
 seeds.clear();
 //trees.clear();     
 real X[6];
 real X0[6];
 real len[3];
 real XC[3];
 uint count=0;
 bool bol=true;

 vector<Tree<N,Nvalue>> tmp;
int cnt=0;
         cout<<" original tree size  "<<trees.size()<<endl;
  for(auto it=trees.begin();it!=trees.end();it++)
  {
    (*it).getCoords(X0);
//    cout<<" ???????????? "<<X0[0]<<" "<<X0[1]<<endl;
    bol=true;
     
    for ( auto it1 = proc.begin(); it1 != proc.end(); it1++ )
    {
            //cout<<" my_rank "<<com.myrank <<"  "<<it->first<<" "<<it->second[0]<<endl;
        if ( it1->second[0] == Com.myrank )
        {
            proc.enclosingBox( it1->first, X );
    
          // cout<<" !!!!! "<<X[0]<<" "<<X[1]<<endl;

            if((fabs(X[0]-X0[0])<1.e-6) && (fabs(X[1]-X0[1])<1.e-6) && (fabs(X[2]-X0[2])<1.e-6) && (fabs(X[3]-X0[3])<1.e-6) && (fabs(X[4]-X0[4])<1.e-6) &&  (fabs(X[5]-X0[5]))<1.e-6 )
            {
//               cout<<" inside"<<endl;
               bol=false;
               cnt++;
               break;              
            }
         }
    } 
    if(bol==true)
    {
     tmp.push_back(*it);
   
    } 

  }
 
  while(tmp.size()!=0)
  {
   trees.erase(std::find(trees.begin(),trees.end(),tmp.at(tmp.size()-1)));
   tmp.pop_back();
   }

//           cout<<" temp list size "<<tmp.size()<<endl;
//           cout<<" cnt "<<cnt<<endl;

    for ( auto it = proc.begin(); it != proc.end(); it++ )
    {
//            cout<<" my_rank "<<com.myrank <<"  "<<it->first<<" "<<it->second[0]<<endl;
        if ( it->second[0] == Com.myrank )
        {
            proc.enclosingBox( it->first, X );

            XC[0] = ( X[0] + X[1] ) * 0.5;
            XC[1] = ( X[2] + X[3] ) * 0.5;
            XC[2] = ( X[4] + X[5] ) * 0.5;

            len[0] = fabs( X[1] - X[0] );
            len[1] = fabs( X[3] - X[2] );
            len[2] = fabs( X[5] - X[4] );

            // regrdless of the topology, we need a tree in tree list

            seeds.push_back( it->first );

//            cout <<BLUE<< it->first << " XC " << XC[0] << " " << XC[1] << " " << XC[2] <<RESET<< endl;
            count++;
        }
    }

    for ( auto it = seeds.begin(); it != seeds.end(); it++ )
    {
       proc.enclosingBox( (*it), X );

       bol=true;

     for ( auto it1 = trees.begin(); it1 != trees.end(); it1++ )
    {
            //cout<<" my_rank "<<com.myrank <<"  "<<it->first<<" "<<it->second[0]<<endl;
             (*it1).getCoords(X0);
    
          // cout<<" !!!!! "<<X[0]<<" "<<X[1]<<endl;

            if((fabs(X[0]-X0[0])<1.e-6) && (fabs(X[1]-X0[1])<1.e-6) && (fabs(X[2]-X0[2])<1.e-6) && (fabs(X[3]-X0[3])<1.e-6) && (fabs(X[4]-X0[4])<1.e-6) &&  (fabs(X[5]-X0[5]))<1.e-6 )
            {
//               cout<<" inside"<<endl;
               bol=false;
//               cnt++;
               break;              
            }
      }

      if(bol==true)
       {
            XC[0] = ( X[0] + X[1] ) * 0.5;
            XC[1] = ( X[2] + X[3] ) * 0.5;
            XC[2] = ( X[4] + X[5] ) * 0.5;

            len[0] = fabs( X[1] - X[0] );
            len[1] = fabs( X[3] - X[2] );
            len[2] = fabs( X[5] - X[4] );
          trees.push_back( Tree<N, Nvalue>( len, XC ) );
        }

    } 
 
            

         //   trees.push_back( Tree<N, Nvalue>( len, XC ) );

//            seeds.push_back( it->first );


 //           cout <<BLUE<< it->first << " XC " << XC[0] << " " << XC[1] << " " << XC[2] <<RESET<< endl;
            count++;
    
   
 




   cout<<" count = "<<count<<endl;
   cout<<" seeds size = "<<seeds.size()<<endl;
   cout<<" trees size = "<<trees.size()<<endl;

}

/*
template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T>
void TemplateForest<N,Nvalue,M,Mvalue,T>::selectkeysToIgnore(T &proc, double xyz[6])
{

double center[6];
bool bol0,bol1,bol2;

for(auto it=proc.begin();proc.end();it++ )
{
   key=it->first;
  
   Proc.centroid( key, center);

   bol0= (center[0] < xyz[1] && center[0]>xyz[0]); 
   bol1= (center[1] < xyz[3] && center[1]>xyz[2]); 
   bol2= (center[2] < xyz[5] && center[2]>xyz[4]); 

  if(bol0 && bol1 && bol2 )   
  {
     
  }
}

}
*/


template class TemplateForest< TREESIZE, real, PROCSIZE, uint, Tree<PROCSIZE, uint>>;
//template class TemplateForest< TREESIZE, real, PROCSIZE, uint, FullTree<PROCSIZE, uint>>;
template class TemplateForest<  TREESIZE, real, WSIZE, uint, FullTree<WSIZE, uint>>;
