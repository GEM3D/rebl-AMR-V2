#include "communicate.h"
#include "datatype.h"

/** communicate class member functions    CommPoint2Point
 *    This Class is a wrapper around MPI functions used in this project
 *
 *  Templating the "Intrinsic Type Conversion Using Template Specialization"
 *  message tag and mpi_communicators are assigned by default if user doesnt assign them
 *  default tag is 0 and default Communicator is MPI_COMM_WORLD
 */

static MPI_Datatype ConvertType( Abstraction::DataType type );

template <class Type>
CommPoint2Point<Type>::CommPoint2Point( void *buff, uint size, uint *tg, uint snd, uint rcv, uint type,
                                        MPI_Comm Comm ) /*!< If communicator is not specified it will be set as MPI_COMM_WORLD by default*/
{
    if ( Comm != MPI_COMM_NULL )
    {
        Com.mpicom = Comm;
    }
    else
    {
        Com.mpicom = MPI_COMM_WORLD;
    }

    MPI_Comm_rank( Com.mpicom, &Com.myrank );
    MPI_Comm_size( Com.mpicom, &Com.comsize );

    if ( size != 0 )
    {
        Msg.count = size;
    }
    else
    {
        MPI_Finalize();
        throw std::runtime_error( "Message Size can not be zero" );
    }

    if ( tg == nullptr )
    {
        Msg.tag
        = 0; /*!< default tag for message send and recieve is "0", MPI_ANY_TAG is onlu OK for recieve commnad so it can not be used here*/
    }
    else
    {
        Msg.tag = *tg;
    }

    Msg.sender = snd;
    Msg.reciever = rcv;

    Msg.datatype = ConvertType( getAbstractionDataType<Type>() );

    // cout<<"Sender "<<Msg.sender<<endl;

    // cout<<"Reciever "<<Msg.reciever<<endl;

    // cout<<"type "<<type<<endl;

    Msg.buf = buff;
}
template <class Type>
int CommPoint2Point<Type>::myRank()
{
    return ( Com.myrank );
}

template <class Type>
int CommPoint2Point<Type>::mySize()
{
    return ( Com.comsize );
}

template <class Type>
CommPoint2Point<Type>::CommPoint2Point( void *buff, uint size, uint snd, uint rcv ) /*!< another constructor for Communicator class*/
{

    Com.mpicom = MPI_COMM_WORLD;

    MPI_Comm_rank( Com.mpicom, &Com.myrank );
    MPI_Comm_size( Com.mpicom, &Com.comsize );

    if ( size != 0 )
    {
        Msg.count = size;
    }
    else
    {
        MPI_Finalize();
        throw std::runtime_error( "Message Size can not be zero" );
    }

    Msg.tag = 0;
    Msg.sender = snd;
    Msg.reciever = rcv;
    Msg.datatype = ConvertType( getAbstractionDataType<Type>() );
    Msg.buf = buff;
}

template <class Type>
CommPoint2Point
<Type>::CommPoint2Point( void *buff, uint size ) /*!< another constructor for Communicator class, deferes assigning sender and reciever*/
{
    Com.mpicom = MPI_COMM_WORLD;

    MPI_Comm_rank( Com.mpicom, &Com.myrank );
    MPI_Comm_size( Com.mpicom, &Com.comsize );

    if ( size != 0 )
    {
        Msg.count = size;
    }
    else
    {
        MPI_Finalize();
        throw std::runtime_error( "Message Size can not be zero" );
    }
    Msg.tag = 0;
    Msg.datatype = ConvertType( getAbstractionDataType<Type>() );
    Msg.buf = buff;
}

template <class Type>
void CommPoint2Point<Type>::assignSender( uint sndr )
{
    Msg.sender = sndr;
}

template <class Type>
void CommPoint2Point<Type>::assignReciever( uint rcvr )
{
    Msg.reciever = rcvr;
}

//===============================================================================
//
//                      Some General Communication functions
//
//===============================================================================
template <class Type>
void CommPoint2Point<Type>::Irecv()
{
    Msg.datatype = ConvertType( getAbstractionDataType<Type>() );
    MPI_Irecv( Msg.buf, Msg.count, Msg.datatype, Msg.sender, Msg.tag, Com.mpicom, &( Msg.request ) );
}

template <class Type>
void CommPoint2Point<Type>::Isend()
{
    MPI_Isend( Msg.buf, Msg.count, ConvertType( getAbstractionDataType<Type>() ), Msg.reciever, Msg.tag, Com.mpicom, &( Msg.request ) );
}

template <class Type>
void CommPoint2Point<Type>::recv()
{
    MPI_Recv( Msg.buf, Msg.count, ConvertType( getAbstractionDataType<Type>() ), Msg.sender, Msg.tag, Com.mpicom, &( Msg.status ) );
}

template <class Type>
void CommPoint2Point<Type>::send()
{
    //cout << "this tag " << Msg.tag << endl;
    MPI_Send( Msg.buf, Msg.count, ConvertType( getAbstractionDataType<Type>() ), Msg.reciever, Msg.tag, Com.mpicom );
}
static MPI_Datatype ConvertType( Abstraction::DataType type )
{
    switch ( type )
    {
        case Abstraction::type_byte:
            return MPI_BYTE;
        case Abstraction::type_char:
            return MPI_CHAR;
        case Abstraction::type_unsigned_char:
            return MPI_UNSIGNED_CHAR;
        case Abstraction::type_short:
            return MPI_SHORT;
        case Abstraction::type_unsigned_short:
            return MPI_UNSIGNED_SHORT;
        case Abstraction::type_int:
            return MPI_INT;
        case Abstraction::type_unsigned_int:
            return MPI_UNSIGNED;
        case Abstraction::type_long:
            return MPI_LONG;
        case Abstraction::type_unsigned_long:
            return MPI_UNSIGNED_LONG;
        case Abstraction::type_float:
            return MPI_FLOAT;
        case Abstraction::type_double: // cout<<"MPI_DOUBLE "<<MPI_DOUBLE<<endl;
            return MPI_DOUBLE;
    };
    throw std::runtime_error( "MPI_Datatype Convert( Abstraction::DataType ) failed" );
}

// I wanted to template this fucntion too but due to the presence of nullptr_t is my MPI_Templitization, this did fail
//  since + sign for nullptr's is undefined, we can come back to this after we decide on what to do on nullptr_t
//  as I reserved that for MPI_TYPE_BYTE which I need for Morton encoding

template <class Type>
void CommPoint2Point<Type>::getOffset( uint myvalue, uint *offset ) /*!<gets the number of elements before the current processor, it is
                                                                         needed for globally defining the elements */
{
    uint recvbuf = 0;
    // CommPoint2Point OffSet(void *buff,uint size,uint type);

    if ( Com.comsize == 1 )
    {
        *offset = 0;
    }
    else
    {

        if ( 1 )
        {
            if ( Com.myrank > 0 )
            {
                Msg.sender = Com.myrank - 1;
                Msg.buf = &recvbuf;
                //  Msg.print();
                recv();
                ( *offset ) = recvbuf;
            }

            if ( Com.myrank < ( Com.comsize - 1 ) )
            {
                Msg.reciever = Com.myrank + 1;
                recvbuf = recvbuf + myvalue;
                Msg.buf = &recvbuf;

                send();
            }
        }
        else
        {
            // kept open to add RMA features
            // printf( "RMA is a boolean set it to 1 for true and 0 for false\n" );
        }
    }
}

template <class Type>
CommCollective<Type>::CommCollective( void *buff, uint size, uint root ) /*!< another constructor for Communicator class*/
{

    Com.mpicom = MPI_COMM_WORLD;

    MPI_Comm_rank( Com.mpicom, &Com.myrank );
    MPI_Comm_size( Com.mpicom, &Com.comsize );

    if ( size != 0 )
    {
        Msg.count = size;
    }
    else
    {
        MPI_Finalize();
        throw std::runtime_error( "Message Size can not be zero" );
    }

    Msg.tag = 0;
    Msg.sender = root;
    Msg.reciever = root;
    Msg.datatype = ConvertType( getAbstractionDataType<Type>() );
    Msg.buf = buff;
}

template <class Type>
CommCollective
<Type>::CommCollective( void *buff,
                        uint size ) /*!< another constructor for Communicator class, this one sets the root as comsize-1 for hdf5*/
{

    Com.mpicom = MPI_COMM_WORLD;

    MPI_Comm_rank( Com.mpicom, &Com.myrank );
    MPI_Comm_size( Com.mpicom, &Com.comsize );

    if ( size != 0 )
    {
        Msg.count = size;
    }
    else
    {
        MPI_Finalize();
        throw std::runtime_error( "Message Size can not be zero" );
    }

    Msg.tag = 0;
    Msg.sender = Com.comsize - 1;
    Msg.reciever = Com.comsize - 1;
    Msg.datatype = ConvertType( getAbstractionDataType<Type>() );
    Msg.buf = buff;
}

template <class Type>
void CommCollective<Type>::bcast()
{
    MPI_Bcast( Msg.buf, Msg.count, ConvertType( getAbstractionDataType<Type>() ), Msg.sender, Com.mpicom );
}

template <class Type>
void CommCollective<Type>::Ibcast()
{
    MPI_Ibcast( Msg.buf, Msg.count, ConvertType( getAbstractionDataType<Type>() ), Msg.sender, Com.mpicom, &( Msg.request ) );
    // dont forget to check the request by MPI_Wait
    // MPI_Wait(&(Msg.request),MPI_STATUS_IGNORE);
}

/*!< I spoecialized this function since this is only needed for unsigned ints */
template <>
void CommCollective<uint>::getTotalNumber( uint *offset, uint *myvalue, uint *totalvalue )
// void CommCollective<uint>::getTotalNumber(uint *totalvalue) /*!< method 0 communicates with MPI-I and method 1 uses one-sided
// communication */
{
    uint c = ( *offset ) + ( *myvalue );
    Msg.buf = &c;

    /*!< the largest offset belongs to the processor with highest rank, add this to its number of cubes
       will give us the total value */
    if ( Com.comsize > 1 )
    {
        //      dispatcher = Com.comsize - 1;
        bcast();
    }
    *totalvalue = c;
}

template <class Type>
void CommCollective<Type>::waitOnRequest()
{
    if ( Com.comsize > 1 )
    {
        MPI_Wait( &( Msg.request ), MPI_STATUS_IGNORE );
    }
}

template <>
void CommCollective<uint>::IgetTotalNumber( uint *offset, uint *myvalue, uint *totalvalue )
{
    *totalvalue = ( *offset ) + ( *myvalue );
    Msg.buf = totalvalue;

    /*!< the largest offset belongs to the processor with highest rank, add this to its number of cubes
       will give us the total value */
    if ( Com.comsize > 1 )
    {
        //      dispatcher = Com.comsize - 1;
        Ibcast();
    }
    //*totalvalue=c;
}

template class CommPoint2Point<int>;
template class CommPoint2Point<uint>;
template class CommPoint2Point<nullptr_t>;
template class CommPoint2Point<double>;
template class CommPoint2Point<float>;

template class CommCollective<int>;
template class CommCollective<uint>;
template class CommCollective<nullptr_t>;
template class CommCollective<double>;
template class CommCollective<float>;
