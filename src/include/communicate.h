#ifndef _COMMUNICATE_H_
#define _COMMUNICATE_H_
#include "definitions.h"

/*!
 * \struct MpiCom
 * \brief class for embedding data related to the communicator
 */

struct MpiCom
{
MPI_Comm mpicom; /*!< Communicator  */
integer myrank; /*!< rank of the processor*/
integer comsize; /*!< size of the communicator*/
MpiCom() {};     /*!< constructor */
MpiCom(MPI_Comm Com); /*!< second constructor */
};

/*!
 * \struct Message
 * \brief struct that embeds data related to the message and envelope
 * */

struct Message  
{
uint count; /*!< count of the message */
int tag; /*!< tag of the message */
uint sender; /*!< who is sending */
uint reciever; /*!< who is the destination*/
void *buf; /*!<pointer to the buffer of the data */
MPI_Datatype datatype; /*<! Templated using template initialization */
MPI_Status status;  
MPI_Request request; 

void print() /*!< printf the message info and is for debugging*/
{
cout<<RED "count "<<count<<endl;
cout<<"tag "<<tag<<endl;
cout<<"sender "<<sender<<endl;
cout<<"reciever " <<reciever<<endl;
cout<<"Datatype " RESET<<datatype<<endl;

};

};

/*!
 *    \class CommPoint2Point 
 *    \brief   A template wrapper around MPI
 functions for point to point communication 
 *    \details class CommPoint2Point abstracts point to point communications 
 *   and functions
 *  
 */


template <class Type>
class CommPoint2Point
{
private:
MpiCom Com;
Message Msg;

public:
CommPoint2Point(void *buff,uint size, uint *tg, uint snd, uint rcv, uint type,MPI_Comm com); /*!< full constructor*/ 
CommPoint2Point(void *buff,uint size, uint snd, uint rcv); /*!< default constructor */
CommPoint2Point(void *buff,uint size); /*!<deferring the assigment of sender and reciever */
void assignSender(uint sndr); /*! assigns the sender of the message*/
void assignReciever(uint rcv);  /*! assigns the destination the message*/ 
integer myRank(); /*!< gets the ranl of the processor*/
integer mySize();
void Irecv();     /*!<non-blocking recieve*/
void Isend();     /*!< non-blocking send*/
void recv();      /*!<blocking recieve*/
void send();      /*!<blocking send*/
//MPI_Request chechRequest();
void getOffset(uint myvalue,uint *offset); /*!< calculates the offset of each processor for variable myvalue*/
~CommPoint2Point() {}; /*!<destructor*/
};

/*!
 *   \class CommCollective 
 *   \brief is a template  wrapper around MPI functions for collective communicatios
 *   \details  This class is specialized for collective  communications 
 *   and functions
 *  
 */

template<class Type>
class CommCollective 
{
    private:
    MpiCom Com;
    Message Msg;

    public:
    CommCollective( void *buf, uint size, uint root); /*!< constructor */
    CommCollective(void *buff,uint size);     /*!< constructor specialized for root=comsize-1*/
    void bcast(); /*!<blockin broadcast   */
    void Ibcast(); /*!<non-blocking broadcast*/
    void getTotalNumber(uint *offset,uint *myvalue,uint* totalvalue); /*!<calculates the world population for a given variable  */ 
    void IgetTotalNumber(uint *offset,uint *myvalue,uint* totalvalue); /*!<calculates the world population for a given variable in a non-blocking fashion */ 
    void waitOnRequest();
    ~CommCollective() {}; /*!<destructor  */
};

#endif





