#ifndef _DATATYPE_H_
#define _DATATYPE_H_
#include "definitions.h"

// this is only to be able to redundantly use the same nonemclature in other namespaces
// does not mean mych really, remove it if not like it

namespace Abstraction
{
    typedef enum
    {   
	type_byte,
        type_char,
        type_unsigned_char,
        type_short,
        type_unsigned_short,
        type_int,
        type_unsigned_int,
        type_long,
        type_unsigned_long,
        type_float,
        type_double
    } DataType;
}

/*!
 * \class DataType 
 * \brief this class abstracts the MPI_Datatypes
 * \details using template initializations one can template MPI functions
 * using template specialization for further details see, 
 *  https://chuckaknight.wordpress.com/2013/03/13/intrinsic-type-conversion-using-template-specialization/
 *  since BYTE is not a native type is C++, for MPI_TYPE, the nullptr_t (null pointer type ) is utilized
*/

template <class T>
Abstraction::DataType getAbstractionDataType()
{
    throw std::runtime_error("Intrinsic type not supported by the abstraction");
}

/*!< Specilizations for the template class */
// the nullptr type has to be used here, which is nullptr_t
// nullptr_t will set the MPI_DATATYPE to MPI_BYTE
// unofrtunately one still has to specify the template parameter as nullptr_t
template <>
inline Abstraction::DataType getAbstractionDataType<nullptr_t>()
{
    return Abstraction::type_byte;
}

template <>
inline Abstraction::DataType getAbstractionDataType<char>()
{
    return Abstraction::type_char;
}

template <>
inline Abstraction::DataType getAbstractionDataType<unsigned char>()
{
    return Abstraction::type_unsigned_char;
}

template <>
inline Abstraction::DataType getAbstractionDataType<short>()
{
    return Abstraction::type_short;
}

template <>
inline Abstraction::DataType getAbstractionDataType<unsigned short>()
{
    return Abstraction::type_unsigned_short;
}

template <>
inline Abstraction::DataType getAbstractionDataType<int>()
{
    return Abstraction::type_int;
}

template <>
inline Abstraction::DataType getAbstractionDataType<unsigned int>()
{
    return Abstraction::type_unsigned_int;
}

template <>
inline Abstraction::DataType getAbstractionDataType<long>()
{
    return Abstraction::type_long;
}

template <>
inline Abstraction::DataType getAbstractionDataType<unsigned long>()
{
    return Abstraction::type_unsigned_long;
}

template <>
inline Abstraction::DataType getAbstractionDataType<float>()
{
    return Abstraction::type_float;
}

template <>
inline Abstraction::DataType getAbstractionDataType<double>()
{
    return Abstraction::type_double;
}


























#endif
