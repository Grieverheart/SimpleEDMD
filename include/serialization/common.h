#ifndef __SERIALIZATION_COMMON_H
#define __SERIALIZATION_COMMON_H

#include "archive.h"
#include <type_traits>

template<typename condition, typename R = void >
using EnableIf = typename std::enable_if<condition::value, R>::type;
template<typename condition, typename R = void >
using EnableIfNot = typename std::enable_if<!condition::value, R>::type;

//Seserialization
template<typename T>
EnableIf<std::has_trivial_copy_constructor<T>,
void> serialize(Archive& ar, const T& val){
    ar.write(&val, sizeof(T));
}

//NOTE: In principle we should also write 'num'.
template<typename T>
EnableIf<std::has_trivial_copy_constructor<T>,
void> serialize(Archive& ar, const T* val, size_t num){
    ar.write(val, num * sizeof(T));
}

//Deserialization
template<typename T>
EnableIf<std::has_trivial_copy_constructor<T>,
void> deserialize(Archive& ar, T* val){
    ar.read(val, sizeof(T));
}

//NOTE: In principle we should also write 'num'.
template<typename T>
EnableIf<std::has_trivial_copy_constructor<T>,
void> deserialize(Archive& ar, T* val, size_t num){
    ar.read(val, num * sizeof(T));
}

#endif
