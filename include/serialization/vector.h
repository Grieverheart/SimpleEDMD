#ifndef __SERIALIZATION_VECTOR_H
#define __SERIALIZATION_VECTOR_H

#include "archive.h"
#include <type_traits>
#include <vector>

template<typename condition, typename R = void >
using EnableIf = typename std::enable_if<condition::value, R>::type;
template<typename condition, typename R = void >
using EnableIfNot = typename std::enable_if<!condition::value, R>::type;

//Seserialization
template<typename T>
EnableIf<std::is_trivially_copyable<T>,
void> serialize(Archive& ar, const std::vector<T>& vec){
    auto size = vec.size();
    serialize(ar, size);
    serialize(ar, vec.data(), size);
}

template<typename T>
EnableIfNot<std::is_trivially_copyable<T>,
void> serialize(Archive& ar, const std::vector<T>& vec){
    auto size = vec.size();
    serialize(ar, size);
    for(const auto& el: vec) serialize(ar, el);
}

//Deserialization
template<typename T>
EnableIf<std::is_trivially_copyable<T>,
void> deserialize(Archive& ar, std::vector<T>* vec){
    using size_type = typename std::vector<T>::size_type;
    size_type size = 0;
    deserialize(ar, &size);
    vec->resize(size);
    ar.read(vec->data(), size * sizeof(T));
}

template<typename T>
EnableIfNot<std::is_trivially_copyable<T>,
void> deserialize(Archive& ar, std::vector<T>* vec){
    using size_type = typename std::vector<T>::size_type;
    size_type size = 0;
    deserialize(ar, &size);
    for(size_type n = 0; n < size; ++n){
        auto el = T();
        deserialize(ar, &el);
        vec->push_back(std::move(el));
    }
}

#endif
