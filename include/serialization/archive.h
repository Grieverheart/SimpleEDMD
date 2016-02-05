#ifndef __ARCHIVE_H
#define __ARCHIVE_H

#include <cstddef>
#include <type_traits>
#include <vector>

class Archive{
public:
    Archive(void);
    Archive(const void* data, size_t size);
    ~Archive(void);
    void write(const void* data, size_t size);
    void read(void* data, size_t size);
    const char* data(void)const;
    size_t size(void)const;
private:
    char* data_;
    size_t position_;
    size_t size_;
};


//TODO: Perhaps move to separate headers.

template<typename condition, typename R = void >
using EnableIf = typename std::enable_if<condition::value, R>::type;

template<typename condition, typename R = void >
using EnableIfNot = typename std::enable_if<!condition::value, R>::type;

template<typename T>
EnableIf<std::is_trivially_copyable<T>,
void> serialize(Archive& ar, const T& val){
    ar.write(&val, sizeof(T));
}

//NOTE: In principle we should also write 'num'.
template<typename T>
EnableIf<std::is_trivially_copyable<T>,
void> serialize(Archive& ar, const T* val, size_t num){
    ar.write(val, num * sizeof(T));
}

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

#endif
