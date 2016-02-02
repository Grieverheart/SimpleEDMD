#include "serialization/archive.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>

Archive::Archive(void):
    data_(nullptr),
    position_(0), size_(0)
{}

Archive::Archive(const void* data, size_t size):
    position_(0), size_(size)
{
    data_ = reinterpret_cast<char*>(malloc(size));
    memcpy(data_, data, size);
}

Archive::~Archive(void){
    free(data_);
}

void Archive::write(const void* data, size_t size){
    if(position_ + size > size_){
        size_ = (size_ > 0)? (position_ + size) * 2: sizeof(double);
        data_ = reinterpret_cast<char*>(realloc(data_, size_));
    }
    memcpy(data_ + position_, data, size);
    position_ += size;
}

//assert(position_ + size < size_)
void Archive::read(void* data, size_t size){
    memcpy(data, data_ + position_, size);
    position_ += size;
}

const char* Archive::data(void)const{
    return data_;
}

size_t Archive::size(void)const{
    return position_;
}
