#ifndef __ARCHIVE_H
#define __ARCHIVE_H

#include <cstddef>

//TODO: We can make this a bit smarter. I.e. we can provide a function to write vector to avoid boilerplate:
//auto vec_size = vec.size();
//ar.write(&vec_size, sizeof(vec_size));
//ar.write(vec.data(), vec_size * sizeof(decltype(vec)::value_type));
//But this would only work for POD value types. Instead, it would be nice if we had a separate serialize function:
//template<typename T>
//void serialize(const T&, Archive&);
//For classes we declare it as friend and continue as earlier. For pods it's simple.

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

#endif
