#ifndef __SERIALIZATION_ARCHIVE_H
#define __SERIALIZATION_ARCHIVE_H

#include <cstddef>

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
