#ifndef __BINARY_HEAP_H
#define __BINARY_HEAP_H

#include <vector>
#include <functional>
#include "serialization/archive.h"

//FUTURE: Add iterators for iterating over tree elements

template< class T, class Compare = std::less<T> >
class BinaryHeap;

template<typename T>
void serialize(Archive&, const BinaryHeap<T>&);

template< class T, class Compare >
class BinaryHeap{
public:
    BinaryHeap(void){
        data_.push_back(T());
    }

    ~BinaryHeap(void){
        clear();
    }
    
    BinaryHeap(const BinaryHeap& other) = delete;

    BinaryHeap(BinaryHeap&& other):
        data_(std::move(other.data_))
    {}

    friend void serialize<T>(Archive& ar, const BinaryHeap<T>&);

    void clear(void){
        data_.clear();
        data_.push_back(T());
    }

    size_t size(void){
        return data_.size() - 1;
    }

    bool empty(void)const{
        return (data_.size() < 2);
    }

    void push(const T& val){
        data_.push_back(val);

        int val_index    = data_.size() - 1;
        int parent_index = val_index >> 1;
        while(parent_index > 0 && comp_(data_[val_index], data_[parent_index])){
            std::swap(data_[val_index], data_[parent_index]);
            val_index = parent_index;
            parent_index >>= 1;
        }
    }

    const T& top(void)const{
        return data_[1];
    }

    T pop(void){
        T ret_val = data_[1];
        data_[1]  = data_.back();

        int min_index = 1;
        int stop = data_.size() >> 1;
        for(int i = 1; i < stop; i = min_index){
            int left  = i << 1;
            int right = (i << 1) + 1;
            min_index = comp_(data_[left], data_[right])? left: right;
            if(comp_(data_[min_index], data_[i])) std::swap(data_[min_index], data_[i]);
            else break;
        }

        data_.erase(data_.end() - 1);

        return ret_val;
    }

private:
    std::vector<T> data_;
    static Compare comp_;
};

template<class T, class Compare>
Compare BinaryHeap<T, Compare>::comp_ = Compare();

template<typename T>
void serialize(Archive& ar, const BinaryHeap<T>& bh){
    serialize(ar, bh.data_);
}


#endif
