#ifndef __BINARY_HEAP_H
#define __BINARY_HEAP_H

#include <vector>

//FUTURE: Add iterators for iterating over tree elements

template<class T, class Compare>
class BinaryHeap{
};

template<class T, class Compare>
class BinaryHeap<T*, Compare>{
public:
    BinaryHeap(void){
        //data_.reserve(6);
        data_.push_back(nullptr);
    }

    ~BinaryHeap(void){
        clear();
    }

    void clear(void){
        for(auto element: data_){
            delete element;
        }
        data_.clear();
    }

    void push(T* val){
        data_.push_back(val);

        int val_index    = data_.size() - 1;
        int parent_index = val_index >> 1;
        while(parent_index > 0 && comp_(data_[val_index], data_[parent_index])){
            swap(val_index, parent_index);
            val_index = parent_index;
            parent_index >>= 1;
        }
    }

    const T* top(void)const{
        return data_[1];
    }

    T* pop(void){
        T* ret_val = data_[1];
        data_[1]   = data_.back();

        int min_index = 1;
        int stop = data_.size() >> 1;
        for(int i = 1; i < stop; i = min_index){
            int left  = i << 1;
            int right = (i << 1) + 1;
            min_index = comp_(data_[left], data_[right])? left: right;
            if(comp_(data_[min_index], data_[i])) swap(min_index, i);
            else break;
        }

        data_.erase(data_.end() - 1);

        return ret_val;
    }

private:
    void swap(int index_a, int index_b){
        T* temp = data_[index_a];
        data_[index_a] = data_[index_b];
        data_[index_b] = temp;
    }

private:
    std::vector<T*> data_;
    static Compare comp_;
};


#endif
