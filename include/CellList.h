#ifndef __CELL_LIST_H
#define __CELL_LIST_H

#include "Vec.h"

#define CLL_EMPTY -1

//NOTE: At some point we should make the box a 3x3 matrix
class CellList{
public:
    CellList(void);
    ~CellList(void);

    void init(int nPart, double boxSize, double minCellSize);
    int add(int pid, const Vec3d& pos);
    int update(int pid, const Vec3d& pos);
    int move(int pid, int coffset);

    int getIndex(int pid)const;
    Vec3d getCellOrigin(int cidx)const;
    Vec3d getCellSize(void)const;

    class NeighbourIterator;
    class DirectionalNeighbourIterator;
    class CellIterator;

    NeighbourIterator getNeighbourIterator(int pid)const;
    DirectionalNeighbourIterator getDirNeighbourIterator(int pid, int dir)const;
    CellIterator getCellIterator(int cid)const;

private:
    int indicesToIndex(const int (&indices)[3])const;
    void indexToIndices(int index, int (&indices)[3])const;
    int cellIndex(const Vec3d& pos)const;

private:
    int    nPart_;
    int    nCells_;         //Number of cells in each dimension
    int*   cell_;           //First particle in cell
    int*   linkedList_;
    int*   pCellIds_;
    int*   cellNeighbours_;
    double cellSize_;       //Cell size in each dimension
};

class CellList::NeighbourIterator{
public:
    NeighbourIterator(const CellList& parent, int pid):
        cellNeighbours_(&parent.cellNeighbours_[27 * parent.pCellIds_[pid]]),
        nidx_(0)
    {}

    NeighbourIterator begin(void){
        NeighbourIterator ret = *this;
        ret.nidx_ = 0;
        return ret;
    }

    NeighbourIterator end(void){
        NeighbourIterator ret = *this;
        ret.nidx_ = 27;
        return ret;
    }

    bool operator!=(const NeighbourIterator& other){
        return nidx_ != other.nidx_;
    }

    NeighbourIterator& operator++(void){
        ++nidx_;
        return *this;
    }

    int operator*(void)const{
        return cellNeighbours_[nidx_];
    }
private:
    const int* cellNeighbours_;
    int nidx_;
};


//Compile-time Directional Neighbour index generator
namespace{
    constexpr bool isDirNeighbour(int n){
        return ((n / 27) / 2 == 0)?  (((n % 27) % 3 - 1) * (2 * ((n / 27) % 2) - 1) > 0):
               ((n / 27) / 2 == 1)? ((((n % 27) / 3) % 3 - 1) * (2 * ((n / 27) % 2) - 1) > 0):
                                    ((((n % 27) / 9) - 1) * (2 * ((n / 27) % 2) - 1) > 0);
    }

    //Decides to keep N % 27 or not
    template<bool B, int N, int...S>
    struct dummy;

    template<int N, int...S>
    struct gen_nums: dummy<isDirNeighbour(N - 1), (N - 1), S...>{
        static_assert(N <= 162, "");
    };

    template<int...S>
    struct gen_nums<0, S...>{
        static constexpr int value[] = {S...};
    };

    template<int...S>
    constexpr int gen_nums<0, S...>::value[];

    template<int N, int...S>
    struct dummy<false, N, S...>: gen_nums<N, S...>{};

    template<int N, int...S>
    struct dummy<true, N, S...>: gen_nums<N, N % 27, S...>{};

    struct DirNeighbours: gen_nums<27 * 6>{};
}

class CellList::DirectionalNeighbourIterator{
public:
    DirectionalNeighbourIterator(const CellList& parent, int pid, int dir):
        cellNeighbours_(&parent.cellNeighbours_[27 * parent.pCellIds_[pid]]),
        nidx_(0), dir_(dir)
    {}

    DirectionalNeighbourIterator begin(void){
        DirectionalNeighbourIterator ret = *this;
        ret.nidx_ = 0;
        return ret;
    }

    DirectionalNeighbourIterator end(void){
        DirectionalNeighbourIterator ret = *this;
        ret.nidx_ = 9;
        return ret;
    }

    bool operator!=(const DirectionalNeighbourIterator& other)const{
        return nidx_ != other.nidx_;
    }

    DirectionalNeighbourIterator& operator++(void){
        ++nidx_;
        return *this;
    }

    int operator*(void)const{
        return cellNeighbours_[DirNeighbours::value[9 * dir_ + nidx_]];
    }
private:
    const int* cellNeighbours_;
    int nidx_, dir_;
};

class CellList::CellIterator{
public:
    CellIterator(const CellList& parent, int cidx):
        first_pid(parent.cell_[cidx]), curr_pid(first_pid),
        linkedList_(parent.linkedList_)
    {}

    CellIterator begin(void){
        CellIterator ret = *this;
        ret.curr_pid = first_pid;
        return ret;
    }

    CellIterator end(void){
        CellIterator ret = *this;
        ret.curr_pid = CLL_EMPTY;
        return ret;
    }

    bool operator!=(const CellIterator& other)const{
        return curr_pid != other.curr_pid;
    }

    CellIterator& operator++(void){
        curr_pid = linkedList_[curr_pid];
        return *this;
    }

    int operator*(void)const{
        return curr_pid;
    }
private:
    int first_pid;
    int curr_pid;
    const int* linkedList_;
};

#endif
