#ifndef __CELL_LIST_H
#define __CELL_LIST_H

#include <cstddef>
#include "Vec.h"

#define CLL_EMPTY -1

constexpr bool isDirNeighbour(int n, int dir){
    return (dir / 2 == 0)? ((n % 3 - 1) * (2 * (dir % 2) - 1) > 0):
           (dir / 2 == 1)? (((n / 3) % 3 - 1) * (2 * (dir % 2) - 1) > 0):
                           (((n / 9) - 1) * (2 * (dir % 2) - 1) > 0);
}

//NOTE: At some point we should make the box a 3x3 matrix
class CellList{
public:
    CellList(void){
        for(int i = 0; i < 6; ++i){
            int nids = 0;
            for(int n = 0; n < 27; ++n){
                if(isDirNeighbour(n, i)) dirNeighbourIds_[9 * i + nids++] = n;
            }
        }
    }

    ~CellList(void){
        delete[] cell_;
        delete[] linkedList_;
        delete[] cellNeighbours_;
        delete[] pCellIds_;
    }

    void init(int nPart, double boxSize, double minCellSize){
        nPart_    = nPart;
        nCells_   = int(boxSize / minCellSize);
        cellSize_ = boxSize / nCells_;

        int nCellsTot = nCells_ * nCells_ * nCells_;

        cell_           = new int[nCellsTot];
        linkedList_     = new int[nPart_];
        pCellIds_       = new int[nPart_]{};
        cellNeighbours_ = new int[27 * nCellsTot];


        for(int i = 0; i < nCellsTot; ++i) cell_[i] = CLL_EMPTY;
        for(int i = 0; i < nPart_; ++i) linkedList_[i] = CLL_EMPTY;

        for(int n = 0; n < nCellsTot; ++n){
            int cellIndices[3];
            indexToIndices(n, cellIndices);
            for(int i = 0; i < 27; ++i){
                int offset[] = {i % 3 - 1, (i / 3) % 3 - 1, i / 9 - 1};
                int cell_idx =  (nCells_ + cellIndices[0] + offset[0]) % nCells_ +
                               ((nCells_ + cellIndices[1] + offset[1]) % nCells_ +
                               ((nCells_ + cellIndices[2] + offset[2]) % nCells_) * nCells_) * nCells_;
                cellNeighbours_[27 * n + i] = cell_idx;
            }
        }
    }

    int add(int pid, const Vec3d& pos){
        int cidx       = cellIndex(pos);
        pCellIds_[pid] = cidx;

        linkedList_[pid] = cell_[cidx];
        cell_[cidx] = pid;

        return cidx;
    }

    int update(int pid, const Vec3d& pos){
        int cidx = cellIndex(pos);
        if(cidx == pCellIds_[pid]) return cidx;

        int cidx_old = pCellIds_[pid];
        if(cell_[cidx_old] == pid) cell_[cidx_old] = linkedList_[pid];
        else{
            for(int ci = cell_[cidx_old]; ci != CLL_EMPTY; ci = linkedList_[ci]){
                if(linkedList_[ci] == pid){
                    linkedList_[ci] = linkedList_[pid];
                    break;
                }
            }
        }

        linkedList_[pid] = cell_[cidx];
        cell_[cidx] = pid;
        pCellIds_[pid] = cidx;

        return cidx;
    }

    int move(int pid, int coffset){
        int cidx_old = pCellIds_[pid];
        int cellIndices[3];
        indexToIndices(cidx_old, cellIndices);
        int dir = coffset / 2;
        cellIndices[dir] = (cellIndices[dir] + nCells_ + 2 * (coffset % 2) - 1) % nCells_;
        int cidx = indicesToIndex(cellIndices);

        if(cell_[cidx_old] == pid) cell_[cidx_old] = linkedList_[pid];
        else{
            for(int ci = cell_[cidx_old]; ci != CLL_EMPTY; ci = linkedList_[ci]){
                if(linkedList_[ci] == pid){
                    linkedList_[ci] = linkedList_[pid];
                    break;
                }
            }
        }

        linkedList_[pid] = cell_[cidx];
        cell_[cidx] = pid;
        pCellIds_[pid] = cidx;

        return cidx;
    }

    int getIndex(int pid)const{
        return pCellIds_[pid];
    }

    Vec3d getCellOrigin(int cidx)const{
        Vec3d cellIndices = Vec3d(
            cidx % nCells_,
           (cidx / nCells_) % nCells_,
            cidx / (nCells_ * nCells_)
        );
        return cellIndices * cellSize_;
    }

    Vec3d getCellSize(void)const{
        return Vec3d(cellSize_);
    }

    class NeighbourIterator{
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

    class DirectionalNeighbourIterator{
    public:
        DirectionalNeighbourIterator(const CellList& parent, int pid, int dir):
            dirNeighbourIds_(&parent.dirNeighbourIds_[9 * dir]),
            cellNeighbours_(&parent.cellNeighbours_[27 * parent.pCellIds_[pid]]),
            nidx_(0)
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
            return cellNeighbours_[dirNeighbourIds_[nidx_]];
        }
    private:
        const int* dirNeighbourIds_;
        const int* cellNeighbours_;
        int nidx_;
    };

    class CellIterator{
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

    NeighbourIterator getNeighbourIterator(int pid)const{
        return NeighbourIterator(*this, pid);
    }

    DirectionalNeighbourIterator getDirNeighbourIterator(int pid, int dir)const{
        return DirectionalNeighbourIterator(*this, pid, dir);
    }

    CellIterator getCellIterator(int cid)const{
        return CellIterator(*this, cid);
    }

private:
    int indicesToIndex(const int (&indices)[3])const{
        return indices[0] + (indices[1] + indices[2] * nCells_) * nCells_;
    }
    void indexToIndices(int index, int (&indices)[3])const{
        indices[0] =  index % nCells_;
        indices[1] = (index / nCells_) % nCells_;
        indices[2] =  index / (nCells_ * nCells_);
    }
    int cellIndex(const Vec3d& pos)const{
        int cellIndices[3];
        for(int i = 0; i < 3; ++i) cellIndices[i] = int(pos[i] / cellSize_);
        return indicesToIndex(cellIndices);
    }

private:
    int    nPart_;
    int    nCells_;         //Number of cells in each dimension
    int*   cell_;           //First particle in cell
    int*   linkedList_;
    int*   pCellIds_;
    int*   cellNeighbours_;
    int    dirNeighbourIds_[54];
    double cellSize_;       //Cell size in each dimension
};

#endif
