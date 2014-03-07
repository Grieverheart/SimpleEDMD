#ifndef __CELL_LIST_H
#define __CELL_LIST_H

#include <cstddef>
#include <cstdio>
#include "Vec.h"

#define CLL_EMPTY -1

//NOTE: At some point we should make the box a 3x3 matrix

class CellList{
public:
    CellList(void){}

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
            int cellIndices[] = {
                n % nCells_,
               (n / nCells_) % nCells_,
                n / (nCells_ * nCells_)
            };
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

    int getIndex(int pid)const{
        return pCellIds_[pid];
    }

    class Neighbours{
    public:
        Neighbours(const CellList& parent, int pid):
            neighIds_(&parent.cellNeighbours_[27 * parent.pCellIds_[pid]]),
            cell_(parent.cell_), linkedList_(parent.linkedList_),
            curr_nidx(0), curr_pid(cell_[neighIds_[curr_nidx]])
        {}

        Neighbours begin(void)const{
            Neighbours ret = *this;
            ret.curr_nidx = 0;
            ret.curr_pid  = cell_[neighIds_[0]];
            return ret;
        }

        Neighbours end(void)const{
            Neighbours ret = *this;
            ret.curr_nidx = 27;
            return ret;
        }

        bool operator<(const Neighbours& other)const{
            return curr_nidx < other.curr_nidx;
        }

        bool operator>(const Neighbours& other)const{
            return curr_nidx > other.curr_nidx;
        }

        bool operator!=(const Neighbours& other)const{
            return curr_nidx != other.curr_nidx;
        }

        Neighbours& operator++(void){
            if(linkedList_[curr_pid] == CLL_EMPTY){
                ++curr_nidx;
                if(curr_nidx < 27) curr_pid = cell_[neighIds_[curr_nidx]];
            }
            else curr_pid = linkedList_[curr_pid];

            return *this;
        }

        int operator*(void)const{
            return curr_pid;
        }

    private:
        const int* neighIds_;
        const int* cell_;
        const int* linkedList_;
        int curr_nidx;
        int curr_pid;
    };

    Neighbours getNeighbours(int pid)const{
        return Neighbours(*this, pid);
    }

private:
    int cellIndex(const Vec3d& pos)const{
        int cellIndices[3];
        for(int i = 0; i < 3; ++i) cellIndices[i] = int(pos[i] / cellSize_);
        return cellIndices[0] + (cellIndices[1] + cellIndices[2] * nCells_) * nCells_;
    }

private:
    int    nPart_;
    int    nCells_;         //Number of cells in each dimension
    int*   cell_;           //First particle in cell
    int*   linkedList_;
    int*   pCellIds_;
    int*   cellNeighbours_;
    double cellSize_;       //Cell size in each dimension
};

#endif
