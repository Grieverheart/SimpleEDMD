#include "CellList.h"

//NOTE: At some point we should make the box a 3x3 matrix

CellList::CellList(void){}

CellList::~CellList(void){
    delete[] cell_;
    delete[] linkedList_;
    delete[] cellNeighbours_;
    delete[] pCellIds_;
}

void CellList::init(int nPart, double boxSize, double minCellSize){
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

int CellList::add(int pid, const clam::Vec3d& pos){
    int cidx       = cellIndex(pos);
    pCellIds_[pid] = cidx;

    linkedList_[pid] = cell_[cidx];
    cell_[cidx] = pid;

    return cidx;
}

int CellList::update(int pid, const clam::Vec3d& pos){
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

int CellList::move(int pid, int coffset){
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

int CellList::getIndex(int pid)const{
    return pCellIds_[pid];
}

clam::Vec3d CellList::getCellOrigin(int cidx)const{
    clam::Vec3d cellIndices = clam::Vec3d(
        cidx % nCells_,
       (cidx / nCells_) % nCells_,
        cidx / (nCells_ * nCells_)
    );
    return cellIndices * cellSize_;
}

clam::Vec3d CellList::getCellSize(void)const{
    return clam::Vec3d(cellSize_);
}

CellList::NeighbourIterator CellList::getNeighbourIterator(int pid)const{
    return NeighbourIterator(*this, pid);
}

CellList::DirectionalNeighbourIterator CellList::getDirNeighbourIterator(int pid, int dir)const{
    return DirectionalNeighbourIterator(*this, pid, dir);
}

CellList::CellIterator CellList::getCellIterator(int cid)const{
    return CellIterator(*this, cid);
}

int CellList::indicesToIndex(const int (&indices)[3])const{
    return indices[0] + (indices[1] + indices[2] * nCells_) * nCells_;
}

void CellList::indexToIndices(int index, int (&indices)[3])const{
    indices[0] =  index % nCells_;
    indices[1] = (index / nCells_) % nCells_;
    indices[2] =  index / (nCells_ * nCells_);
}

int CellList::cellIndex(const clam::Vec3d& pos)const{
    int cellIndices[3];
    for(int i = 0; i < 3; ++i) cellIndices[i] = int(pos[i] / cellSize_);
    return indicesToIndex(cellIndices);
}

