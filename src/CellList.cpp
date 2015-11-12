#include "CellList.h"

//NOTE: At some point we should make the box a 3x3 matrix

CellList::CellList(void):
    cell_(nullptr),
    linkedList_(nullptr),
    pCellIds_(nullptr),
    cellNeighbours_(nullptr),
    cellDirNeighbours_(nullptr),
    cellVolNeighbours_(nullptr)
{}

CellList::~CellList(void){
    delete[] cell_;
    delete[] linkedList_;
    delete[] cellNeighbours_;
    delete[] cellDirNeighbours_;
    delete[] cellVolNeighbours_;
    delete[] pCellIds_;
}

void CellList::init(int nPart, double boxSize, double minCellSize){
    nPart_    = nPart;
    nCells_   = int(boxSize / minCellSize);
    cellSize_ = boxSize / nCells_;

    nCellsTot_ = nCells_ * nCells_ * nCells_;

    cell_              = new int[nCellsTot_];
    linkedList_        = new int[nPart_];
    pCellIds_          = new int[nPart_]{};
    cellNeighbours_    = new int[27 * nCellsTot_];
    cellDirNeighbours_ = new int[54 * nCellsTot_];
    cellVolNeighbours_ = new int[14 * nCellsTot_];


    for(int i = 0; i < nCellsTot_; ++i) cell_[i] = CLL_EMPTY;
    for(int i = 0; i < nPart_; ++i) linkedList_[i] = CLL_EMPTY;

    static const int vol_neigh[] = {11, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
    int k = 0;
    for(int n = 0; n < nCellsTot_; ++n){
        int cellIndices[3];
        index_to_indices(n, cellIndices);
        int ds[6] = {0};
        for(int i = 0; i < 27; ++i){
            int offset[] = {i % 3 - 1, (i / 3) % 3 - 1, i / 9 - 1};
            int cell_idx =  (nCells_ + cellIndices[0] + offset[0]) % nCells_ +
                           ((nCells_ + cellIndices[1] + offset[1]) % nCells_ +
                           ((nCells_ + cellIndices[2] + offset[2]) % nCells_) * nCells_) * nCells_;
            for(int d = 0; d < 6; ++d){
                if(offset[d / 2] == 1 - 2 * (d % 2)){
                    cellDirNeighbours_[9 * (6 * n + d) + ds[d]++] = cell_idx;
                }
            }
            cellNeighbours_[27 * n + i] = cell_idx;
            if(i == vol_neigh[k % 14]) cellVolNeighbours_[k++] = cell_idx;
        }
    }
}

int CellList::add(int pid, const clam::Vec3d& pos){
    int cidx       = cell_index(pos);
    pCellIds_[pid] = cidx;

    linkedList_[pid] = cell_[cidx];
    cell_[cidx] = pid;

    return cidx;
}

int CellList::update(int pid, const clam::Vec3d& pos){
    int cidx = cell_index(pos);
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
    index_to_indices(cidx_old, cellIndices);
    int dir = coffset / 2;
    cellIndices[dir] = (cellIndices[dir] + nCells_ + 2 * (coffset % 2) - 1) % nCells_;
    int cidx = indices_to_index(cellIndices);

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

int CellList::cell_index(int pid)const{
    return pCellIds_[pid];
}

clam::Vec3d CellList::cell_origin(int cidx)const{
    clam::Vec3d cellIndices = clam::Vec3d(
        cidx % nCells_,
       (cidx / nCells_) % nCells_,
        cidx / (nCells_ * nCells_)
    );
    return cellIndices * cellSize_;
}

clam::Vec3d CellList::cell_size(void)const{
    return clam::Vec3d(cellSize_);
}

CellList::CellIterator CellList::cells(void)const{
    return CellIterator(nCellsTot_);
}

CellList::CellNeighbourIterator CellList::cell_nbs(int cid)const{
    return CellNeighbourIterator(cellNeighbours_ + 27 * cid);
}

CellList::CellNeighbourIterator CellList::cell_vol_nbs(int cid)const{
    return CellNeighbourIterator(cellVolNeighbours_ + 14 * cid, 14);
}

CellList::CellNeighbourIterator CellList::particle_cell_nbs(int pid)const{
    return CellNeighbourIterator(cellNeighbours_ + 27 * pCellIds_[pid]);
}

CellList::CellNeighbourIterator CellList::particle_cell_nbs(const clam::Vec3d& pos)const{
    return CellNeighbourIterator(cellNeighbours_ + 27 * cell_index(pos));
}

CellList::CellNeighbourIterator CellList::cell_dir_nbs(int cid, int dir)const{
    return CellNeighbourIterator(cellDirNeighbours_ + 54 * cid + 9 * dir, 9);
}

CellList::CellContentIterator CellList::cell_content(int cid)const{
    return CellContentIterator(*this, cid);
}

int CellList::indices_to_index(const int (&indices)[3])const{
    return indices[0] + (indices[1] + indices[2] * nCells_) * nCells_;
}

void CellList::index_to_indices(int index, int (&indices)[3])const{
    indices[0] =  index % nCells_;
    indices[1] = (index / nCells_) % nCells_;
    indices[2] =  index / (nCells_ * nCells_);
}

int CellList::cell_index(const clam::Vec3d& pos)const{
    int cellIndices[3];
    for(int i = 0; i < 3; ++i) cellIndices[i] = int(pos[i] / cellSize_);
    return indices_to_index(cellIndices);
}

