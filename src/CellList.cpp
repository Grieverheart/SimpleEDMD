#include "CellList.h"

//NOTE: At some point we should make the box a 3x3 matrix
template<class T> 
static inline T max(T a, T b){
    return (a < b) ? b : a;
}

CellList::CellList(void):
    cell_(nullptr),
    linked_list_(nullptr),
    p_cell_ids_(nullptr),
    cell_neighbours_(nullptr),
    cell_dir_neighbours_(nullptr),
    cell_vol_neighbours_(nullptr)
{}

CellList::~CellList(void){
    delete[] cell_;
    delete[] linked_list_;
    delete[] cell_neighbours_;
    delete[] cell_dir_neighbours_;
    delete[] cell_vol_neighbours_;
    delete[] p_cell_ids_;
}

void serialize(Archive& ar, const CellList& cll){
    serialize(ar, cll.n_part_);
    serialize(ar, cll.n_cells_, 3);
    serialize(ar, cll.n_cells_tot_);
    serialize(ar, cll.cell_, cll.n_cells_tot_);
    serialize(ar, cll.linked_list_, cll.n_part_);
    serialize(ar, cll.p_cell_ids_, cll.n_part_);
    serialize(ar, cll.cell_neighbours_, 27 * cll.n_cells_tot_);
    serialize(ar, cll.cell_dir_neighbours_, 54 * cll.n_cells_tot_);
    serialize(ar, cll.cell_vol_neighbours_, 14 * cll.n_cells_tot_);
    serialize(ar, cll.cell_size_);
}

void CellList::init(int nPart, const clam::Vec3d& box_bounds, double minCellSize){
    n_part_ = nPart;

    n_cells_tot_ = 1;
    auto nCells = box_bounds / minCellSize;
    for(int i = 0; i < 3; ++i){
        // Handle the case where n_cells_ < 3. We just need to make it 3!
        n_cells_[i]   = max(3, int(nCells[i]));
        cell_size_[i] = box_bounds[i] / n_cells_[i];
        n_cells_tot_ *= n_cells_[i];
    }

    cell_              = new int[n_cells_tot_];
    linked_list_        = new int[n_part_];
    p_cell_ids_          = new int[n_part_]{};
    cell_neighbours_    = new int[27 * n_cells_tot_];
    cell_dir_neighbours_ = new int[54 * n_cells_tot_];
    cell_vol_neighbours_ = new int[14 * n_cells_tot_];


    for(int i = 0; i < n_cells_tot_; ++i) cell_[i] = CLL_EMPTY;
    for(int i = 0; i < n_part_; ++i) linked_list_[i] = CLL_EMPTY;

    static const int vol_neigh[] = {11, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
    int k = 0;
    for(int n = 0; n < n_cells_tot_; ++n){
        int cellIndices[3];
        index_to_indices(n, cellIndices);
        int ds[6] = {0};
        for(int i = 0; i < 27; ++i){
            int offset[] = {i % 3 - 1, (i / 3) % 3 - 1, i / 9 - 1};
            int cell_idx =  (n_cells_[0] + cellIndices[0] + offset[0]) % n_cells_[0] +
                           ((n_cells_[1] + cellIndices[1] + offset[1]) % n_cells_[1] +
                           ((n_cells_[2] + cellIndices[2] + offset[2]) % n_cells_[2]) * n_cells_[1]) * n_cells_[0];
            for(int d = 0; d < 6; ++d){
                if(offset[d / 2] == 1 - 2 * (d % 2)){
                    cell_dir_neighbours_[9 * (6 * n + d) + ds[d]++] = cell_idx;
                }
            }
            cell_neighbours_[27 * n + i] = cell_idx;
            if(i == vol_neigh[k % 14]) cell_vol_neighbours_[k++] = cell_idx;
        }
    }
}

int CellList::add(int pid, const clam::Vec3d& pos){
    int cidx       = cell_index(pos);
    p_cell_ids_[pid] = cidx;

    linked_list_[pid] = cell_[cidx];
    cell_[cidx] = pid;

    return cidx;
}

int CellList::update(int pid, const clam::Vec3d& pos){
    int cidx = cell_index(pos);
    if(cidx == p_cell_ids_[pid]) return cidx;

    int cidx_old = p_cell_ids_[pid];
    if(cell_[cidx_old] == pid) cell_[cidx_old] = linked_list_[pid];
    else{
        for(int ci = cell_[cidx_old]; ci != CLL_EMPTY; ci = linked_list_[ci]){
            if(linked_list_[ci] == pid){
                linked_list_[ci] = linked_list_[pid];
                break;
            }
        }
    }

    linked_list_[pid] = cell_[cidx];
    cell_[cidx] = pid;
    p_cell_ids_[pid] = cidx;

    return cidx;
}

int CellList::move(int pid, int coffset){
    int cidx_old = p_cell_ids_[pid];
    int cellIndices[3];
    index_to_indices(cidx_old, cellIndices);
    int dir = coffset / 2;
    cellIndices[dir] = (cellIndices[dir] + n_cells_[dir] + 2 * (coffset % 2) - 1) % n_cells_[dir];
    int cidx = indices_to_index(cellIndices);

    if(cell_[cidx_old] == pid) cell_[cidx_old] = linked_list_[pid];
    else{
        for(int ci = cell_[cidx_old]; ci != CLL_EMPTY; ci = linked_list_[ci]){
            if(linked_list_[ci] == pid){
                linked_list_[ci] = linked_list_[pid];
                break;
            }
        }
    }

    linked_list_[pid] = cell_[cidx];
    cell_[cidx] = pid;
    p_cell_ids_[pid] = cidx;

    return cidx;
}

int CellList::cell_index(int pid)const{
    return p_cell_ids_[pid];
}

clam::Vec3d CellList::cell_origin(int cidx)const{
    return cell_size_ * clam::Vec3d(
        cidx % n_cells_[0],
       (cidx / n_cells_[0]) % n_cells_[1],
        cidx / (n_cells_[0] * n_cells_[1])
    );
}

clam::Vec3d CellList::cell_size(void)const{
    return clam::Vec3d(cell_size_);
}

CellList::CellIterator CellList::cells(void)const{
    return CellIterator(n_cells_tot_);
}

CellList::CellNeighbourIterator CellList::cell_nbs(int cid)const{
    return CellNeighbourIterator(cell_neighbours_ + 27 * cid);
}

CellList::CellNeighbourIterator CellList::cell_vol_nbs(int cid)const{
    return CellNeighbourIterator(cell_vol_neighbours_ + 14 * cid, 14);
}

CellList::CellNeighbourIterator CellList::particle_cell_nbs(int pid)const{
    return CellNeighbourIterator(cell_neighbours_ + 27 * p_cell_ids_[pid]);
}

CellList::CellNeighbourIterator CellList::particle_cell_nbs(const clam::Vec3d& pos)const{
    return CellNeighbourIterator(cell_neighbours_ + 27 * cell_index(pos));
}

CellList::CellNeighbourIterator CellList::cell_dir_nbs(int cid, int dir)const{
    return CellNeighbourIterator(cell_dir_neighbours_ + 54 * cid + 9 * dir, 9);
}

CellList::CellContentIterator CellList::cell_content(int cid)const{
    return CellContentIterator(*this, cid);
}

int CellList::indices_to_index(const int (&indices)[3])const{
    return indices[0] + (indices[1] + indices[2] * n_cells_[1]) * n_cells_[0];
}

void CellList::index_to_indices(int index, int (&indices)[3])const{
    indices[0] =  index % n_cells_[0];
    indices[1] = (index / n_cells_[0]) % n_cells_[1];
    indices[2] =  index / (n_cells_[0] * n_cells_[1]);
}

int CellList::cell_index(const clam::Vec3d& pos)const{
    int cellIndices[3];
    for(int i = 0; i < 3; ++i) cellIndices[i] = int(pos[i] / cell_size_[i]);
    return indices_to_index(cellIndices);
}

