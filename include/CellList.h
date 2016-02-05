#ifndef __CELL_LIST_H
#define __CELL_LIST_H

#include "clam.h"
#include "serialization/archive.h"

#define CLL_EMPTY -1

//NOTE: At some point we should make the box a 3x3 matrix
class CellList{
public:
    CellList(void);
    ~CellList(void);

    void serialize(Archive&)const;

    void init(int nPart, const clam::Vec3d& boxSize, double minCellSize);
    int add(int pid, const clam::Vec3d& pos);
    int update(int pid, const clam::Vec3d& pos);
    int move(int pid, int coffset);

    int cell_index(int pid)const;
    clam::Vec3d cell_origin(int cidx)const;
    clam::Vec3d cell_size(void)const;

    class CellIterator;
    class CellNeighbourIterator;
    class CellContentIterator;

    CellIterator          cells(void)const;
    CellNeighbourIterator cell_vol_nbs(int pid)const;
    CellNeighbourIterator cell_nbs(int cid)const;
    CellNeighbourIterator particle_cell_nbs(int pid)const;
    CellNeighbourIterator particle_cell_nbs(const clam::Vec3d& pos)const;
    CellNeighbourIterator cell_dir_nbs(int cid, int dir)const;
    CellContentIterator   cell_content(int cid)const;

private:
    int indices_to_index(const int (&indices)[3])const;
    void index_to_indices(int index, int (&indices)[3])const;
    int cell_index(const clam::Vec3d& pos)const;

private:
    int    n_part_;
    int    n_cells_[3]; //Number of cells in each dimension
    int    n_cells_tot_;
    int*   cell_; //First particle in cell
    int*   linked_list_;
    int*   p_cell_ids_;
    int*   cell_neighbours_;
    int*   cell_dir_neighbours_;
    int*   cell_vol_neighbours_;
    clam::Vec3d cell_size_; //Cell size in each dimension
};

class CellList::CellNeighbourIterator{
    friend class CellList;
public:
    constexpr CellNeighbourIterator(const int* neighbour_list, int n_neighbours = 27):
        cell_neighbours_(neighbour_list),
        nidx_(0), end_(n_neighbours)
    {}

    constexpr CellNeighbourIterator begin(void)const{
        return CellNeighbourIterator(cell_neighbours_, end_);
    }

    CellNeighbourIterator end(void)const{
        CellNeighbourIterator ret = *this;
        ret.nidx_ = end_;
        return ret;
    }

    constexpr bool operator!=(const CellNeighbourIterator& other)const{
        return nidx_ != other.nidx_;
    }

    CellNeighbourIterator& operator++(void){
        ++nidx_;
        return *this;
    }

    constexpr int operator*(void)const{
        return cell_neighbours_[nidx_];
    }
private:
    const int* cell_neighbours_;
    int nidx_;
    int end_;
};

class CellList::CellIterator{
public:
    constexpr CellIterator(int ncells):
        ncells_(ncells), cid_(0)
    {}

    constexpr CellIterator begin(void)const{
        return CellIterator(ncells_);
    }

    CellIterator end(void)const{
        CellIterator ret = *this;
        ret.cid_ = ncells_;
        return ret;
    }

    constexpr bool operator!=(const CellIterator& other)const{
        return cid_ != other.cid_;
    }

    CellIterator& operator++(void){
        ++cid_;
        return *this;
    }

    constexpr int operator*(void)const{
        return cid_;
    }
private:
    int ncells_;
    int cid_;
};

class CellList::CellContentIterator{
public:
    constexpr CellContentIterator(const CellList& parent, int cidx):
        first_pid(parent.cell_[cidx]), curr_pid(parent.cell_[cidx]),
        linked_list_(parent.linked_list_)
    {}

    CellContentIterator begin(void)const{
        CellContentIterator ret = *this;
        ret.curr_pid = first_pid;
        return ret;
    }

    CellContentIterator end(void)const{
        CellContentIterator ret = *this;
        ret.curr_pid = CLL_EMPTY;
        return ret;
    }

    constexpr bool operator!=(const CellContentIterator& other)const{
        return curr_pid != other.curr_pid;
    }

    CellContentIterator& operator++(void){
        curr_pid = linked_list_[curr_pid];
        return *this;
    }

    CellContentIterator operator+(int num)const{
        auto ret = *this;
        for(int i = 0; i < num; ++i) ret.curr_pid = ret.linked_list_[ret.curr_pid];
        return ret;
    }

    constexpr int operator*(void)const{
        return curr_pid;
    }
private:
    int first_pid;
    int curr_pid;
    const int* linked_list_;
};

#endif
