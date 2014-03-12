#ifndef __EVENT_H
#define __EVENT_H

#include <cstddef>

//NOTE: Consider using boost::variant, it should give better performance
enum ParticleEventType{
    PE_NONE = 0,
    PE_COLLISION,
    PE_CELLCROSS
};

class ParticleEvent{
public:
    ParticleEvent(void):
        id_(0)
    {}
    ParticleEvent(double time, int pid, int id, int optional = -1):
        time_(time), pid_(pid), id_(id), optional_(optional)
    {}
    ~ParticleEvent(void){};

    ParticleEventType getType(int npart){
        return (id_ == 0)? PE_NONE: (id_ >= 1 && id_ < npart + 1)? PE_COLLISION: PE_CELLCROSS;
    }

    double time_;
    int pid_;
    int id_, optional_;
};

struct ParticleEventLess{
    bool operator()(const ParticleEvent& a, const ParticleEvent& b)const{
        return a.time_ < b.time_;
    }
};

#endif
