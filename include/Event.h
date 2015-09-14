#ifndef __EVENT_H
#define __EVENT_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include "clam.h"

enum ParticleEventType{
    PE_NONE               = 0,
    PE_COLLISION          = 1,
    PE_POSSIBLE_COLLISION = 2,
    PE_CELLCROSS          = 3
};

class ParticleEvent{
private:
    ParticleEvent(ParticleEventType type, double time, int32_t pid, int32_t id):
        time_(time), pid_(pid), id_((type << 29) | id)
    {}

public:
    ParticleEvent(void):
        id_(0)
    {}

    static inline ParticleEvent None(void){
        return ParticleEvent();
    }

    static inline ParticleEvent Collision(double time, int pid, int id, int optional, clam::Vec3d normal){
        ParticleEvent event(PE_COLLISION, time, pid, id);
        event.optional_ = optional;
        event.normal_ = normal;
        return event;
    }

    static inline ParticleEvent PossibleCollision(double time, int pid, int id, int optional){
        ParticleEvent event(PE_POSSIBLE_COLLISION, time, pid, id);
        event.optional_ = optional;
        return event;
    }

    static inline ParticleEvent CellCross(double time, int pid, int id){
        return ParticleEvent(PE_CELLCROSS, time, pid, id);
    }

    int get_type(void){
        return (id_ >> 29);
    }

    int32_t get_id(void)const{
        return (id_ & 0x1FFFFFFF);
    }

    bool operator<(const ParticleEvent& b)const{
        return time_ < b.time_;
    }

    //Essential data
    double time_;
    uint32_t pid_, id_;

    //Data needed for collision events.
    uint32_t optional_;
    //TODO: Make this a pointer for saving space
    clam::Vec3d normal_;
};

#endif
