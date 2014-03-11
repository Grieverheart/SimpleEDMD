#ifndef __EVENT_H
#define __EVENT_H

#include <cstddef>

//Will eventually change this to a hashed string
enum EventType{
    EVT_COLLISION,
    EVT_CELLCROSS
};


//NOTE: Consider using boost::variant, it should give better performance
class Event{
public:
    Event(double time):
        time_(time)
    {}
    virtual ~Event(void) = 0;
    virtual EventType getType(void)const = 0;

    double time_;
};
inline Event::~Event(void){}

struct CollisionEvent: public Event{
    CollisionEvent(double time, size_t particleA, size_t particleB, size_t numCollisionsB):
        Event(time), pA(particleA), pB(particleB), nBCollisions(numCollisionsB)
    {}

    ~CollisionEvent(void){}

    EventType getType(void)const{
        return EVT_COLLISION;
    }

    size_t pA, pB;
    size_t nBCollisions;
};

struct CellCrossEvent: public Event{
    CellCrossEvent(double time, size_t particleID, int cellOffset):
        Event(time), pid(particleID), coffset(cellOffset)
    {}

    ~CellCrossEvent(void){}

    EventType getType(void)const{
        return EVT_CELLCROSS;
    }

    size_t pid;
    int coffset;
};

struct EventPtrLess{
    bool operator()(Event* a, Event* b)const{
        return a->time_ < b->time_;
    }
};

#endif
