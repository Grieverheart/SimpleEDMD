#ifndef __EVENT_H
#define __EVENT_H

#include "Time.h"
#include <cstddef>

//Will eventually change this to a hashed string
enum EventType{
    EVT_COLLISION
};


//NOTE: Consider using boost::variant, it should give better performance
class Event{
public:
    Event(Time time):
        time_(time)
    {}
    virtual ~Event(void) = 0;
    virtual EventType getType(void)const = 0;

    Time time_;
};
inline Event::~Event(void){}

struct CollisionEvent: public Event{
    CollisionEvent(Time time, size_t particleA, size_t particleB):
        Event(time), pA(particleA), pB(particleB)
    {}

    ~CollisionEvent(void){}

    EventType getType(void)const{
        return EVT_COLLISION;
    }

    size_t pA, pB;
};

struct EventPtrLess{
    bool operator()(Event* a, Event* b)const{
        return a->time_ < b->time_;
    }
};

#endif
