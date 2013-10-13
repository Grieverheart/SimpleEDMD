#ifndef __EVENT_H
#define __EVENT_H

#include "Time.h"

typedef int EventType; //Will eventually change this to a hashed string

class EventData{
public:
    virtual ~EventData(void) = 0;
};

struct Event{
    Event(EventType type, Time time, EventData *data = nullptr):
        type_(type), time_(time), data_(data)
    {}
    ~Event(void){
        delete data_;
    }
        
    const EventType  type_;
    const Time       time_;
    const EventData *data_;
};

#endif
