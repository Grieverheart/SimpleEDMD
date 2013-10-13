#ifndef __EVENT_MANAGER_H
#define __EVENT_MANAGER_H

#include "Event.h"

typedef int EventRef;

struct EventItem;

//Consider collecting statistics to improve estimates of scaleFactor and llSize after clears
class EventManager{
public:
    explicit EventManager(size_t nMaxEvents, double scaleFactor, int llSize); //Consider making these template parameters
    ~EventManager(void);

    void resize(size_t newSize);

    EventRef queueEvent(Event* event);
    Event*   getNextEvent(void);
    void     deleteEvent(EventRef ref);
    void     clear(void);
//private Functions
private:
    void cbtUpdate(EventRef eRef);
    void cbtDelete(EventRef eRef);
    void cbtInsert(EventRef eRef);
//Members
private:
    //Consider implementing these as fixed size and picking an id from an id_manager
    //NEW: Even better, use the id_manager, and when no id is available push_back
    std::vector<Event>     events_;
    std::vector<EventItem> eventItems_; //Same size as events_. In principle it is a linked-list/binary tree node.

    size_t size_; //The number of events to hold. Should be at least equal to the number of particles.

    //Priority queue vars
    std::vector<EventRef>  llQueue_;    //The array of linked-lists

    int    currentIndex_;
    int    baseIndex_;

    int    llSize_;
    double scaleFactor_;

    //Complete binary tree vars
    std::vector<EventRef> nodes_;
    std::vector<size_t>   leafs_;

    int nEvents_; //Current number of events in binary tree
};

#endif
