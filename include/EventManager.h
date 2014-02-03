#ifndef __EVENT_MANAGER_H
#define __EVENT_MANAGER_H

#include "Event.h"
#include <set>
#include <vector>
#include <cstddef>

typedef int EventRef;

struct EventItem;

//Consider collecting statistics to improve estimates of scaleFactor and llSize after clears
class EventManager{
public:
    explicit EventManager(double scaleFactor, int llSize); //Consider making these template parameters
    EventManager(void); //NOTE: temporary  function
    ~EventManager(void);

    EventRef     queueEvent(Event* event);
    void         updateEvent(EventRef ref, Event* event);
    const Event* getNextEvent(void);
    void         deleteEvent(EventRef ref);
    void         clear(void);
//private Functions
private:
    void cbtUpdate(EventRef eRef);
    void cbtDelete(EventRef eRef);
    void cbtInsert(EventRef eRef);

    void insertInEventQ(EventRef eRef);
    void processOverflowList(void);
    void releaseEventRef(EventRef eRef);
//Members
private:
    //NEW: Even better, use the id_manager, and when no id is available push_back
    std::vector<Event*>    events_;

    //Priority queue vars
    std::vector<EventItem> eventItems_; //Same size as events_. In principle it is a linked-list/binary tree node.
    std::vector<EventRef>  llQueue_;    //The array of linked-lists. Holds the index of the first element each linked list.

    size_t currentIndex_;
    size_t baseIndex_;

    size_t llSize_;
    double scaleFactor_;

    //Complete binary tree vars
    std::vector<EventRef> nodes_;
    std::vector<size_t>   leafs_;

    int nCBTEvents_; //Current number of events in binary tree

    //EventRef manager vars
    int eRefCount_; //Current number of total events
    std::set<EventRef> available_;
};

#endif
