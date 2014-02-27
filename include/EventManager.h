#ifndef __EVENT_MANAGER_H
#define __EVENT_MANAGER_H

#include "Event.h"
#include "BinaryHeap.h"
#include <vector>
#include <cstddef>

typedef int EventRef;

struct EventItem;

//Consider collecting statistics to improve estimates of scaleFactor and llSize after clears
class EventManager{
public:
     EventManager(void); //Consider making these template parameters
    ~EventManager(void);

    void    resize(size_t nPart);
    void    init(void);
    void    push(size_t pID, Event* event);
    void    update(size_t pID);
    void    clear(void);
    void    clear(size_t pID);
    bool    empty(size_t pID)const;
    Event*  getNextEvent(void);
//private Functions
private:
    void cbtUpdate(EventRef eRef);
    void cbtDelete(EventRef eRef);
    void cbtInsert(EventRef eRef);

    void insertInEventQ(size_t eRef);
    void deleteFromEventQ(size_t eRef);
    void processOverflowList(void);
//Members
private:
    typedef BinaryHeap<Event*, EventPtrLess> PEL;
    std::vector<PEL> events_;
    //We use this as a temporary measure, since with cell list, all particles will have at least one event
    std::vector<bool> isInserted_; 

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
};

#endif
