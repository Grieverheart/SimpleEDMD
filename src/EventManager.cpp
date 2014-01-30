#include "include/EventManager.h"
#include <cstdio>
#define EMPTY -1

struct EventItem{
    EventRef previous_;
    EventRef next_;
    size_t   qIndex_;
};

EventManager::EventManager(double scaleFactor, int llSize):
    currentIndex_(0), baseIndex_(0),
    llSize_(llSize), scaleFactor_(scaleFactor),
    nCBTEvents_(0),
    eRefCount_(0)
{
    llQueue_.resize(llSize_ + 1, EMPTY);
    ////Fiddle with numbers
    //events_.reserve(1000);
    //eventItems_.reserve(1000);
    //leafs_.reserve(1000);
    //nodes_.reserve(1001);
}

EventManager::~EventManager(void){
    for(auto event: events_){
        delete event;
    }
    currentIndex_ = 0;
    baseIndex_    = 0;
    nCBTEvents_   = 0;
    eRefCount_    = 0;

    available_.clear();
}

void EventManager::clear(void){
    for(auto event: events_){
        delete event;
    }
}

EventRef EventManager::queueEvent(Event* event){
    EventRef  eRef;
    EventItem eItem;
    if(available_.empty()){
        eRef = eRefCount_++;

        events_.push_back(event);
        eventItems_.push_back(eItem);
        //Resize binary tree like this
        size_t nEvents = events_.size();
        leafs_.resize(nEvents + 1);
        nodes_.resize(2 * nEvents);
    }
    else{
        auto first = available_.begin();
        eRef       = *first;
        available_.erase(first);
        //We need to delete the previous event now, because when returning
        //the next event even if we were to return a copy and destroy the local
        //one, the Event's destructor would delete the data member
        delete events_[eRef];
        events_[eRef] = event;
    }

    insertInEventQ(eRef);

    return eRef;
}

void EventManager::insertInEventQ(EventRef eRef){
    size_t index = (size_t)(scaleFactor_ * events_[eRef]->time_ - baseIndex_); 

    if(index > llSize_ - 1){
        index -= llSize_;
        if(index + 1 >= currentIndex_) index = llSize_;
    }

    EventItem& eItem = eventItems_[eRef];

    eItem.qIndex_ = index;

    if(index == currentIndex_) cbtInsert(eRef); //Insert in binary tree
    else{                                       //Insert in linked-list
        EventRef oldFirst = llQueue_[index];

        eItem.previous_ = EMPTY;
        eItem.next_     = oldFirst;
        llQueue_[index] = eRef;

        if(oldFirst != EMPTY){
            eventItems_[oldFirst].previous_ = eRef;
        }
    }
}

void EventManager::processOverflowList(void){
    size_t index    = llSize_;
    EventRef eRef   = llQueue_[index];
    llQueue_[index] = EMPTY;

    while(eRef != EMPTY){
        EventRef next = eventItems_[eRef].next_;
        insertInEventQ(eRef);
        eRef = next;
    }
}

void EventManager::deleteEvent(EventRef eRef){
    //For safety, we might have to first check if the event occurs
    //if(available_.find(eRef) != available_.end()) return; //Event does not exist
    EventItem& eItem = eventItems_[eRef];
    size_t index     = eItem.qIndex_;
    if(index == currentIndex_) cbtDelete(eRef);
    else{
        EventRef previous = eItem.previous_;
        EventRef next     = eItem.next_;
        if(previous == EMPTY) llQueue_[index] = next;
        else eventItems_[previous].next_ = next;

        if(next != EMPTY) eventItems_[next].previous_ = previous;
    }
    releaseEventRef(eRef);
}

//CAUTION: Experimental
void EventManager::updateEvent(EventRef eRef, Event* event){
    EventItem& eItem = eventItems_[eRef];
    size_t index     = eItem.qIndex_;

    delete events_[eRef];
    events_[eRef] = event;

    if(index == currentIndex_) cbtUpdate(eRef);
    else{
        //Delete old event from linked list
        EventRef previous = eItem.previous_;
        EventRef next     = eItem.next_;
        if(previous == EMPTY) llQueue_[index] = next;
        else eventItems_[previous].next_ = next;

        if(next != EMPTY) eventItems_[next].previous_ = previous;

        //Insert new event in linked list
        size_t newIndex = (size_t)(scaleFactor_ * events_[eRef]->time_ - baseIndex_); 
        EventRef oldFirst = llQueue_[newIndex];

        eItem.previous_    = EMPTY;
        eItem.next_        = oldFirst;
        eItem.qIndex_      = newIndex;
        llQueue_[newIndex] = eRef;

        if(oldFirst != EMPTY){
            eventItems_[oldFirst].previous_ = eRef;
        }
    }
    releaseEventRef(eRef);
}

//We might have to change this function at some point to return an EventRef instead
//and add a function const Event& getEvent(EventRef ref); deleting multiple events
//associated with one particle can then skip the event that just occured
const Event* EventManager::getNextEvent(void){
    while(nCBTEvents_ == 0){
        ++currentIndex_;
        if(currentIndex_ == llSize_){
            currentIndex_  = 0;
            baseIndex_    += llSize_;
            processOverflowList();
        }

        //populate binary tree
        for(EventRef eRef = llQueue_[currentIndex_]; eRef != EMPTY; eRef = eventItems_[eRef].next_)
            cbtInsert(eRef);

        llQueue_[currentIndex_] = EMPTY;
    }

    EventRef eRef = nodes_[1]; //The root contains the earliest event
    cbtDelete(eRef);
    releaseEventRef(eRef);

    return events_[eRef];
}

void EventManager::releaseEventRef(EventRef eRef){
    if(eRef < eRefCount_) available_.insert(eRef);
}

void EventManager::cbtUpdate(EventRef eRef){
    size_t father;
    for(father = leafs_[eRef] / 2; (nodes_[father] == eRef) && (father > 0); father /= 2){
        EventRef leftNode  = nodes_[2 * father];
        EventRef rightNode = nodes_[2 * father + 1];

        nodes_[father] = (events_[leftNode]->time_ < events_[rightNode]->time_)? leftNode : rightNode;
    }

    for( ; father > 0; father /= 2){
        EventRef oldWinner = nodes_[father];

        EventRef leftNode  = nodes_[2 * father];
        EventRef rightNode = nodes_[2 * father + 1];

        nodes_[father] = (events_[leftNode]->time_ < events_[rightNode]->time_)? leftNode : rightNode;

        if(nodes_[father] == oldWinner) return;
    }
}

void EventManager::cbtInsert(EventRef eRef){
    if(nCBTEvents_ == 0){
        nodes_[1] = eRef;
        ++nCBTEvents_;
        return;
    }

    EventRef lastNode           = nodes_[nCBTEvents_];
    nodes_[2 * nCBTEvents_]     = lastNode;
    nodes_[2 * nCBTEvents_ + 1] = eRef;

    leafs_[lastNode] = 2 * nCBTEvents_;
    leafs_[eRef]     = 2 * nCBTEvents_ + 1;

    ++nCBTEvents_;
    cbtUpdate(lastNode);
}

void EventManager::cbtDelete(EventRef eRef){
    if(nCBTEvents_ < 2){
        nodes_[1] = 0;
        leafs_[0] = 1;
        --nCBTEvents_;
        return;
    }

    size_t lastLeaf = 2 * nCBTEvents_ - 1;
    if(nodes_[lastLeaf - 1] == eRef){
        leafs_[nodes_[lastLeaf]] = lastLeaf / 2;
        nodes_[lastLeaf / 2] = nodes_[lastLeaf];
        cbtUpdate(nodes_[lastLeaf]);
        --nCBTEvents_;
        return;
    }

    leafs_[nodes_[lastLeaf - 1]] = lastLeaf / 2;
    nodes_[lastLeaf / 2] = nodes_[lastLeaf - 1];
    cbtUpdate(nodes_[lastLeaf - 1]);

    if(nodes_[lastLeaf] != eRef){
        nodes_[leafs_[eRef]] = nodes_[lastLeaf];
        leafs_[nodes_[lastLeaf]] = leafs_[eRef];
        cbtUpdate(nodes_[lastLeaf]);
    }

    --nCBTEvents_;
}

