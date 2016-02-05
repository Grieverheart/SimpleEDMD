#include "EventManager.h"
#include <cstdio>
#include <limits>
#include <cassert>

static const size_t NO_LINK = std::numeric_limits<size_t>::max();

struct EventItem{
    EventItem(void):
        previous_(NO_LINK), next_(NO_LINK), qIndex_(NO_LINK)
    {}

    EventRef previous_;
    EventRef next_;
    size_t   qIndex_; //Index in linked lists.
};

EventManager::EventManager(void):
    currentIndex_(0), baseIndex_(0),
    nCBTEvents_(0)
{}

EventManager::~EventManager(void){
}

//TODO: Perhaps we don't have to save everything, we could save just enough
//and use insertInEventQ() for all particles.
void serialize(Archive& ar, const EventManager& evt_mgr){
    serialize(ar, evt_mgr.events_);
    serialize(ar, evt_mgr.eventItems_);
    serialize(ar, evt_mgr.llQueue_);
    serialize(ar, evt_mgr.currentIndex_);
    serialize(ar, evt_mgr.baseIndex_);
    serialize(ar, evt_mgr.llSize_);
    serialize(ar, evt_mgr.scaleFactor_);
    serialize(ar, evt_mgr.nodes_);
    serialize(ar, evt_mgr.leafs_);
    serialize(ar, evt_mgr.nCBTEvents_);
}

void EventManager::resize(size_t nPart){
    events_.resize(nPart);
    eventItems_.resize(nPart);
    //Check these numbers
    leafs_.resize(nPart + 1);
    nodes_.resize(2 * nPart);
}

//TODO: Change to 'optimize'.
void EventManager::init(void){
    //Instrument queue
    {
        double tmax = 0.0;
        double tmin = std::numeric_limits<double>::max();
        int n_events = 0;
        for(size_t i = 0; i < events_.size(); ++i){
            if(events_[i].empty()) continue;

            double time = events_[i].top().time_;
            tmin = std::min(tmin, time);
            tmax = std::max(tmax, time);
            ++n_events;
        }
        double dtmax = tmax - tmin;

        scaleFactor_ = n_events / dtmax;
        llSize_      = scaleFactor_ * dtmax;
    }

    llQueue_.resize(llSize_ + 1, NO_LINK);

    for(size_t i = 0; i < events_.size(); ++i){
        if(!events_[i].empty()) insertInEventQ(i);
    }
}

//TODO: Need to properly clear
void EventManager::clear(void){
    for(size_t i = 0; i < events_.size(); ++i) clear(i);
}

void EventManager::clear(size_t pID){
    events_[pID].clear();
}

bool EventManager::empty(size_t pID)const{
    return events_[pID].empty();
}

void EventManager::push(size_t pID, const ParticleEvent& event){
    events_[pID].push(event);
}

void EventManager::update(size_t pID){
    insertInEventQ(pID);
}

void EventManager::insertInEventQ(size_t eRef){
    assert(!events_[eRef].empty());
    const ParticleEvent& event = events_[eRef].top();
    size_t index = (size_t)(scaleFactor_ * event.time_ - baseIndex_); 

    if(index > llSize_ - 1){
        //wraparound
        index -= llSize_;
        //If still too large, put into overflow list
        if(index + 1 >= currentIndex_) index = llSize_;
    }

    EventItem& eItem = eventItems_[eRef];

    if(eItem.qIndex_ != NO_LINK){
        deleteFromEventQ(eRef);
    }

    eItem.qIndex_ = index;

    if(index == currentIndex_) cbtInsert(eRef); //Insert in binary tree
    else{                                       //Insert in linked-list (push front)
        EventRef oldFirst = llQueue_[index];

        eItem.previous_ = NO_LINK;
        eItem.next_     = oldFirst;
        llQueue_[index] = eRef;

        if(oldFirst != NO_LINK) eventItems_[oldFirst].previous_ = eRef;
    }
}

void EventManager::processOverflowList(void){
    EventRef eRef   = llQueue_[llSize_];
    llQueue_[llSize_] = NO_LINK;

    while(eRef != NO_LINK){
        EventRef next = eventItems_[eRef].next_;
        insertInEventQ(eRef);
        eRef = next;
    }
}

void EventManager::deleteFromEventQ(size_t eRef){
    EventItem& eItem = eventItems_[eRef];
    size_t index = eItem.qIndex_;
    if(index == currentIndex_) cbtDelete(eRef); //Delete from binary tree
    else if(index != NO_LINK){                  //Delete from linked-list
        EventRef previous = eItem.previous_;
        EventRef next     = eItem.next_;
        if(previous == NO_LINK) llQueue_[index] = next; //If first node in linked-list
        else eventItems_[previous].next_ = next;

        if(next != NO_LINK) eventItems_[next].previous_ = previous;
    }
    eItem.qIndex_ = NO_LINK;
}

ParticleEvent EventManager::getNextEvent(void){
    while(nCBTEvents_ == 0){
        if(++currentIndex_ == llSize_){
            currentIndex_  = 0;
            baseIndex_    += llSize_;
            processOverflowList();
        }

        //populate binary tree
        for(EventRef eRef = llQueue_[currentIndex_]; eRef != NO_LINK; eRef = eventItems_[eRef].next_){
            cbtInsert(eRef);
        }

        llQueue_[currentIndex_] = NO_LINK;
    }

    EventRef eRef = nodes_[1]; //The root contains the earliest event

    return events_[eRef].pop();
}

void EventManager::cbtUpdate(EventRef eRef){
    size_t father = leafs_[eRef] >> 1;
    for(; (nodes_[father] == eRef) && (father > 0); father >>= 1){
        EventRef leftNode  = nodes_[2 * father];
        EventRef rightNode = nodes_[2 * father + 1];

        nodes_[father] = (events_[leftNode].top().time_ < events_[rightNode].top().time_)? leftNode : rightNode;
    }

    for(; father > 0; father >>= 1){
        EventRef oldWinner = nodes_[father];

        EventRef leftNode  = nodes_[2 * father];
        EventRef rightNode = nodes_[2 * father + 1];

        nodes_[father] = (events_[leftNode].top().time_ < events_[rightNode].top().time_)? leftNode : rightNode;

        if(nodes_[father] == oldWinner) return;
    }
}

void EventManager::cbtInsert(EventRef eRef){
    if(nCBTEvents_ == 0){
        nodes_[1] = eRef;
        leafs_[eRef] = 1;
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
        leafs_[nodes_[lastLeaf]] = lastLeaf >> 1;
        nodes_[lastLeaf >> 1] = nodes_[lastLeaf];
        cbtUpdate(nodes_[lastLeaf]);
        --nCBTEvents_;
        return;
    }

    leafs_[nodes_[lastLeaf - 1]] = lastLeaf >> 1;
    nodes_[lastLeaf >> 1] = nodes_[lastLeaf - 1];
    cbtUpdate(nodes_[lastLeaf - 1]);

    if(nodes_[lastLeaf] != eRef){
        nodes_[leafs_[eRef]] = nodes_[lastLeaf];
        leafs_[nodes_[lastLeaf]] = leafs_[eRef];
        cbtUpdate(nodes_[lastLeaf]);
    }

    --nCBTEvents_;
}

