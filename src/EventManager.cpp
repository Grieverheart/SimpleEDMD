#include "include/EventManager.h"

EventRef EventManager::getEventRef(void){
    if(available_.empty()) return eRefCount_++;

    auto first       = available_.begin();
    EventRef ret_val = *first;

    available_.erase(first);

    return ret_val;
}

void EventManager::releaseEventRef(EventRef eRef){
    if(eRef < eRefCount_) available_.insert(eRef);
}

void EventManager::cbtUpdate(EventRef eRef){
    size_t father;
    for(father = leafs_[eRef] / 2; father > 0; father /= 2){
        if(nodes_[father] != eRef) break;

        EventRef leftNode  = nodes_[2 * father];
        EventRef rightNode = nodes_[2 * father + 1];
        if(events_[leftNode].time < events_[rightNode].rime){
            nodes_[father] = leftNode;
        }
        else nodes_[father] = rightNode;
    }

    for( ; father > 0; father /= 2){
        EventRef oldWinner = nodes_[father];

        EventRef leftNode  = nodes_[2 * father];
        EventRef rightNode = nodes_[2 * father + 1];
        if(events_[leftNode].time < events_[rightNode].time){
            nodes_[father] = leftNode;
        }
        else nodes_[father] = rightNode;

        if(nodes_[father] == oldWinner) return;
    }
}

void EventManager::cbtInsert(EventRef eRef){
    if(nCBTEvents_ == 0){
        nodes_[1] = eRef;
        ++nCBTEvents_;
        return;
    }

    eRef lastNode               = nodes_[nCBTEvents_];
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
    if(nodes_[lastLeaf - 1] != eRef){
        leafs_[nodes_[lastLeaf - 1]] = l / 2;
        nodes_[lastLeaf / 2] = nodes_[lastLeaf - 1];
        cbtUpdate(nodes_[lastLeaf - 1]);
    }
    else{
        leafs_[nodes_[lastLeaf]] = lastLeaf / 2;
        nodes_[lastLeaf / 2] = nodes_[lastLeaf];
        cbtUpdate(nodes_[lastLeaf]);
        --nCBTEvents_;
        return;
    }

    if(nodes_[lastLeaf] != eRef){
        nodes_[leafs_[eRef]] = nodes_[lastLeaf];
        leafs_[nodes_[lastLeaf]] = leafs_[eRef];
        cbtUpdate(nodes_[lastLeaf]);
    }

    --nCBTEvents_;
}

