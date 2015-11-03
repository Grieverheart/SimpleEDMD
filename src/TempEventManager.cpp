#include "TempEventManager.h"
#include <cassert>

TempEventManager::TempEventManager(void){
}

TempEventManager::~TempEventManager(void){
}

void TempEventManager::resize(size_t nPart){
    events_.resize(nPart);
}

void TempEventManager::init(void){
    //queue_.sort();
}

void TempEventManager::push(size_t pID, const ParticleEvent& event){
    queue_.push_back(event);
    EventRef event_ref = --queue_.end();
    assert(event.time_ == event_ref->time_);
    events_[pID].push_back(event_ref);
}

void TempEventManager::update(size_t pID){
    //queue_.sort();
}

void TempEventManager::clear(size_t pID){
    for(auto event_ref: events_[pID]){
        queue_.erase(event_ref);
    }
    events_[pID].clear();
}

ParticleEvent TempEventManager::getNextEvent(void){
    queue_.sort();
    auto event_ref = queue_.begin();
    auto event = *event_ref;
    for(auto event_itr = events_[event.pid_].begin(); event_itr != events_[event.pid_].end(); ++event_itr){
        if(*event_itr == event_ref){
            events_[event.pid_].erase(event_itr);
            break;
        }
    }
    queue_.pop_front();
    return event;
}
