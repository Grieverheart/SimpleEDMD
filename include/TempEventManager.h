#ifndef __TEMP_EVENT_MANAGER_H
#define __TEMP_EVENT_MANAGER_H

#include "Event.h"
#include <vector>
#include <list>
#include <cstddef>

//Consider collecting statistics to improve estimates of scaleFactor and llSize after clears
class TempEventManager{
public:
     TempEventManager(void); //Consider making these template parameters
    ~TempEventManager(void);

    void resize(size_t nPart);
    void init(void);
    void push(size_t pID, const ParticleEvent& event);
    void update(size_t pID);
    void clear(size_t pID);

    ParticleEvent getNextEvent(void);
//Members
private:
    using EventList = std::list<ParticleEvent>;
    using EventRef = EventList::iterator;
    EventList queue_;
    std::vector<std::vector<EventRef>> events_;
};

#endif
