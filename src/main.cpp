#include <iostream>
#include "include/Box.h"
#include "include/Time.h"
#include "include/Sphere.h"
#include "include/EventManager.h"


Time timeOfCollision(Sphere &sphA, Sphere &sphB){
    return 0.0;
}

Time timeOfCollision(Sphere &sphere, Box &box){
    return 0.0;
}

int main(int argc, char *argv[]){

    EventManager manager(10.0, 50);

    std::vector<EventRef> refs;
    for(int i = 0; i < 100; ++i){
        EventRef ref = manager.queueEvent(new Event(0, 0.1 * double(i)));
        refs.push_back(ref);
    }

    for(int i = 0; i < 100; ++i){
        const Event* event = manager.getNextEvent();
        std::cout << event->time_ << std::endl;
    }

    return 0;
}
