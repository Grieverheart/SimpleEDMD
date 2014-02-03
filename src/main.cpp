#include <iostream>
#include "include/Time.h"
#include "include/EventManager.h"
#include "include/Simulation.h"

//Test the event Manager
void test1(void){
    EventManager manager(10.0, 50);

    std::vector<EventRef> refs;
    for(int i = 0; i < 100; ++i){
        EventRef ref = manager.queueEvent(new CollisionEvent(0.1 * double(i), i, i + 1));
        refs.push_back(ref);
    }

    for(int i = 0; i < 100; ++i){
        const Event* event = manager.getNextEvent();
        std::cout << event->time_ << std::endl;
    }
}

int main(int argc, char *argv[]){

    Simulation sim;

    sim.readConfig(argv[1]);

    test1();

    return 0;
}
