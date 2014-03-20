#ifndef __PERIODIC_CONDITION_H
#define __PERIODIC_CONDITION_H

#include <functional>

class PeriodicCondition{
public:
    PeriodicCondition(double time):
        type_(TIME), time_(time)
    {}

    PeriodicCondition(unsigned int nEvents):
        type_(EVENT), nEvents_(nEvents)
    {}

    ~PeriodicCondition(void){}

    bool operator()(double time, unsigned int nEvents)const{
        return (type_ == TIME)? (time >= time_): (nEvents >= nEvents_);
    }

    typedef double (*time_function)(double);
    typedef unsigned int (*event_function)(unsigned int);
    typedef std::function<double(double)> time_functor;
    typedef std::function<unsigned int(unsigned int)> event_functor;

    void setNextFunction(const time_function& func){
        timeNext_ = func;
    }
    void setNextFunction(const event_function& func){
        eventNext_ = func;
    }

    void next(void){
        switch(type_){
        case TIME:
            time_ = timeNext_(time_);
            break;
        case EVENT:
            nEvents_ = eventNext_(nEvents_);
            break;
        default:
            break;
        }
    };

private:
    enum type{TIME, EVENT};
    type type_;
    union{
        double       time_;
        unsigned int nEvents_;
    };
    union{
        time_functor  timeNext_;
        event_functor eventNext_;
    };
};

#endif
