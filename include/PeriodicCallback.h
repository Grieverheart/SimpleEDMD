#ifndef __PERIODIC_CALLBACK_H
#define __PERIODIC_CALLBACK_H

#include <functional>

class PeriodicCallback{
public:
    PeriodicCallback(double t0):
        time_(t0),
        timeNext_([](double time){return time + 10000.0;}),
        caller_([](double) -> bool {return true;})
    {}

    ~PeriodicCallback(void){}

    void operator()(double time){
        if(time >= time_ && caller_(time_)) time_ = timeNext_(time_);
    }

    typedef std::function<double(double)> time_functor;
    typedef std::function<bool(double)> call_functor;

    void setNextFunction(const time_functor& func){
        timeNext_ = func;
    }

    void setCallback(const call_functor& func){
        caller_ = func;
    }

    double getTime(void)const{
        return time_;
    }

    void reset(void){
        time_ = 0.11;
    }

private:
    double       time_;
    time_functor timeNext_;
    call_functor caller_;
};

#endif
