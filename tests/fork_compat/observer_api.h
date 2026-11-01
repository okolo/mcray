#pragma once

struct ObserverInput {
    int mode; // 1: strict inside, 2: strict outside, -1: no observer.
    double center[3], radius;
    long double position[3], direction[3], initialTime, duration;
};
struct ObserverPoint {
    bool reached;
    long double step;
};
struct ObserverFlight {
    bool reached, pastEnd;
    long double position[3], elapsed, endError, endErrorUlps;
    unsigned long steps;
};
ObserverPoint HistoricalObserverPoint(const ObserverInput&);
ObserverFlight HistoricalObserverFlight(const ObserverInput&);
