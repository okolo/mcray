//
// Created by ok on 10/05/18.
//

#ifndef CRBEAM_LOGGER_H
#define CRBEAM_LOGGER_H
#include <iostream>
#include "Particle.h"

namespace MCPropag {

    class ILogger {
    public:
        // should be thread-safe function
        virtual void log(const Particle &aParticle) = 0;
    };

    class Logger : public ILogger {
        std::ostream& fLog;
        bool fHeaderPrinted;
    public:
        Logger(std::ostream& aLog);
        // should be thread-safe function
        virtual void log(const Particle &aParticle);
        virtual void print_data(const Particle &p, std::ostream& aOut);
        virtual void print_header(std::ostream& aOut);
    };

}

#endif //CRBEAM_LOGGER_H
