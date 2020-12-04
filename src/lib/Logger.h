/*
 * Logger.h
 *
 * Author:
 *       Oleg Kalashev
 *
 * Copyright (c) 2020 Institute for Nuclear Research, RAS
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


#ifndef CRBEAM_LOGGER_H
#define CRBEAM_LOGGER_H
#include <iostream>
#include "Particle.h"

namespace mcray {

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
