/*
 * Logger.cpp
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


#include "Logger.h"
#include <iomanip>

namespace mcray {
    Logger::Logger(std::ostream& aOut) :
            fLog(aOut),
            fHeaderPrinted(false)
    {
    }

    void Logger::log(const Particle &p) {
#pragma omp critical (ParticleStack)
        {
            if(!fHeaderPrinted){
                fLog << "# id\tparent_id\t";
                print_header(fLog);
                fLog << std::endl;
                fHeaderPrinted = true;
            }
            fLog << p.id << "\t" << p.fPrevId << "\t";
            print_data(p, fLog);
            fLog << std::endl;
        }
    }

    void Logger::print_data(const Particle &p, std::ostream& aOut){
        double t = p.Time.t()-p.SourceParticle->Time.t();
        double t_1 = 1.;
        double zDif = 0.;
        if(t>0){
            t_1 = 1./t;
            zDif = 1.-t_1*p.X[2];
        }
        aOut << p.Type << "\t";
        std::streamsize orig_pr = aOut.precision(12);
        aOut << t/units.Mpc << std::setprecision(orig_pr) << "\t" << p.ElectricCharge()
             << "\t" << p.Energy/units.eV << "\t" << p.Pdir[0] << "\t" << p.Pdir[1] << "\t" << t_1*p.X[0]
             << "\t" << t_1*p.X[1]  << "\t" << zDif;
    }

    void Logger::print_header(std::ostream& aOut){
        aOut << "type\tt/Mpc\tcharge\tE/eV\tp[0]/p\tp[1]/p\tx/t\ty/t\t1-z/t";
    }
}