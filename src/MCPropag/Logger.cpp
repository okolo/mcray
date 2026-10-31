//
// Created by ok on 10/05/18.
//
//#include <math.h>
#include "Logger.h"
#include <iomanip>

namespace MCPropag {
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