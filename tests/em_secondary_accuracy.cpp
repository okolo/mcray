#include <gtest/gtest.h>
#include "ICS.h"
#include "EMSecondarySampling.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <gsl/gsl_integration.h>

// Independent positive differential kernels and fixed-order quadrature from
// the validated correction branch. They do not call production primitives.
namespace secondary {
using Real = long double;
// Direct differential laws: arXiv:1106.5508v2, equations (11), (12).
// Overall factors cancel in a distribution conditional on s.
inline Real pairDensity(Real r, Real s) {
    const Real u=1-r, a=4/s;
    return (r*r/u+u+a/u-a*a/(4*r*u*u))/r/(1+2*(1-a)*a);
}
inline Real comptonDensity(Real y, Real s) {
    const Real a=1/s, b=1-a;
    return ((1+y*y)/2-2*a*(y-a)*(1-y)/(y*b*b))/y;
}
inline Real pairMinimum(Real s) {return 2/(s*(1+std::sqrt(1-4/s)));}
struct Quadrature {
    gsl_integration_glfixed_table* table;
    explicit Quadrature(int order):table(gsl_integration_glfixed_table_alloc(order)) {if(!table)throw std::runtime_error("Cannot allocate quadrature");}
    ~Quadrature(){gsl_integration_glfixed_table_free(table);}
    template<class F> Real logarithmic(F density,Real lo,Real hi) const {
        if(hi<=lo)return 0;
        if(!(lo>0))throw std::runtime_error("Nonpositive quadrature bound");
        const Real a=std::log(lo),b=std::log(hi);
        const int pieces=std::max(1,int(std::ceil((b-a)/.25L)));
        Real value=0;
        for(int j=0;j<pieces;++j)for(size_t i=0;i<table->n;++i) {
            double t,w;gsl_integration_glfixed_point(double(a+(b-a)*j/pieces),double(a+(b-a)*(j+1)/pieces),i,&t,&w,table);
            const Real x=std::exp(Real(t)),f=density(x);
            if(!std::isfinite(f)||f<0)throw std::runtime_error("Invalid independent differential density");
            value+=w*x*f;
        }
        return value;
    }
};

} // namespace secondary

using namespace Interactions;
using namespace mcray;
using secondary::Real;

TEST(SecondaryAccuracy, InvalidCDFIsFatalWithoutResampling) {
    try {
        detail::InvertSecondaryCDF([](Real){return 2.L;},0.L,1.L,.5L);
        FAIL()<<"Invalid CDF was accepted";
    } catch(Utils::Exception* error) {
        EXPECT_FALSE(error->Retryable());
        delete error;
    }
}
namespace {
template<class F> Real integrate(const secondary::Quadrature& q,F density,Real lo,Real hi) {
    if(hi<=lo)return 0;
    // The shared log quadrature passes bounds through double GSL interfaces.
    // Use a unit interval for very narrow ranges to preserve their width.
    if((hi-lo)/lo>=.1L)return q.logarithmic(density,lo,hi);
    Real sum=0;
    for(size_t i=0;i<q.table->n;++i) {
        double t,w;gsl_integration_glfixed_point(0,1,i,&t,&w,q.table);
        sum+=w*density(lo+(hi-lo)*t);
    }
    return (hi-lo)*sum;
}
}

TEST(SecondaryAccuracy, ComptonIntegralMatchesPositiveDensityAcrossSeriesBoundary) {
    secondary::Quadrature q(32);
    for(Real s:{1.00000001L,1.001L,1.01L,2.L,10.L,1e4L,1e8L,1e14L}) {
        const Real a=1/s;
        for(Real fraction:{1e-6L,.001L,.01L,.1L,.5L,.99L,1.L}) {
            const Real y=a+(1-a)*fraction;
            const Real direct=integrate(q,[&](Real x){return secondary::comptonDensity(x,s);},a,y);
            const Real analytic=detail::ComptonIntegral(a,y)*(1-a)/(2*a);
            SCOPED_TRACE(::testing::Message()<<"S="<<double(s)<<" fraction="<<double(fraction));
            ASSERT_GT(analytic,0);
            EXPECT_NEAR(double(analytic/direct),1,s<1.0001L?2e-8:2e-11);
        }
    }
}

TEST(SecondaryAccuracy, PairIntegralMatchesPositiveDensityAtBothEndpoints) {
    secondary::Quadrature q(32);
    for(Real s:{4.0001L,5.L,10.L,1e4L,1e8L,1e14L}) {
        const Real lo=secondary::pairMinimum(s);
        for(Real fraction:{1e-6L,.001L,.5L,.999999L,1.L}) {
            const Real r=lo+(.5L-lo)*fraction;
            const Real direct=integrate(q,[&](Real x){return secondary::pairDensity(x,s);},lo,r);
            const Real analytic=detail::PairIntegral(r,s);
            SCOPED_TRACE(::testing::Message()<<"S="<<double(s)<<" fraction="<<double(fraction));
            ASSERT_GT(analytic,0);
            EXPECT_NEAR(double(analytic/direct),1,2e-9);
        }
    }
}

TEST(SecondaryAccuracy, ProductionComptonQuantilesResolveRecoilAndPhotonTails) {
    secondary::Quadrature q(32);
    const double energy=1e20*units.eV, mass=Particle::Mass(Electron);
    for(double input:{1.0001,1.01,2.,10.,1e4,1e8,1e14}) {
        const double invariant=input*mass*mass;
        const Real s=Real(invariant)/(Real(mass)*mass),lo=1/s,hi=1-1e-8L;
        const Real norm=q.logarithmic([&](Real x){return secondary::comptonDensity(x,s);},lo,hi);
        ICS process(nullptr,energy*1e-8);
        Particle primary(Electron,0);primary.Energy=energy;
        for(unsigned seed=1;seed<=24;++seed) {
            Randomizer draw(seed),sample(seed);
            const Real uniform=draw.Rand();
            std::vector<Particle> products;
            process.SampleSecondaries(primary,products,invariant,sample);
            ASSERT_EQ(products.size(),2u);
            const Real photon=Real(products[0].Energy)/energy, electron=Real(products[1].Energy)/energy;
            const Real y=electron<.5L?electron:1-photon;
            const Real probability=q.logarithmic([&](Real x){return secondary::comptonDensity(x,s);},lo,y)/norm;
            SCOPED_TRACE(::testing::Message()<<"S="<<input<<" seed="<<seed);
            EXPECT_NEAR(double(probability),double(uniform),2e-10);
            EXPECT_NEAR(double(photon+electron),1,4e-16);
            EXPECT_GE(y,lo*(1-1e-14L));
            if(input>=1e8 && seed==1) {EXPECT_LT(y,1e-4L);}
        }
    }
}

TEST(SecondaryAccuracy, ProductionPairQuantilesAndEqualSplittingTail) {
    secondary::Quadrature q(32);
    const double energy=1e20*units.eV,mass=Particle::Mass(Electron);
    GammaPP process(nullptr);
    Particle primary(Photon,0);primary.Energy=energy;
    for(double input:{4.0001,5.,10.,1e4,1e8,1e14}) {
        const double invariant=input*mass*mass;
        const Real s=Real(invariant)/(Real(mass)*mass),lo=secondary::pairMinimum(s);
        const Real norm=q.logarithmic([&](Real x){return secondary::pairDensity(x,s);},lo,.5L);
        for(unsigned seed=1;seed<=24;++seed) {
            Randomizer draw(seed),sample(seed);draw.Rand();const Real uniform=draw.Rand();
            std::vector<Particle> products;process.SampleSecondaries(primary,products,invariant,sample);
            ASSERT_EQ(products.size(),2u);
            const Real r=Real(products[0].Energy)/energy;
            const Real probability=q.logarithmic([&](Real x){return secondary::pairDensity(x,s);},lo,r)/norm;
            SCOPED_TRACE(::testing::Message()<<"S="<<input<<" seed="<<seed);
            EXPECT_NEAR(double(probability),double(uniform),2e-10);
            EXPECT_GE(r,lo*(1-1e-14L));EXPECT_LE(r,.5L);
        }
    }
    Randomizer sample(1);double nearest=0;
    for(int i=0;i<5000;++i) {
        std::vector<Particle> products;process.SampleSecondaries(primary,products,10*mass*mass,sample);
        nearest=std::max(nearest,products[0].Energy/energy);
    }
    EXPECT_GT(nearest,.4999);
}

TEST(SecondaryAccuracy, InversionResolvesProbabilitiesNearBothLimits) {
    for(Real uniform:{1e-12L,1e-8L,.001L,.5L,.999L,1-1e-8L,1-1e-12L}) {
        for(Real s:{10.L,1e6L,1e14L}) {
            const Real lo=1/s,norm=detail::ComptonIntegral(lo,1-1e-8L);
            const Real value=detail::InvertSecondaryCDF([&](Real t){return detail::ComptonIntegral(lo,std::exp(t))/norm;},std::log(lo),std::log(1-1e-8L),uniform);
            const Real obtained=detail::ComptonIntegral(lo,std::exp(value))/norm;
            const double tolerance=double(1e-12L*std::min(uniform,1-uniform))+2e-15;
            EXPECT_NEAR(double(obtained),double(uniform),tolerance);
            const Real rmin=detail::PairMinimum(s),pnorm=detail::PairIntegral(.5L,s);
            const Real pair=detail::InvertSecondaryCDF([&](Real t){return detail::PairIntegral(std::exp(t),s)/pnorm;},std::log(rmin),std::log(.5L),uniform);
            EXPECT_NEAR(double(detail::PairIntegral(std::exp(pair),s)/pnorm),double(uniform),tolerance);
            if(uniform>1-1e-8L) {EXPECT_GT(std::exp(pair),.499999999L);}
        }
    }
}
