#include <gtest/gtest.h>
#include "BoundedLogEnvelope.h"
#include <gsl/gsl_errno.h>
#include "ICS.h"
#include <functional>
using namespace Utils;

namespace {
class Curved : public Function {
public:
    double f(double x) const override {return (1+.9*std::sin(7*std::log(x)))/x;}
    bool HasIntervalUpperBound() const override {return true;}
    double LogDensityUpperBound(double a,double b,double fa,double fb) const override {
        const long double h=std::log(static_cast<long double>(b)/a);
        return SamplingUpperRound(std::max(fa,fb)+.9L*49*h*h/8);
    }
};
class RareTail : public Function {
public:
    double f(double x) const override {return (x<2?1.:1e-12)/x;}
    bool HasIntervalUpperBound() const override {return true;}
    double LogDensityUpperBound(double a,double,double,double) const override {return SamplingUpperRound(a<2?1.:1e-12);}
};
class Bad : public Curved {
    double LogDensityUpperBound(double,double,double,double) const override {return 0;}
};
}

TEST(BoundedRejection, CurvatureEnvelopeCoversEntireOscillatorySupport) {
    Curved f;BoundedLogEnvelope envelope(f,1,std::exp(8.));
    ASSERT_TRUE(envelope.Supported());
    EXPECT_LE(envelope.Upper(),1.01*envelope.ObservedMaximum());
    for(int i=0;i<=100000;++i) {
        const auto p=envelope.Proposal(i/100000.);
        ASSERT_LE(p.y,envelope.Upper());
        EXPECT_TRUE(envelope.Accept(p,p.y/envelope.Upper()/2));
        if(p.y/envelope.Upper()<.99)EXPECT_FALSE(envelope.Accept(p,1.01*p.y/envelope.Upper()));
    }
}

TEST(BoundedRejection, RetainsTinyTailAndBothEndpoints) {
    RareTail f;BoundedLogEnvelope envelope(f,1,1000);
    ASSERT_TRUE(envelope.Supported());
    for(double u:{0.,1e-12,.1,.5,.99,1-1e-12,1.}) {
        const auto p=envelope.Proposal(u);
        EXPECT_GE(p.x,1);EXPECT_LE(p.x,1000);
        EXPECT_TRUE(envelope.Accept(p,0));
    }
    EXPECT_NEAR(envelope.Proposal(1).y,1e-12,1e-26);
    EXPECT_DOUBLE_EQ(envelope.Proposal(0).x,1);
    EXPECT_DOUBLE_EQ(envelope.Proposal(1).x,1000);
}

TEST(BoundedRejection, CollapsedAbsoluteLogCoordinatesCannotReturnAnEmptyInteraction) {
    gsl_set_error_handler_off();RareTail f;MathUtils math;mcray::Randomizer rng(12345);
    const double lo=std::exp(10.),hi=std::nextafter(lo,INFINITY);
    ASSERT_EQ(std::log(lo),std::log(hi));MathUtils::BoundedSamplingStats stats;
    try {
        double x,total;EXPECT_TRUE(math.SampleLogDistributionBounded(f,rng,x,total,lo,hi,.001,&stats));
        EXPECT_GT(total,0);EXPECT_GE(x,lo);EXPECT_LE(x,hi);
    } catch(Exception* e){EXPECT_FALSE(e->Retryable());delete e;}
    EXPECT_EQ(stats.fallback,1u);
}

TEST(BoundedRejection, FailsOnInvalidEnvelopeInsteadOfAdjustingIt) {
    Bad f;bool failed=false;
    try {BoundedLogEnvelope envelope(f,1,100);}catch(Exception* e){failed=true;delete e;}
    EXPECT_TRUE(failed);
}

TEST(BoundedRejection, DoesNotSilentlyIncreaseEnvelopeAfterConstruction) {
    struct Mutable : RareTail {mutable bool changed=false;double f(double x) const override{return changed?100/x:RareTail::f(x);}} f;
    BoundedLogEnvelope envelope(f,1,100);ASSERT_TRUE(envelope.Supported());f.changed=true;
    bool failed=false;
    try {envelope.Accept(envelope.Proposal(.5),0);}catch(Exception* e){failed=true;delete e;}
    EXPECT_TRUE(failed);
}

TEST(BoundedRejection, UnsupportedDensityUsesCDFBeforeDrawing) {
    struct Power : Function {double f(double x) const override{return x;}} f;
    MathUtils math;mcray::Randomizer a(12345),b(12345);double x,y,ix,iy;
    ASSERT_TRUE(math.SampleLogDistributionBounded(f,a,x,ix,1,10,.001));
    ASSERT_TRUE(math.SampleLogDistributionCDF(f,b,y,iy,1,10,.001));
    EXPECT_DOUBLE_EQ(x,y);EXPECT_DOUBLE_EQ(ix,iy);EXPECT_DOUBLE_EQ(a.Rand(),b.Rand());
}

TEST(BoundedRejection, QuadratureFailureCannotDiscardOnlyAnUnluckyHistory) {
    gsl_set_error_handler_off();Curved f;MathUtils math;mcray::Randomizer rng(12345);
    bool failed=false;
    try {double x,total;math.SampleLogDistributionBounded(f,rng,x,total,1,100,1e-30);}
    catch(Exception* e){failed=true;EXPECT_FALSE(e->Retryable());delete e;}
    EXPECT_TRUE(failed);
}

TEST(BoundedRejection, SamplesWithoutNarrowingOrPointMassFallback) {
    gsl_set_error_handler_off();Curved f;MathUtils math;mcray::Randomizer rng(12345);
    int below=0,above=0;
    for(int i=0;i<2000;++i){double x,total;ASSERT_TRUE(math.SampleLogDistributionBounded(f,rng,x,total,1,std::exp(8.),.001));below+=x<2;above+=x>1000;}
    EXPECT_GT(below,20);EXPECT_GT(above,20);
}

TEST(BoundedRejection, PhotonTableSlopeDropsIncludeInterpolationKnotsAndClones) {
    std::vector<double> x{1,2,4,8},y{1,.5,.015625,.00390625};
    GSLTableFunc table(x,y,new LogScale(),new LogScale());table.SetAutoLimits();
    EXPECT_GE(table.LogSlopeDrop(1.1,3.9),4);
    EXPECT_LT(table.LogSlopeDrop(1.1,1.9),1e-9);
    EXPECT_GE(table.LogSlopeDrop(.5,1.9),1);
    SafePtr<Function> copy(table.Clone());EXPECT_DOUBLE_EQ(copy->LogSlopeDrop(1.1,3.9),table.LogSlopeDrop(1.1,3.9));
    EXPECT_TRUE(std::isnan(table.LogSlopeDrop(1,10)));
}

TEST(BoundedRejection, AnalyticCrossSectionBoundsAndLogCurvatureCoverPhysicalIntervals) {
    using namespace mcray;using namespace Interactions;
    struct Inspect : BackgroundIntegral {
        std::function<void(const Function&)> inspect;
        BackgroundIntegral* Clone() const override {return new Inspect(*this);}
        double GetRateAndSampleS(const Function& f,const Particle&,Randomizer&,double& s,double) override {inspect(f);s=f.Xmin();return 1;}
        double GetRateS(const Function&,const Particle&,double) override {return 0;}
        double GetRateAndSampleK(const Function&,const Particle&,Randomizer&,double&,double) override {return 0;}
        double GetRateK(const Function&,const Particle&,double) override {return 0;}
    };
    auto* background=new Inspect;PBackgroundIntegral owner=background;
    background->inspect=[](const Function& f) {
        ASSERT_TRUE(f.HasIntervalUpperBound());
        for(int i=0;i<100;++i)for(double h:{.001,.03,.3,1.}) {
            const double a=f.Xmin()*std::exp(i*.18),b=a*std::exp(h),fa=f(a),fb=f(b);
            const double upper=f.IntervalUpperBound(a,b,fa,fb),k=f.LogConcavityBound(a,b,fa,fb);
            for(int j=0;j<=16;++j) {
                const double q=j/16.,s=a*std::exp(q*h),actual=f(s);
                ASSERT_LE(actual,upper*(1+1e-11));
                if(fa>0 && fb>0 && actual>0 && k>=0)
                    ASSERT_LE(std::log(actual),(1-q)*std::log(fa)+q*std::log(fb)+k*h*h*q*(1-q)/2+1e-9);
            }
        }
    };
    Randomizer rng(1234);double s;
    for(double energy:{1.01e12,1e13,1e14,1e15,1e18,1e20}) {
        ICS ics(background,1e12*units.eV);Particle electron(Electron,0);electron.Energy=energy*units.eV;
        ASSERT_TRUE(ics.SampleS(electron,s,rng));
    }
    GammaPP pp(background);Particle photon(Photon,0);photon.Energy=1e15*units.eV;
    ASSERT_TRUE(pp.SampleS(photon,s,rng));
}

TEST(BoundedRejection, ContinuousRedshiftMixturesFiniteBetaAndMonochromaticBackgrounds) {
    using namespace mcray;using namespace Interactions;
    cosmology.Init(.1,.692,67.8);const double t=2.7255*units.K;
    PlankBackground photons(t,t*.001,t*1000,0,.1);
    PBackgroundIntegral continuous=new ContinuousBackgroundIntegral(photons,.01,std::pow(10.,.01),.1,.001);
    struct Constant : Function {double value;explicit Constant(double v):value(v){}double f(double) const override{return value;}};
    // The monochromatic photons must permit the hard-ICS energy threshold.
    Constant energy(100*t),density(1);
    PBackgroundIntegral mono=new MonochromaticBackgroundIntegral(density,energy,1.1,.001);
    Randomizer rng(20261020);
    for(auto* background:{continuous.operator->(),mono.operator->()})for(double z:{0.,.003,.037,.099}) {
        ICS ics(background,1e12*units.eV);GammaPP pp(background);
        Particle electron(Electron,z),photon(Photon,z);electron.Energy=1e13*units.eV;photon.Energy=1e16*units.eV;
        for(int i=0;i<100;++i) {
            double a,b;ASSERT_TRUE(ics.SampleS(electron,a,rng));ASSERT_TRUE(pp.SampleS(photon,b,rng));
            EXPECT_GT(a,Particle::Mass(Electron)*Particle::Mass(Electron));
            EXPECT_GT(b,4*Particle::Mass(Electron)*Particle::Mass(Electron));
        }
    }
}
