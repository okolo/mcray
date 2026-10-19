#include <gtest/gtest.h>
#include "NumericalLogCDF.h"
#include "Randomizer.h"
#include <functional>

namespace {
struct Density : Utils::Function {
    std::function<double(double)> value;
    std::vector<double> knots;
    explicit Density(std::function<double(double)> f):value(std::move(f)){}
    double f(double x) const override {return value(x);}
    std::vector<double> IntegrationBreakpoints() const override {return knots;}
};
template<class F> void fatal(F action) {
    try {action();FAIL()<<"Invalid distribution was accepted";}
    catch (Utils::Exception* error) {EXPECT_FALSE(error->Retryable());delete error;}
}
}

TEST(LogCDF, AnalyticPowerLawQuantilesResolveBothProbabilityTails) {
    for (double alpha:{-3.,-.1,0.,.1,3.}) for (double scale:{1e-200,1.,1e200}) {
        const double width=alpha==0?80.:20.;
        Density f([=](double x){return scale*std::exp(alpha*std::log(x))/x;});
        const double xmin=std::exp(-width/2),xmax=std::exp(width/2);
        const long double lower=std::log(static_cast<long double>(xmin)),upper=std::log(static_cast<long double>(xmax));
        Utils::NumericalLogCDF cdf(f,xmin,xmax,1e-3);
        const long double norm=alpha==0?width:
            (std::exp(static_cast<long double>(alpha)*upper)-std::exp(static_cast<long double>(alpha)*lower))/alpha;
        EXPECT_NEAR(cdf.IntegralValue()/double(norm*scale),1,3e-7);
        for (double u:{1e-12,1e-9,1e-6,.01,.1,.5,.9,.99,1-1e-6,1-1e-9,1-1e-12}) {
            const double x=cdf.Quantile(u);
            const long double t=std::log(static_cast<long double>(x));
            long double mass;
            if (u<=.5) mass=alpha==0?t-lower:
                std::exp(static_cast<long double>(alpha)*lower)*std::expm1(alpha*(t-lower))/alpha;
            else mass=alpha==0?upper-t:
                -std::exp(static_cast<long double>(alpha)*upper)*std::expm1(alpha*(t-upper))/alpha;
            const long double probability=u<=.5?u:1-u;
            SCOPED_TRACE(::testing::Message()<<"alpha="<<alpha<<" scale="<<scale<<" u="<<u);
            // Returning a binary64 energy has an unavoidable probability step
            // at a nonzero-density endpoint. Include that representational
            // floor explicitly instead of claiming sub-ulp tail accuracy.
            const long double floor=4*std::numeric_limits<double>::epsilon()*
                std::exp(static_cast<long double>(alpha)*t)/norm/probability;
            EXPECT_NEAR(double(mass/norm/probability),1,3e-5+double(floor));
        }
    }
}

TEST(LogCDF, SmallSecondaryModeIsRetainedWithoutTailTruncation) {
    // Two separated log-normal modes; the left mode has only 1e-7 mass.
    const double sigma=.12,weight=1e-7;
    Density f([=](double x){const double t=std::log(x);return
        (weight*std::exp(-.5*std::pow((t+5)/sigma,2))+
         std::exp(-.5*std::pow((t-3)/sigma,2)))/x;});
    Utils::NumericalLogCDF cdf(f,std::exp(-10.),std::exp(10.),1e-3);
    const double norm=1+weight;
    for (double u:{1e-12,1e-9,4e-8,.01,.5,.99,1-1e-9,1-1e-12}) {
        const double t=std::log(cdf.Quantile(u));
        const double probability=u<=.5?
            (weight*.5*std::erfc(-(t+5)/(sigma*std::sqrt(2.)))+.5*std::erfc(-(t-3)/(sigma*std::sqrt(2.))))/norm:
            (weight*.5*std::erfc((t+5)/(sigma*std::sqrt(2.)))+.5*std::erfc((t-3)/(sigma*std::sqrt(2.))))/norm;
        EXPECT_NEAR(probability/(u<=.5?u:1-u),1,3e-5)<<u;
    }
}

TEST(LogCDF, ZeroDensityIntervalsAndZeroIntegral) {
    Density f([](double x){const double t=std::log(x);return t<0?0:1/x;});
    Utils::NumericalLogCDF cdf(f,std::exp(-4.),std::exp(4.),1e-3);
    EXPECT_NEAR(cdf.IntegralValue(),4,1e-10);
    for (double u:{1e-9,.01,.5,.99,1-1e-9}) EXPECT_NEAR(std::log(cdf.Quantile(u)),4*u,4e-5*std::min(u,1-u));
    Density empty([](double){return 0.;});
    Utils::NumericalLogCDF zero(empty,1,10,1e-3);
    EXPECT_EQ(zero.IntegralValue(),0);
    fatal([&]{zero.Quantile(.5);});
    Utils::MathUtils math;mcray::Randomizer random(1);double x=9,integral=1;
    EXPECT_FALSE(math.SampleLogDistributionCDF(empty,random,x,integral,1,10,1e-3));
    EXPECT_EQ(integral,0);EXPECT_EQ(x,0);
}

TEST(LogCDF, SuppliedFeatureKnotsResolveAPeakBetweenBroadQuadratureNodes) {
    const double center=.123456,width=1e-6;
    Density f([=](double x){return std::max(0.,1-std::abs(std::log(x)-center)/width)/x;});
    for (double t:{center-width,center,center+width}) f.knots.push_back(std::exp(t));
    Utils::NumericalLogCDF cdf(f,std::exp(-10.),std::exp(10.),1e-3);
    EXPECT_NEAR(cdf.IntegralValue()/width,1,1e-7);
    for (double u:{1e-4,.1,.5,.9,1-1e-4}) {
        const double t=std::log(cdf.Quantile(u));
        const double actual=u<=.5?.5*std::pow((t-center+width)/width,2):
            .5*std::pow((center+width-t)/width,2);
        EXPECT_NEAR(actual/(u<=.5?u:1-u),1,4e-5);
    }
}

TEST(LogCDF, InvalidDensitySupportToleranceAndProbabilityAreFatal) {
    Density negative([](double){return -1.;}),nan([](double){return NAN;}),ok([](double x){return 1/x;});
    fatal([&]{Utils::NumericalLogCDF cdf(negative,1,10,1e-3);});
    fatal([&]{Utils::NumericalLogCDF cdf(nan,1,10,1e-3);});
    fatal([&]{Utils::NumericalLogCDF cdf(ok,0,10,1e-3);});
    fatal([&]{Utils::NumericalLogCDF cdf(ok,1,10,1e-16);});
    Utils::NumericalLogCDF cdf(ok,1,10,1e-3);
    fatal([&]{cdf.Quantile(-.1);});fatal([&]{cdf.Quantile(NAN);});
    EXPECT_EQ(cdf.Quantile(0),1);EXPECT_EQ(cdf.Quantile(1),10);
}

TEST(LogCDF, FiniteSupportRatioMayExceedBinary64Range) {
    Density f([](double x){return 1/x;});
    f.knots={1e-200,1,1e200};
    Utils::NumericalLogCDF cdf(f,1e-250,1e250,1e-3);
    const long double lo=std::log(1e-250L),span=std::log(1e250L)-lo;
    EXPECT_NEAR(cdf.IntegralValue()/double(span),1,1e-12);
    for (double u:{1e-6,.01,.5,.99,1-1e-6}) {
        const double x=cdf.Quantile(u);
        EXPECT_NEAR(double((std::log(static_cast<long double>(x))-lo)/span),u,1e-10);
    }
}
