#ifndef MCRAY_NUCLEAR_PPP_H
#define MCRAY_NUCLEAR_PPP_H

#include "PPP.h"
#include "ProtonPP.h"
#include <cmath>

namespace Interactions {

// Native proton kernel evaluated at equal Lorentz factor. The surrogate is
// a copy so secondary metadata retain the nuclear parent's history.
inline mcray::Particle NuclearPairProton(const mcray::Particle& nucleus) {
    nucleus.ValidateIdentity();
    if (nucleus.Type != mcray::Nucleus || nucleus.MassNumber() > 56 ||
        nucleus.MassModel() != mcray::Particle::NuclearMassModel::TALYS)
        Utils::Exception::Throw("Nuclear PPP requires A<=56 and a tabulated nuclear mass");
    if (!std::isfinite(nucleus.Energy) || nucleus.Energy < nucleus.Mass())
        Utils::Exception::Throw("Nuclear PPP requires finite energy >= mass");
    mcray::Particle proton = nucleus;
    proton.SetType(mcray::Proton);
    proton.Energy = nucleus.Energy * (proton.Mass() / nucleus.Mass());
    return proton;
}

class NuclearPPPcel : public mcray::CELInteraction {
public:
    explicit NuclearPPPcel(mcray::BackgroundIntegral* background)
        : fBackground(background), fProton(background) {}
    double Rate(const mcray::Particle& p) const override {
        if (p.Type != mcray::Nucleus) return 0.;
        const auto proton = NuclearPairProton(p);
        const double z = p.AtomicNumber();
        return z*z * (proton.Mass()/p.Mass()) * fProton.Rate(proton);
    }
    mcray::CELInteraction* Clone() const override {
        return new NuclearPPPcel(fBackground->Clone());
    }
private:
    Utils::SmartPtr<mcray::BackgroundIntegral> fBackground;
    ProtonPPcel fProton;
};

inline double NuclearPairScale(double scale) {
    if (!std::isfinite(scale) || scale <= 0.)
        Utils::Exception::Throw("Nuclear PPP thinning must be finite and positive");
    return scale;
}

class NuclearPPP : public mcray::RandomInteractionS {
public:
    NuclearPPP(mcray::BackgroundIntegral* background, double scale = 1.)
        : fBackground(background), fProton(background, NuclearPairScale(scale)), fScale(scale) {}
    double Rate(const mcray::Particle& p) const override {
        if (p.Type != mcray::Nucleus) return 0.;
        const auto proton = NuclearPairProton(p);
        const double z = p.AtomicNumber();
        return z*z * fProton.Rate(proton);
    }
    bool KillsPrimary() const override { return false; }
    bool SampleS(const mcray::Particle& p, double& s, mcray::Randomizer& random) const override {
        return fProton.SampleS(NuclearPairProton(p), s, random);
    }
    void SampleSecondaries(mcray::Particle& p, std::vector<mcray::Particle>& out,
                           double s, mcray::Randomizer& random) const override {
        auto proton = NuclearPairProton(p);
        fProton.SampleSecondaries(proton, out, s, random);
    }
    mcray::RandomInteraction* Clone() const override {
        return new NuclearPPP(fBackground->Clone(), fScale);
    }
private:
    Utils::SmartPtr<mcray::BackgroundIntegral> fBackground;
    PPP fProton;
    double fScale;
};
}
#endif
