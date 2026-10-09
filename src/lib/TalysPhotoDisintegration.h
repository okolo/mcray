#ifndef MCRAY_TALYS_PHOTODISINTEGRATION_H
#define MCRAY_TALYS_PHOTODISINTEGRATION_H

#include "Background.h"
#include "Interaction.h"
#include "TableFunction.h"
#include <array>
#include <map>
#include <memory>
#include <set>

namespace Interactions {

// Immutable raw cross sections and tabulated masses shared by process clones.
class TalysData {
public:
    using Identity = std::pair<unsigned int, unsigned int>; // A,Z
    struct Channel {
        std::string code;
        std::array<unsigned int, 6> multiplicity{}; // n,p,d,t,He3,He4
        Identity residual;
        std::shared_ptr<const Utils::LinearFunc> sigma;
    };
    struct Isotope {
        double mass = 0.;
        bool terminal = false;
        bool approximateMass = false;
        std::shared_ptr<const Utils::LinearFunc> total;
        std::vector<Channel> channels;
    };
    explicit TalysData(const std::string& directory, bool loadCrossSections = true,
                       bool crpropaFallback = false);
    const Isotope& Get(unsigned int a, unsigned int z) const;
    const Isotope* Find(unsigned int a, unsigned int z) const;
    bool CRPropaFallback() const { return fCRPropaFallback; }
    bool GzkBundle() const { return fGzkBundle; }
    void SetParticle(mcray::Particle& p, unsigned int a, unsigned int z) const;
    std::set<Identity> ValidateNetwork(unsigned int a, unsigned int z) const;
    const std::string& Metadata() const { return fMetadata; }
    static const std::array<Identity, 6> Emitted;
private:
    std::map<Identity, Isotope> fIsotopes;
    std::string fMetadata;
    bool fCRPropaFallback = false, fGzkBundle = false;
};

class TalysPhotoDisintegration : public mcray::RandomInteraction {
public:
    struct Rates {
        double total = 0., exclusiveSum = 0.;
        std::vector<double> exclusive;
    };
    TalysPhotoDisintegration(mcray::BackgroundIntegral* background, std::shared_ptr<const TalysData> data);
    double Rate(const mcray::Particle& p) const override;
    Rates ChannelRates(const mcray::Particle& p) const;
    bool GetSecondaries(mcray::Particle& p, std::vector<mcray::Particle>& products,
                        mcray::Randomizer& random) const override;
    mcray::RandomInteraction* Clone() const override;
private:
    const TalysData::Isotope& Validate(const mcray::Particle& p) const;
    Utils::SmartPtr<mcray::BackgroundIntegral> fBackground;
    std::shared_ptr<const TalysData> fData;
};
}
#endif
