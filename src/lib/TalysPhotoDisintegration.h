#ifndef MCRAY_TALYS_PHOTODISINTEGRATION_H
#define MCRAY_TALYS_PHOTODISINTEGRATION_H

#include "Background.h"
#include "Interaction.h"
#include "TableFunction.h"
#include "NuclearPhoton.h"
#include <array>
#include <map>
#include <memory>
#include <set>

namespace Interactions {

// Immutable raw cross sections and tabulated masses shared by process clones.
class TalysData {
public:
    enum class MissingPDPolicy { Strict, CRPropa };
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
                       bool crpropaFallback = false, bool loadPhotons = false,
                       MissingPDPolicy missingPD = MissingPDPolicy::Strict);
    const Isotope& Get(unsigned int a, unsigned int z) const;
    const Isotope* Find(unsigned int a, unsigned int z) const;
    bool PDInactive(const Isotope& isotope) const;
    std::string MissingPDMetadata(const std::set<Identity>& reached) const;
    bool CRPropaFallback() const { return fCRPropaFallback; }
    bool PhotonsLoaded() const { return fPhotonsLoaded; }
    bool GzkBundle() const { return fGzkBundle; }
    void SetParticle(mcray::Particle& p, unsigned int a, unsigned int z) const;
    std::set<Identity> ValidateNetwork(unsigned int a, unsigned int z) const;
    const std::string& Metadata() const { return fMetadata; }
    static const std::array<Identity, 6> Emitted;
    struct PhotonTransition {
        struct Line { double energy=0.; std::shared_ptr<const Utils::LinearFunc> sigma; };
        std::shared_ptr<const Utils::LinearFunc> total;
        std::vector<Line> lines;
    };
    void RequirePhotonCoverage(unsigned int a,unsigned int z) const;
    const PhotonTransition* Photons(Identity parent,Identity residual) const;

private:
    std::map<Identity, Isotope> fIsotopes;
    std::string fMetadata;
    MissingPDPolicy fMissingPD;
    bool fCRPropaFallback = false, fGzkBundle = false, fPhotonsLoaded = false;
    std::map<std::pair<Identity,Identity>,PhotonTransition> fPhotons;
};

class TalysPhotoDisintegration : public mcray::RandomInteraction {
public:
    struct Rates {
        double total = 0., exclusiveSum = 0.;
        std::vector<double> exclusive;
    };
    TalysPhotoDisintegration(mcray::BackgroundIntegral* background, std::shared_ptr<const TalysData> data, bool photons = false, bool strictPhotonYields = false);
    double Rate(const mcray::Particle& p) const override;
    std::vector<NuclearPhotonLine> PhotonYields(const mcray::Particle& p, TalysData::Identity residual) const;
    Rates ChannelRates(const mcray::Particle& p) const;
    bool GetSecondaries(mcray::Particle& p, std::vector<mcray::Particle>& products,
                        mcray::Randomizer& random) const override;
    // Conditional kernel only; validation callers supply the original branch law.
    bool GetChannelSecondaries(const mcray::Particle& p, size_t channel,
                               std::vector<mcray::Particle>& products, mcray::Randomizer& random) const;
    mcray::RandomInteraction* Clone() const override;
private:
    bool GenerateSecondaries(const mcray::Particle& p, const TalysData::Channel& channel,
                             std::vector<mcray::Particle>& products, mcray::Randomizer& random) const;
    const TalysData::Isotope& Validate(const mcray::Particle& p) const;
    Utils::SmartPtr<mcray::BackgroundIntegral> fBackground;
    std::shared_ptr<const TalysData> fData;
    bool fPhotons, fStrictPhotonYields;
};
}
#endif
