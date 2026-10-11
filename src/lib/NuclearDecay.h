#ifndef MCRAY_NUCLEAR_DECAY_H
#define MCRAY_NUCLEAR_DECAY_H

#include "TalysPhotoDisintegration.h"
#include "NuclearPhoton.h"

namespace Interactions {
// Partial mean lifetimes already include branching and ionization corrections.
// Absence of a channel is stable only with an explicit stable-state record.
class NuclearDecayData {
public:
    struct Channel {
        std::string code; // beta-, beta+, alpha, proton, neutron multiplicities
        double restRate = 0.;
        std::vector<NuclearPhotonLine> photons;
    };
    NuclearDecayData(const std::string& directory, std::shared_ptr<const TalysData> nuclei, bool loadPhotons = false);
    const std::vector<Channel>& Get(unsigned int a, unsigned int z) const;
    std::set<TalysData::Identity> ValidateNetwork(unsigned int a, unsigned int z,
                                                bool photodisintegration, bool gzk = false) const;
    void AppendResidual(const mcray::Particle& parent, unsigned int a, unsigned int z,
                        double energy, std::vector<mcray::Particle>& products) const;
    bool PhotonsLoaded() const { return fPhotonsLoaded; }
    const TalysData& Nuclei() const { return *fNuclei; }
    const std::string& Metadata() const { return fMetadata; }
private:
    std::shared_ptr<const TalysData> fNuclei;
    std::map<TalysData::Identity, std::vector<Channel>> fChannels;
    std::set<TalysData::Identity> fUnsupported;
    std::string fMetadata;
    bool fPhotonsLoaded = false;
};

class NuclearDecay : public mcray::RandomInteraction {
public:
    explicit NuclearDecay(std::shared_ptr<const NuclearDecayData> data, bool photons = false);
    double Rate(const mcray::Particle& p) const override;
    bool GetSecondaries(mcray::Particle& p, std::vector<mcray::Particle>& products,
                        mcray::Randomizer& random) const override;
    mcray::RandomInteraction* Clone() const override;
private:
    void BetaDecay(mcray::Particle& residual, bool plus,
                   std::vector<mcray::Particle>& products, mcray::Randomizer& random) const;
    std::shared_ptr<const NuclearDecayData> fData;
    bool fPhotons;
};
}
#endif
