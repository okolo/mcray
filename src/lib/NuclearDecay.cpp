#include "NuclearDecay.h"
#include <algorithm>
#include <fstream>
#include <sstream>

namespace Interactions {
using namespace mcray;
using Utils::Exception;
namespace {
std::vector<std::string> Lines(const std::string& file) {
    std::ifstream input(file);
    if (!input) Exception::Throw("Cannot read nuclear decay data: " + file);
    std::vector<std::string> result;
    for (std::string line; std::getline(input,line);) {
        const auto start = line.find_first_not_of(" \t\r");
        if (start != std::string::npos && line[start] != '#') result.push_back(line);
    }
    return result;
}
TalysData::Identity Identity(std::istringstream& row) {
    unsigned int a=0,z=0;
    if (!(row>>a>>z) || a<2 || a>56 || !z || z>a)
        Exception::Throw("Invalid nuclear decay identity");
    return {a,z};
}
void End(std::istringstream& row) {
    std::string extra;
    if (row>>extra) Exception::Throw("Unexpected nuclear decay data column");
}
// Enumerate every intermediate identity: compound channels follow CRPropa order.
std::vector<TalysData::Identity> Products(TalysData::Identity parent, const std::string& code) {
    int a=parent.first,z=parent.second;
    std::vector<TalysData::Identity> products;
    auto add=[&] {
        if (a<1 || z<0 || z>a)
            Exception::Throw("Invalid residual in nuclear decay channel");
        products.emplace_back(a,z);
    };
    for (int i=0;i<code[0]-'0';++i) { ++z; add(); }
    for (int i=0;i<code[1]-'0';++i) { --z; add(); }
    const std::array<TalysData::Identity,3> emitted={{{4,2},{1,1},{1,0}}};
    for (size_t i=0;i<emitted.size();++i)
        for (int count=0;count<code[i+2]-'0';++count) {
            a-=emitted[i].first; z-=emitted[i].second; add();
            products.push_back(emitted[i]);
        }
    return products;
}
}

NuclearDecayData::NuclearDecayData(const std::string& directory, std::shared_ptr<const TalysData> nuclei)
    : fNuclei(std::move(nuclei)) {
    if (!fNuclei) Exception::Throw("Nuclear decays require tabulated nuclear masses");
    const auto provenance=Lines(directory+"/decay_provenance.txt");
    if (provenance.size()<2 || provenance[0]!="mcray-decay-v1" ||
        provenance[1]!="34b10e27375c0e9894d7c984ac15fb96e12975ef")
        Exception::Throw("Unsupported nuclear decay format or source revision");
    for (const auto& line:provenance) fMetadata+="# NuclearDecay "+line+'\n';
    for (const auto& line:Lines(directory+"/decay_manifest.json"))
        fMetadata+="# NuclearDecay manifest "+line+'\n';
    if (fNuclei->GzkBundle()) for (const auto& line:Lines(directory+"/unsupported_decays.tsv")) {
        std::istringstream row(line); const auto key=Identity(row); std::string mode;
        if (!(row>>mode)) Exception::Throw("Invalid unsupported nuclear decay record");
        End(row); fUnsupported.insert(key);
    }
    for (const auto& line:Lines(directory+"/decay_stable.tsv")) {
        std::istringstream row(line); auto key=Identity(row); End(row);
        fNuclei->Get(key.first,key.second);
        if (!fChannels.emplace(key,std::vector<Channel>{}).second)
            Exception::Throw("Duplicate nuclear decay stable state");
    }
    auto decayLines=Lines(directory+"/decays.tsv");
    if (fNuclei->GzkBundle() && fNuclei->CRPropaFallback()) {
        const auto fallback=Lines(directory+"/fallback_decays.tsv");
        for (const auto& line:fallback) {
            std::istringstream row(line); auto key=Identity(row); std::string code; double tau=0.;
            if (!(row>>code>>tau) || tau!=1e-99 ||
                code!=(key.second>key.first-key.second ? "00010" : "00001"))
                Exception::Throw("Invalid CRPropa fallback decay");
            End(row);
        }
        decayLines.insert(decayLines.end(),fallback.begin(),fallback.end());
        fMetadata+="# NuclearDecay artificial-dripping=CRPropa; tau_s=1e-99\n";
    }
    for (const auto& line:decayLines) {
        std::istringstream row(line); auto key=Identity(row); Channel channel; double tau=0.;
        if (!(row>>channel.code>>tau) || channel.code.size()!=5 ||
            channel.code.find_first_not_of("0123456789")!=std::string::npos ||
            channel.code=="00000" || !std::isfinite(tau) || tau<=0.)
            Exception::Throw("Invalid nuclear decay channel or partial lifetime");
        End(row);
        if (tau<1e-90 && !fNuclei->CRPropaFallback())
            Exception::Throw("Artificial nuclear dripping requires CRPropa residual policy");
        if (fNuclei->GzkBundle() && !fNuclei->Find(key.first,key.second)) continue;
        fNuclei->Get(key.first,key.second);
        auto found=fChannels.find(key);
        if (found!=fChannels.end() && found->second.empty())
            Exception::Throw("Stable nuclear decay state also has a decay channel");
        auto& channels=fChannels[key];
        for (const auto& existing:channels) if (existing.code==channel.code)
            Exception::Throw("Duplicate nuclear decay channel");
        channel.restRate=1./(tau*units.second);
        if (!std::isfinite(channel.restRate) || channel.restRate<=0.)
            Exception::Throw("Unrepresentable nuclear decay rate");
        // A compatibility bundle contains states outside the selected root's
        // network. Validate their daughters only when the state is reached.
        if (fNuclei->GzkBundle()) { channels.push_back(channel); continue; }
        for (const auto& id:Products(key,channel.code)) {
            if (id.first>1 && (id.second==0 || id.second==id.first)) continue;
            Particle product(Nucleus,0.); fNuclei->SetParticle(product,id.first,id.second);
        }
        Particle previous(Nucleus,0.); fNuclei->SetParticle(previous,key.first,key.second);
        for (size_t beta=0;beta<2;++beta)
            for (int count=0;count<channel.code[beta]-'0';++count) {
                Particle next=previous;
                fNuclei->SetParticle(next,key.first,previous.AtomicNumber()+(beta?-1:1));
                if (previous.Mass()-next.Mass()<=Particle::Mass(Electron))
                    Exception::Throw("Nonpositive nuclear beta decay endpoint");
                previous=next;
            }
        channels.push_back(channel);
    }
}
const std::vector<NuclearDecayData::Channel>& NuclearDecayData::Get(unsigned int a,unsigned int z) const {
    if (fUnsupported.count({a,z})) Exception::Throw("Unsupported source decay mode for A"+
        Utils::ToString(a)+"Z"+Utils::ToString(z));
    const auto found=fChannels.find({a,z});
    if (found==fChannels.end()) Exception::Throw("Missing nuclear decay or stable-state data for A"+
        Utils::ToString(a)+"Z"+Utils::ToString(z));
    if (fNuclei->GzkBundle()) for (const auto& channel:found->second) {
        for (const auto& id:Products({a,z},channel.code))
            if (id.first>1 && id.second && id.second<id.first) fNuclei->Get(id.first,id.second);
        Particle previous(Nucleus,0.); fNuclei->SetParticle(previous,a,z);
        for (size_t beta=0;beta<2;++beta) for (int count=0;count<channel.code[beta]-'0';++count) {
            Particle next=previous;
            fNuclei->SetParticle(next,a,previous.AtomicNumber()+(beta?-1:1));
            if (previous.Mass()-next.Mass()<=Particle::Mass(Electron))
                Exception::Throw("Nonpositive nuclear beta decay endpoint");
            previous=next;
        }
    }
    return found->second;
}
std::set<TalysData::Identity> NuclearDecayData::ValidateNetwork(unsigned int a,unsigned int z,bool pd,bool gzk) const {
    std::set<TalysData::Identity> reached;
    std::vector<TalysData::Identity> pending{{a,z}};
    while (!pending.empty()) {
        const auto key=pending.back(); pending.pop_back();
        if (!reached.insert(key).second || key.first==1) continue;
        if (key.second==0 || key.second==key.first) {
            pending.emplace_back(1,key.second ? 1 : 0);
            continue;
        }
        const auto& isotope=fNuclei->Get(key.first,key.second);
        if (gzk) {
            if (key.second) pending.emplace_back(key.first-1,key.second-1);
            if (key.first>key.second) pending.emplace_back(key.first-1,key.second);
        }
        for (const auto& channel:Get(key.first,key.second)) {
            const auto ids=Products(key,channel.code); pending.insert(pending.end(),ids.begin(),ids.end());
        }
        if (pd && !isotope.terminal) {
            if (!isotope.total || isotope.channels.empty()) Exception::Throw("Missing reachable photodisintegration data with nuclear decays");
            for (const auto& channel:isotope.channels) {
                pending.push_back(channel.residual);
                for (size_t i=0;i<TalysData::Emitted.size();++i)
                    if (channel.multiplicity[i]) pending.push_back(TalysData::Emitted[i]);
            }
        }
    }
    return reached;
}
void NuclearDecayData::AppendResidual(const Particle& parent,unsigned int a,unsigned int z,
                                      double energy,std::vector<Particle>& products) const {
    if (a<1 || z>a || !std::isfinite(energy) || energy<=0.)
        Exception::Throw("Invalid nuclear spectator identity or energy");
    if (z==0 || z==a) {
        Particle nucleon=parent; nucleon.SetType(z ? Proton : Neutron);
        nucleon.Energy=energy/a;
        for (unsigned int i=0;i<a;++i) products.push_back(nucleon);
    } else {
        Particle residual=parent; fNuclei->SetParticle(residual,a,z);
        residual.Energy=energy; products.push_back(residual);
    }
}
NuclearDecay::NuclearDecay(std::shared_ptr<const NuclearDecayData> data):fData(std::move(data)) {
    if (!fData) Exception::Throw("NuclearDecay requires immutable decay data");
}
RandomInteraction* NuclearDecay::Clone() const { return new NuclearDecay(fData); }
double NuclearDecay::Rate(const Particle& p) const {
    if (p.Type!=Nucleus) return 0.;
    p.ValidateIdentity();
    if (p.MassModel()!=Particle::NuclearMassModel::TALYS ||
        p.Mass()!=fData->Nuclei().Get(p.MassNumber(),p.AtomicNumber()).mass)
        Exception::Throw("NuclearDecay requires the pinned tabulated nuclear mass");
    double rate=0.;
    for (const auto& channel:fData->Get(p.MassNumber(),p.AtomicNumber())) rate+=channel.restRate;
    rate*=p.Mass()/p.Energy; // Physical time: no comoving (1+z) conversion in mcray.
    if (!std::isfinite(rate) || rate<0.) Exception::Throw("Invalid nuclear decay lab rate");
    return rate;
}
void NuclearDecay::BetaDecay(Particle& residual,bool plus,std::vector<Particle>& products,Randomizer& random) const {
    const double gamma=residual.Energy/residual.Mass(), oldMass=residual.Mass();
    fData->Nuclei().SetParticle(residual,residual.MassNumber(),residual.AtomicNumber()+(plus?-1:1));
    residual.Energy=gamma*residual.Mass();
    const double me=Particle::Mass(Electron),q=oldMass-residual.Mass()-me;
    if (!(q>0.)) Exception::Throw("Nonpositive nuclear beta decay endpoint");
    // CRPropa's 51-node allowed beta spectrum, without Coulomb or recoil corrections.
    // Dimensionless energies keep the cumulative sum well scaled in natural units.
    std::array<double,51> cdf{};
    const double endpoint=q/me;
    for (size_t i=1;i<cdf.size();++i) {
        const double e=1.+i/50.*endpoint;
        cdf[i]=cdf[i-1]+e*std::sqrt(std::max(0.,e*e-1.))*std::pow(endpoint+1.-e,2);
    }
    const double u=random.Rand()*cdf.back();
    size_t hi=std::lower_bound(cdf.begin()+1,cdf.end(),u)-cdf.begin();
    if (hi>=cdf.size()) hi=cdf.size()-1;
    const double interval=cdf[hi]-cdf[hi-1];
    const double fraction=interval>0. ? (u-cdf[hi-1])/interval : 0.;
    const double e=me+(hi-1+fraction)/50.*q;
    const double cosine=2.*random.Rand()-1.;
    Particle electron=residual; electron.SetType(plus?Positron:Electron);
    electron.Energy=gamma*(e-std::sqrt(std::max(0.,e*e-me*me))*cosine);
    Particle neutrino=residual; neutrino.SetType(plus?NeutrinoE:NeutrinoAE);
    neutrino.Energy=gamma*(q+me-e)*(1.+cosine);
    products.push_back(electron); products.push_back(neutrino);
}
bool NuclearDecay::GetSecondaries(Particle& p,std::vector<Particle>& products,Randomizer& random) const {
    if (Rate(p)==0.) return false;
    const auto& channels=fData->Get(p.MassNumber(),p.AtomicNumber());
    double sum=0.; for (const auto& channel:channels) sum+=channel.restRate;
    double u=random.Rand()*sum;
    const auto* selected=&channels.back();
    for (const auto& channel:channels) { u-=channel.restRate; if (u<0.) { selected=&channel; break; } }
    Particle residual=p;
    for (int i=0;i<selected->code[0]-'0';++i) BetaDecay(residual,false,products,random);
    for (int i=0;i<selected->code[1]-'0';++i) BetaDecay(residual,true,products,random);
    const std::array<TalysData::Identity,3> emitted={{{4,2},{1,1},{1,0}}};
    for (size_t i=0;i<emitted.size();++i)
        for (int count=0;count<selected->code[i+2]-'0';++count) {
            const auto a=residual.MassNumber(),z=residual.AtomicNumber();
            const double energyPerA=residual.Energy/a;
            Particle fragment=residual; fData->Nuclei().SetParticle(fragment,emitted[i].first,emitted[i].second);
            fragment.Energy=energyPerA*emitted[i].first; products.push_back(fragment);
            const auto nextA=a-emitted[i].first,nextZ=z-emitted[i].second;
            if (nextA>1 && (nextZ==0 || nextZ==nextA)) {
                fData->AppendResidual(residual,nextA,nextZ,energyPerA*nextA,products);
                return true;
            }
            fData->Nuclei().SetParticle(residual,a-emitted[i].first,z-emitted[i].second);
            residual.Energy=energyPerA*(a-emitted[i].first);
        }
    products.push_back(residual);
    return true;
}
}
