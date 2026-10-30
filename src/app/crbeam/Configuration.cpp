#include "Configuration.h"
#include "MagneticConfiguration.h"
#include "SpatialConfiguration.h"
#include "Particle.h"
#include "Sha256.h"
#include "Utils.h"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <yaml-cpp/eventhandler.h>
#include <yaml-cpp/parser.h>

namespace fs = std::filesystem;
using cors::cmdline::CmdInfo;
using cors::cmdline::CmdLine;
namespace {
struct Binding {
  const char *path;
  const char *option;
  const char *kind;
  const char *choices;
};
const Binding bindings[] = {
    {"source.redshift", "--z", "positive", ""},
    {"source.particle.type", "--particle", "enum",
     "electron|positron|photon|neutrinoE|neutrinoM|neutrinoT|neutrinoAE|"
     "neutrinoAM|neutrinoAT|neutron|proton|nucleus"},
    {"source.particle.a", "--A", "unsigned", ""},
    {"source.particle.z", "--Z", "unsigned", ""},
    {"simulation.primaries", "--nparticles", "integer", ""},
    {"simulation.minimum_energy_eV", "--emin", "positive", ""},
    {"propagation.filament_stop_redshift", "--z-filament", "nonnegative", ""},
    {"backgrounds.cmb.model", "--monocmb", "switch", "blackbody|monochromatic"},
    {"backgrounds.cmb.fixed_temperature", "--fixedcmb", "boolean", ""},
    {"backgrounds.ebl.model", "--background", "indexed",
     "none|kneiske_elmag_best_fit|kneiske_elmag_minimal|inoue12_baseline|"
     "inoue12_lower|inoue12_upper|franceschini08|kneiske0309|kneiske1001|"
     "stecker05|stecker16_lower|stecker16_upper|franceschini17"},
    {"backgrounds.ebl.normalization", "--backgr-mult", "nonnegative", ""},
    {"backgrounds.ebl.maximum_photon_energy_eV", "--backgr-cut", "positive",
     ""},
    {"backgrounds.ebl.high_redshift_extension.power_law_index", "--eblPowerLow",
     "number", ""},
    {"backgrounds.ebl.high_redshift_extension.exponential_delta_z",
     "--eblDeltaZexp", "number", ""},
    {"backgrounds.additional.file", "--add-backgr", "path", ""},
    {"backgrounds.additional.density", "--add-backgr-phys", "switch",
     "comoving|physical"},
    {"nuclear.photodisintegration.model", "--photodisintegration", "mapped",
     "off|psb|talys"},
    {"nuclear.photodisintegration.missing_data", "--pd-missing-data", "mapped",
     "strict|crpropa"},
    {"nuclear.photodisintegration.photons.enabled",
     "--photodisintegration-photons", "onoff", ""},
    {"nuclear.photodisintegration.photons.yield_policy", "--pd-photon-yields",
     "mapped", "strict|crpropa"},
    {"nuclear.decay.enabled", "--nuclear-decay", "onoff", ""},
    {"nuclear.decay.photons", "--decay-photons", "onoff", ""},
    {"nuclear.pair_production.mode", "--nuclear-ppp", "mapped",
     "off|losses|full"},
    {"nuclear.photopion.enabled", "--nuclear-gzk", "onoff", ""},
    {"nuclear.photopion.residual_policy", "--nuclear-gzk-residuals", "mapped",
     "strict|crpropa"},
    {"electromagnetic.cascade", "--noEMcascade", "inverse", ""},
    {"magnetic_field.model", "--turbulentEGMF", "switch",
     "monochromatic|kolmogorov"},
    {"magnetic_field.strength_G", "--EGMF", "nonnegative", ""},
    {"magnetic_field.minimum_scale_Mpc", "--lminEGMF", "positive", ""},
    {"magnetic_field.maximum_scale_Mpc", "--lmaxEGMF", "positive", ""},
    {"magnetic_field.modes", "--nmodesEGMF", "integer", ""},
    {"magnetic_field.realization", "--randomEGMF", "switch",
     "shared|per_primary"},
    {"magnetic_field.deflection_accuracy", "--accEGMF", "integer", ""},
    {"sampling.thinning.alpha", "--thinning", "thinning_alpha", ""},
    {"sampling.pair_production_resampling", "--ppp-thinning", "positive", ""},
    {"output.directory", "--output", "path", ""},
    {"output.overwrite", "--overwrite", "boolean", ""},
    {"output.name_suffix", "--output-suffix", "string", ""},
    {"output.trajectory_log", "--tlog", "boolean", ""},
    {"logging.enabled", "--log", "boolean", ""},
    {"logging.level", "--log-level", "indexed",
     "error|warning|message|verbose"},
    {"diagnostics.optical_depth", "--print-tau", "boolean", ""}};
std::vector<std::string> Split(const std::string &s, char separator) {
  std::vector<std::string> out;
  std::istringstream stream(s);
  std::string item;
  while (std::getline(stream, item, separator))
    out.push_back(item);
  return out;
}
YAML::Node Get(const YAML::Node &root, const std::string &path) {
  YAML::Node current = root;
  for (const auto &key : Split(path, '.')) {
    if (!current.IsMap())
      return YAML::Node(YAML::NodeType::Undefined);
    const YAML::Node next = static_cast<const YAML::Node &>(current)[key];
    if (!next.IsDefined())
      return YAML::Node(YAML::NodeType::Undefined);
    current.reset(next);
  }
  return current;
}
void Put(YAML::Node root, const std::string &path, const YAML::Node &value) {
  auto keys = Split(path, '.');
  YAML::Node current = root;
  for (size_t i = 0; i + 1 < keys.size(); ++i)
    current.reset(current[keys[i]]);
  current[keys.back()] = YAML::Clone(value);
}
template <class T>
void Put(YAML::Node root, const std::string &path, const T &value) {
  Put(root, path, YAML::Node(value));
}
std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}
std::string Legacy(const std::string &s) {
  if (s == "psb")
    return "PSB";
  if (s == "talys")
    return "TALYS";
  if (s == "crpropa")
    return "CRPropa";
  std::string out = s;
  out[0] = std::toupper(out[0]);
  return out;
}
std::string Number(double value) {
  std::ostringstream s;
  s << std::setprecision(17) << value;
  return s.str();
}
[[noreturn]] void Invalid(const std::string &field, const std::string &why) {
  throw std::runtime_error(field + ": " + why);
}
unsigned long Unsigned(const YAML::Node &n, const std::string &path) {
  if (!n.IsScalar() || n.Tag() == "!")
    Invalid(path, "expected an unsigned integer");
  const auto s = n.Scalar();
  unsigned long out = 0;
  const auto r = std::from_chars(s.data(), s.data() + s.size(), out);
  if (r.ec != std::errc() || r.ptr != s.data() + s.size() || s.empty())
    Invalid(path, "invalid unsigned integer");
  return out;
}
double Numeric(const YAML::Node &n, const std::string &path) {
  if (!n.IsScalar() || n.Tag() == "!")
    Invalid(path, "expected a number");
  try {
    size_t used = 0;
    const double value = std::stod(n.Scalar(), &used);
    if (used != n.Scalar().size() || !std::isfinite(value))
      Invalid(path, "expected a finite number");
    return value;
  } catch (const std::exception &) {
    Invalid(path, "expected a finite number");
  }
}
bool Boolean(const YAML::Node &n, const std::string &path) {
  if (!n.IsScalar() || n.Tag() == "!" ||
      (n.Scalar() != "true" && n.Scalar() != "false"))
    Invalid(path, "expected true or false");
  return n.Scalar() == "true";
}
std::string StringValue(const YAML::Node &node, const std::string &path) {
  if (!node.IsScalar())
    Invalid(path, "expected a string");
  const auto text = node.Scalar();
  // yaml-cpp preserves implicit/quoted scalar tags but does not expose a
  // resolved scalar type. Reject implicit numbers/booleans in string fields.
  if (node.Tag() == "?") {
    const auto lower = Lower(text);
    if (lower == "true" || lower == "false" || lower == ".nan" ||
        lower == ".inf" || lower == "+.inf" || lower == "-.inf")
      Invalid(path, "expected a string; quote literal names when necessary");
    bool number = false;
    try {
      size_t used = 0;
      std::stod(text, &used);
      number = used == text.size();
    } catch (const std::exception &) {
    }
    if (number)
      Invalid(path, "expected a string; quote numeric file names");
  }
  return text;
}
std::set<std::string> Fields(bool versionTwo = false) {
  std::set<std::string> fields = {
      "schema_version",
      "source.distribution",
      "source.evolution_index",
      "source.spectrum.type",
      "source.spectrum.energy_eV",
      "source.spectrum.minimum_energy_eV",
      "source.spectrum.maximum_energy_eV",
      "source.spectrum.index",
      "simulation.seed",
      "simulation.threads",
      "simulation.batch_size",
      "logging.thread",
      "backgrounds.ebl.high_redshift_extension.enabled",
      "backgrounds.ebl.high_redshift_extension.constant_delta_z"};
  for (const auto &b : bindings)
    fields.insert(b.path);
  if (versionTwo)
    for (const char *field :
         {"source.position_Mpc", "source.angular_distribution.model",
          "source.angular_distribution.axis",
          "source.angular_distribution.half_angle_deg", "source.angular_distribution.kappa",
          "propagation.end_redshift", "observer.model", "observer.center_Mpc",
          "observer.radius_Mpc", "observer.stop_when",
          "output.primary_direction", "sampling.thinning.model",
          "sampling.thinning.threshold_eV"})
      fields.insert(field);
  return fields;
}
void CheckKeys(const YAML::Node &n, const std::string &path,
               const std::set<std::string> &fields, bool magneticV2 = false) {
  if (!n.IsMap())
    Invalid(path, "expected a mapping");
  std::set<std::string> seen;
  for (const auto &entry : n) {
    if (!entry.first.IsScalar())
      Invalid(path, "mapping keys must be strings");
    const auto key = entry.first.Scalar();
    const std::string child = path.empty() ? key : path + "." + key;
    if (!seen.insert(key).second)
      Invalid(child, "duplicate key");
    if (magneticV2 &&
        (child == "magnetic_field" || child == "source.position_Mpc" ||
         child == "source.angular_distribution.axis" ||
         child == "observer.center_Mpc"))
      continue;
    if (fields.count(child)) {
      if (entry.second.IsMap() || entry.second.IsSequence())
        Invalid(child, "expected a scalar");
    } else {
      bool parent = false;
      for (const auto &field : fields)
        if (field.rfind(child + ".", 0) == 0)
          parent = true;
      if (!parent)
        Invalid(child, "unknown configuration key");
      CheckKeys(entry.second, child, fields, magneticV2);
    }
  }
}
YAML::Node VectorNode(const mcray::SpatialVector &vector) {
  YAML::Node node(YAML::NodeType::Sequence);
  for (auto v : vector)
    node.push_back(static_cast<double>(v));
  return node;
}
mcray::SpatialVector SpatialVectorValue(const YAML::Node &input,
                                        const std::string &path,
                                        mcray::SpatialVector fallback) {
  const auto node = Get(input, path);
  if (!node)
    return fallback;
  if (!node.IsSequence() || node.size() != 3)
    Invalid(path, "expected exactly three finite numbers");
  for (size_t i = 0; i < 3; ++i)
    fallback[i] = Numeric(node[i], path);
  return fallback;
}
void NormalizeSpatial(YAML::Node effective, const YAML::Node &input) {
  const auto position = SpatialVectorValue(input, "source.position_Mpc", {});
  auto axis =
      SpatialVectorValue(input, "source.angular_distribution.axis", {0, 0, 1});
  const auto scale =
      std::max({std::abs(axis[0]), std::abs(axis[1]), std::abs(axis[2])});
  if (!(scale > 0))
    Invalid("source.angular_distribution.axis", "must not be zero");
  // Keep already-unit serialized directions stable when replaying config.yaml.
  const auto length = std::hypot(axis[0], axis[1], axis[2]);
  if (!std::isfinite(length) ||
      std::abs(length - 1) > 8 * std::numeric_limits<double>::epsilon())
    axis = mcray::SourceGeometry::NormalizeAxis(axis);
  auto text = [&](const char *path, const char *fallback) {
    const auto node = Get(input, path);
    return node ? StringValue(node, path) : std::string(fallback);
  };
  auto number = [&](const char *path, double fallback) {
    const auto node = Get(input, path);
    return node ? Numeric(node, path) : fallback;
  };
  const auto model = text("source.angular_distribution.model", "pencil");
  if (model != "pencil" && model != "cone" && model != "vmf")
    Invalid("source.angular_distribution.model", "choices: pencil|cone|vmf");
  const auto kappaNode = Get(input, "source.angular_distribution.kappa");
  double kappa = 0;
  if (model == "vmf") {
    if (!Get(input, "source.angular_distribution.half_angle_deg"))
      Invalid("source.angular_distribution.half_angle_deg", "required for vmf");
    if (!kappaNode)
      Invalid("source.angular_distribution.kappa", "required for vmf");
    kappa = Numeric(kappaNode, "source.angular_distribution.kappa");
    if (kappa < 0)
      Invalid("source.angular_distribution.kappa", "must be nonnegative");
  } else if (kappaNode)
    Invalid("source.angular_distribution.kappa", "requires vmf model");
  const double half = number("source.angular_distribution.half_angle_deg", 0);
  if (half < 0 || half > 180 || (model == "pencil" && half != 0))
    Invalid("source.angular_distribution.half_angle_deg",
            "require 0..180 for cone/vmf, zero for pencil");
  const double end = number("propagation.end_redshift", 0);
  const auto sourceZ = Get(effective, "source.redshift").as<double>();
  if (end < 0 || (sourceZ > 0 && end >= sourceZ) || (sourceZ <= 0 && end != 0))
    Invalid("propagation.end_redshift",
            "must be nonnegative and below source redshift");
  const auto observer = text("observer.model", "none");
  if (observer != "none" && observer != "sphere")
    Invalid("observer.model", "choices: none|sphere");
  YAML::Node observerNode;
  observerNode["model"] = observer;
  if (observer == "none") {
    if (Get(input, "observer.center_Mpc") ||
        Get(input, "observer.radius_Mpc") || Get(input, "observer.stop_when"))
      Invalid("observer", "center/radius/stop_when require sphere");
  } else {
    const auto center = SpatialVectorValue(input, "observer.center_Mpc", {});
    if (!Get(input, "observer.radius_Mpc"))
      Invalid("observer.radius_Mpc", "required for sphere");
    const double radius = number("observer.radius_Mpc", 0);
    if (!(radius > 0))
      Invalid("observer.radius_Mpc", "must be positive");
    const auto stop = text("observer.stop_when", "inside");
    if (stop != "inside" && stop != "outside")
      Invalid("observer.stop_when", "choices: inside|outside");
    observerNode["center_Mpc"] = VectorNode(center);
    observerNode["radius_Mpc"] = radius;
    observerNode["stop_when"] = stop;
  }
  Put(effective, "source.position_Mpc", VectorNode(position));
  Put(effective, "source.angular_distribution.model", model);
  Put(effective, "source.angular_distribution.axis", VectorNode(axis));
  Put(effective, "source.angular_distribution.half_angle_deg", half);
  if (model == "vmf")
    Put(effective, "source.angular_distribution.kappa", kappa);
  Put(effective, "propagation.end_redshift", end);
  effective["observer"] = observerNode;
  const auto direction = Get(input, "output.primary_direction");
  Put(effective, "output.primary_direction",
      direction ? Boolean(direction, "output.primary_direction") : false);
  // Construct the same runtime settings now, so invalid unit conversions fail
  // before output.
  const SpatialConfiguration checked(effective);
}
// Keep the version-1 effective schema unchanged. Version 2 explicitly names
// the policy, so config.yaml replay cannot accidentally select another model.
void NormalizeThinning(YAML::Node effective) {
  const std::string prefix = "sampling.thinning.";
  const auto modelNode = Get(effective, prefix + "model");
  const std::string model = modelNode ? StringValue(modelNode, prefix + "model") : "energy";
  if (model != "energy" && model != "photon")
    Invalid(prefix + "model", "choices: energy|photon");
  const double alpha = Numeric(Get(effective, prefix + "alpha"), prefix + "alpha");
  if (alpha < 0 || alpha > (model == "photon" ? 10 : 1))
    Invalid(prefix + "alpha", model == "photon" ? "require 0..10" : "require 0..1");
  const auto threshold = Get(effective, prefix + "threshold_eV");
  if (model == "photon") {
    if (!threshold)
      Invalid(prefix + "threshold_eV", "required for photon model");
    const double value = Numeric(threshold, prefix + "threshold_eV");
    const double internal = value * mcray::units.eV;
    if (!(value > 0) || !(internal > 0) || !std::isfinite(internal))
      Invalid(prefix + "threshold_eV", "must be positive and representable in internal units");
    Put(effective, prefix + "threshold_eV", value);
  } else if (threshold)
    Invalid(prefix + "threshold_eV", "requires photon model");
  Put(effective, prefix + "model", model);
  Put(effective, prefix + "alpha", alpha);
}
bool ThinningCli(CmdLine& cli, YAML::Node effective) {
  const bool extended = cli.has_param("--thinning-model") || cli.has_param("--thinning-threshold-eV");
  for (const auto* option : {"--thinning-model", "--thinning-threshold-eV", "--thinning"}) {
    if (!cli.has_param(option))
      continue;
    const auto parameter = cli.get_param(option);
    if (!parameter.value || !*parameter.value)
      Invalid(option, "requires a nonempty argument");
    if (std::string(option) == "--thinning-model")
      Put(effective, "sampling.thinning.model", std::string(parameter.value));
    else
      Put(effective, std::string(option) == "--thinning" ? "sampling.thinning.alpha" :
          "sampling.thinning.threshold_eV", Numeric(YAML::Node(parameter.value), option));
  }
  if (!extended)
    return false;
  NormalizeThinning(effective);
  return true;
}
bool SpatialCli(CmdLine &cli, YAML::Node effective) {
  struct Option {
    const char *option;
    const char *field;
    const char *kind;
  };
  const Option options[] = {
      {"--source-position-Mpc", "source.position_Mpc", "vector"},
      {"--source-angular-model", "source.angular_distribution.model", "text"},
      {"--source-axis", "source.angular_distribution.axis", "vector"},
      {"--source-half-angle-deg", "source.angular_distribution.half_angle_deg",
       "number"},
      {"--source-vmf-kappa", "source.angular_distribution.kappa", "number"},
      {"--end-redshift", "propagation.end_redshift", "number"},
      {"--observer", "observer.model", "text"},
      {"--observer-center-Mpc", "observer.center_Mpc", "vector"},
      {"--observer-radius-Mpc", "observer.radius_Mpc", "number"},
      {"--observer-stop-when", "observer.stop_when", "text"},
      {"--output-primary-direction", "output.primary_direction", "flag"}};
  YAML::Node input(YAML::NodeType::Map);
  bool any = false;
  for (const auto &option : options) {
    if (!cli.has_param(option.option))
      continue;
    any = true;
    const auto parameter = cli.get_param(option.option);
    if (std::string(option.kind) == "flag") {
      Put(input, option.field, true);
      continue;
    }
    if (!parameter.value || !*parameter.value)
      Invalid(option.option, "requires a nonempty argument");
    const std::string raw = parameter.value;
    if (std::string(option.kind) == "vector") {
      if (std::count(raw.begin(), raw.end(), ',') != 2)
        Invalid(option.option, "expected x,y,z");
      YAML::Node vector(YAML::NodeType::Sequence);
      for (auto item : Split(raw, ',')) {
        const auto first = item.find_first_not_of(" \t\r\n"),
                   last = item.find_last_not_of(" \t\r\n");
        if (first == std::string::npos)
          Invalid(option.option, "empty vector component");
        item = item.substr(first, last - first + 1);
        vector.push_back(Numeric(YAML::Node(item), option.option));
      }
      if (vector.size() != 3)
        Invalid(option.option, "expected exactly three coordinates");
      Put(input, option.field, vector);
    } else if (std::string(option.kind) == "number")
      Put(input, option.field, Numeric(YAML::Node(raw), option.option));
    else
      Put(input, option.field, raw);
  }
  if (!any)
    return false;
  if (!cli.has_param("--observer") && (cli.has_param("--observer-center-Mpc") ||
                                       cli.has_param("--observer-radius-Mpc") ||
                                       cli.has_param("--observer-stop-when")))
    Invalid("observer", "sphere parameters require --observer sphere");
  NormalizeSpatial(effective, input);
  return true;
}
// Check events before constructing nodes, so aliases cannot create recursive
// trees.
class PlainYaml : public YAML::EventHandler {
  void Check(const std::string &tag, YAML::anchor_t anchor) {
    if (anchor || (tag != "?" && tag != "!" && !tag.empty()))
      Invalid("YAML", "anchors and explicit tags are unsupported");
  }

public:
  void OnDocumentStart(const YAML::Mark &) override {}
  void OnDocumentEnd() override {}
  void OnNull(const YAML::Mark &, YAML::anchor_t a) override { Check("", a); }
  void OnAlias(const YAML::Mark &, YAML::anchor_t) override {
    Invalid("YAML", "aliases are unsupported");
  }
  void OnScalar(const YAML::Mark &, const std::string &t, YAML::anchor_t a,
                const std::string &) override {
    Check(t, a);
  }
  void OnSequenceStart(const YAML::Mark &, const std::string &t,
                       YAML::anchor_t a, YAML::EmitterStyle::value) override {
    Check(t, a);
  }
  void OnSequenceEnd() override {}
  void OnMapStart(const YAML::Mark &, const std::string &t, YAML::anchor_t a,
                  YAML::EmitterStyle::value) override {
    Check(t, a);
  }
  void OnMapEnd() override {}
};
YAML::Node Load(const std::string &path) {
  std::ifstream input(path);
  if (!input)
    Invalid(path, "cannot read configuration");
  YAML::Parser parser(input);
  PlainYaml events;
  if (!parser.HandleNextDocument(events) || parser.HandleNextDocument(events))
    Invalid(path, "expected one YAML document");
  return YAML::LoadFile(path);
}
void Overlay(YAML::Node target, const YAML::Node &source) {
  for (const auto &item : source) {
    const auto key = item.first.as<std::string>();
    if (item.second.IsMap())
      Overlay(target[key], item.second);
    else
      target[key] = YAML::Clone(item.second);
  }
}
std::string Absolute(const std::string &path, const fs::path &base) {
  fs::path p(path);
  return fs::absolute(p.is_absolute() ? p : base / p)
      .lexically_normal()
      .string();
}
void Emit(const fs::path &path, const YAML::Node &node) {
  YAML::Emitter emitter;
  emitter.SetDoublePrecision(17);
  emitter << node;
  if (!emitter.good())
    throw std::runtime_error("Cannot serialize YAML: " +
                             emitter.GetLastError());
  const auto tmp = path.string() + ".partial";
  std::ofstream out(tmp);
  out << emitter.c_str() << '\n';
  out.close();
  if (!out)
    throw std::runtime_error("Cannot write " + tmp);
  fs::rename(tmp, path);
}
std::string Now() {
  const auto now = std::time(nullptr);
  std::tm utc{};
#ifdef _WIN32
  gmtime_s(&utc, &now);
#else
  gmtime_r(&now, &utc);
#endif
  std::ostringstream out;
  out << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
  return out.str();
}
void Inventory(YAML::Node result, const fs::path &root,
               const fs::path &output) {
  if (!fs::is_directory(root))
    return;
  const auto excluded = fs::weakly_canonical(output);
  const auto base = fs::absolute(root).lexically_normal();
  std::set<fs::path> visited{fs::weakly_canonical(base)};
  fs::recursive_directory_iterator it(
      base, fs::directory_options::follow_directory_symlink),
      end;
  for (; it != end; ++it) {
    if (it->is_directory()) {
      const auto real = fs::weakly_canonical(it->path());
      if (real == excluded || !visited.insert(real).second)
        it.disable_recursion_pending();
    } else if (it->is_regular_file()) {
      result[it->path().lexically_relative(base).generic_string()] =
          Configuration::HashFile(it->path().string());
    }
  }
}
} // namespace

Configuration::Configuration(int argc, char **argv, CmdInfo *commands) {
  try {
    for (int i = 0; i < argc; ++i)
      fOriginal.emplace_back(argv[i]);
    int file = -1;
    for (int i = 1; i < argc; ++i) {
      const CmdInfo *known = nullptr;
      for (const auto *c = commands; c->code; ++c)
        if (std::string(argv[i]) == c->name ||
            std::string(argv[i]) == c->short_name) {
          known = c;
          break;
        }
      if (known) {
        if (known->flags & CmdInfo::FLAG_ARGUMENT)
          ++i;
        continue;
      }
      const auto ext = fs::path(argv[i]).extension().string();
      if (ext == ".yaml" || ext == ".yml") {
        file = i;
        break;
      }
    }
    if (file >= 0 && argc != 2)
      Invalid("YAML", "cannot mix a YAML file and CLI options");
    const bool yaml = file >= 0;
    std::vector<char *> defaults{argv[0], nullptr};
    // Preserve legacy CLI interpretation; YAML emits validated equivalents
    // into the same command model consumed by CRbeam process setup.
    CmdLine cli(yaml ? 1 : argc, yaml ? defaults.data() : argv, commands);
    fEffective["schema_version"] = 1;
    for (const auto &b : bindings) {
      auto p = cli.get_param(b.option);
      std::string value = static_cast<const char *>(p);
      const std::string kind = b.kind;
      if (kind == "boolean" || kind == "inverse")
        Put(fEffective, b.path, kind == "inverse" ? !bool(p) : bool(p));
      else if (kind == "onoff")
        Put(fEffective, b.path, value == "On");
      else if (kind == "switch")
        Put(fEffective, b.path, Split(b.choices, '|')[bool(p) ? 1 : 0]);
      else if (kind == "mapped")
        Put(fEffective, b.path, Lower(value));
      else if (kind == "indexed") {
        const int index = static_cast<int>(p);
        const auto choices = Split(b.choices, '|');
        if (index < 0 || static_cast<size_t>(index) >= choices.size())
          Invalid(b.path, "invalid legacy enum index");
        Put(fEffective, b.path, choices[index]);
      } else if (kind == "path")
        Put(fEffective, b.path,
            value.empty() ? YAML::Node(YAML::NodeType::Null)
                          : YAML::Node(value));
      else if (kind == "string" || kind == "enum")
        Put(fEffective, b.path, value);
      else if (kind == "unsigned")
        Put(fEffective, b.path, static_cast<unsigned int>(p));
      else if (kind == "integer")
        Put(fEffective, b.path, static_cast<int>(p));
      else
        Put(fEffective, b.path, static_cast<double>(p));
    }
    if (Get(fEffective, "source.particle.type").as<std::string>().empty()) {
      const int type = static_cast<int>(cli.get_param("--primary"));
      if (type < mcray::Electron || type > mcray::Proton)
        Invalid("--primary", "invalid primary particle type");
      Put(fEffective, "source.particle.type",
          mcray::Particle::Name(static_cast<mcray::ParticleType>(type)));
    }
    const double power = cli.get_param("--power"),
                 evolution = cli.get_param("--m_z");
    Put(fEffective, "source.distribution",
        evolution == 1000 ? "point" : "continuous");
    Put(fEffective, "source.evolution_index",
        evolution == 1000 ? 0. : evolution);
    Put(fEffective, "source.spectrum.type",
        power > 0 ? "power_law" : "monoenergetic");
    if (power > 0) {
      Put(fEffective, "source.spectrum.index", power);
      Put(fEffective, "source.spectrum.minimum_energy_eV",
          std::max(static_cast<double>(cli.get_param("--emin_source")),
                   static_cast<double>(cli.get_param("--emin"))));
      Put(fEffective, "source.spectrum.maximum_energy_eV",
          static_cast<double>(cli.get_param("--emax")));
    } else
      Put(fEffective, "source.spectrum.energy_eV",
          static_cast<double>(cli.get_param("--emax")));
    fAutomaticSeed = !cli.has_param("--seed");
    Put(fEffective, "simulation.seed",
        cli.has_param("--seed")
            ? YAML::Node(static_cast<unsigned int>(cli.get_param("--seed")))
            : YAML::Node("auto"));
    const int batch = cli.get_param("--batch"),
              thread = cli.get_param("--log-thread");
    Put(fEffective, "simulation.batch_size",
        batch ? YAML::Node(batch) : YAML::Node("auto"));
    Put(fEffective, "simulation.threads", "auto");
    Put(fEffective, "logging.thread",
        thread < 0 ? YAML::Node(YAML::NodeType::Null) : YAML::Node(thread));
    const double dz = cli.get_param("--eblDeltaZconst");
    Put(fEffective, "backgrounds.ebl.high_redshift_extension.enabled", dz >= 0);
    Put(fEffective, "backgrounds.ebl.high_redshift_extension.constant_delta_z",
        std::max(0., dz));
    if (!yaml) {
      fArguments = fOriginal;
      for (const char* option : {"--magnetic-field","--magnetic-field-file","--magnetic-field-boundary",
          "--source-position-Mpc","--source-angular-model","--source-axis","--source-half-angle-deg","--source-vmf-kappa",
          "--end-redshift","--observer","--observer-center-Mpc","--observer-radius-Mpc",
          "--observer-stop-when","--output-primary-direction",
          "--thinning-model","--thinning-threshold-eV"})
        if (std::count(fOriginal.begin()+1,fOriginal.end(),std::string(option)) > 1)
          Invalid(std::string(option).rfind("--magnetic-", 0) == 0 ? "magnetic_field" :
                  std::string(option).rfind("--thinning-", 0) == 0 ? "sampling.thinning" : "spatial_geometry",
                  std::string("duplicate option ") + option);
      const auto field = MagneticFieldFromCli(cli, fEffective["magnetic_field"]);
      if (field.IsDefined()) {
        fEffective["schema_version"] = 2;
        fEffective["magnetic_field"] = field;
        MagneticLegacyArguments(field,fArguments);
      }
      const bool spatialCli = SpatialCli(cli, fEffective);
      const bool thinningCli = ThinningCli(cli, fEffective);
      if ((spatialCli || thinningCli) && fEffective["schema_version"].as<int>()==1) {
        const YAML::Node missing(YAML::NodeType::Undefined);
        fEffective["magnetic_field"]=NormalizeMagneticField(missing,fEffective["magnetic_field"],fs::current_path().string());
        fEffective["schema_version"]=2;
      }
      if (fEffective["schema_version"].as<int>() == 2)
        NormalizeThinning(fEffective);
      else if (cli.has_param("--thinning")) {
        // Enforce the legacy policy's existing range in Release too, without
        // adding version-2 fields to an ordinary legacy configuration.
        const double alpha = Get(fEffective, "sampling.thinning.alpha").as<double>();
        if (alpha < 0 || alpha > 1)
          Invalid("--thinning", "require 0..1 for energy model");
      }
    } else {
      fInputFile = Absolute(argv[file], fs::current_path());
      fInputHash = HashFile(fInputFile);
      auto input = Load(fInputFile);
      if (!input["schema_version"]) Invalid("schema_version", "supported versions are 1 and 2");
      const auto version = Unsigned(input["schema_version"], "schema_version");
      if (version != 1 && version != 2) Invalid("schema_version", "supported versions are 1 and 2");
      CheckKeys(input,"",Fields(version == 2),version == 2);
      const auto legacyMagnetic = YAML::Clone(fEffective["magnetic_field"]);
      YAML::Node magnetic;
      if (version == 2) {
        magnetic = NormalizeMagneticField(static_cast<const YAML::Node&>(input)["magnetic_field"],
            legacyMagnetic,fs::path(fInputFile).parent_path().string());
        input.remove("magnetic_field");
      }
      CheckKeys(input, "", Fields(version == 2),version == 2);
      const auto spectrum = Get(input, "source.spectrum");
      if (spectrum && spectrum["type"] &&
          spectrum["type"].as<std::string>() == "power_law")
        fEffective["source"]["spectrum"] = YAML::Node(YAML::NodeType::Map);
      Overlay(fEffective, input);
      if (version == 2) {
        fEffective["magnetic_field"] = magnetic;
        NormalizeSpatial(fEffective,input);
        NormalizeThinning(fEffective);
      }
      fArguments.push_back(argv[0]);
      auto add = [&](const std::string &key, const std::string &value) {
        fArguments.push_back(key);
        fArguments.push_back(value);
      };
      const auto base = fs::path(fInputFile).parent_path();
      for (const auto &b : bindings) {
        if (version == 2 && std::string(b.path).rfind("magnetic_field.",0) == 0) continue;
        auto n = Get(fEffective, b.path);
        const std::string kind = b.kind;
        if (kind == "boolean" || kind == "inverse" || kind == "onoff") {
          bool enabled = Boolean(n, b.path);
          if (kind == "inverse")
            enabled = !enabled;
          if (kind == "onoff")
            add(b.option, enabled ? "On" : "Off");
          else if (enabled)
            fArguments.push_back(b.option);
        } else if (kind == "enum" || kind == "indexed" || kind == "mapped" ||
                   kind == "switch") {
          if (!n.IsScalar())
            Invalid(b.path, "expected an enum string");
          const auto value = StringValue(n, b.path);
          const auto choices = Split(b.choices, '|');
          const auto it = std::find(choices.begin(), choices.end(), value);
          if (it == choices.end())
            Invalid(b.path, std::string("choices: ") + b.choices);
          if (kind == "switch") {
            if (it != choices.begin())
              fArguments.push_back(b.option);
          } else
            add(b.option, kind == "indexed"
                              ? std::to_string(it - choices.begin())
                          : kind == "mapped" ? Legacy(value)
                                             : value);
        } else if (kind == "path") {
          if (!n.IsNull()) {
            const auto name = StringValue(n, b.path);
            if (name.empty())
              Invalid(b.path, "expected a nonempty path or null");
            const auto path = Absolute(name, base);
            Put(fEffective, b.path, path);
            add(b.option, path);
          }
        } else if (kind == "string") {
          add(b.option, StringValue(n, b.path));
        } else if (kind == "unsigned" || kind == "integer") {
          const auto value = Unsigned(n, b.path);
          const auto limit =
              kind == "unsigned"
                  ? std::numeric_limits<unsigned int>::max()
                  : static_cast<unsigned long>(std::numeric_limits<int>::max());
          if (value > limit || (kind == "integer" && !value))
            Invalid(b.path, "integer outside supported range");
          add(b.option, std::to_string(value));
        } else {
          const double v = Numeric(n, b.path);
          if ((kind == "positive" && v <= 0) ||
              (kind == "nonnegative" && v < 0) ||
              (kind == "thinning_alpha" && version == 1 && (v < 0 || v > 1)))
            Invalid(b.path, "outside allowed numeric range");
          add(b.option, Number(v));
        }
      }
      const auto particle =
          Get(fEffective, "source.particle.type").as<std::string>();
      if (particle == "nucleus") {
        if (!Get(input, "source.particle.a") ||
            !Get(input, "source.particle.z") ||
            Unsigned(Get(input, "source.particle.a"), "source.particle.a") <
                2 ||
            !Unsigned(Get(input, "source.particle.z"), "source.particle.z"))
          Invalid("source.particle", "nucleus requires explicit a>=2 and z>0");
      } else if (Get(input, "source.particle.a") ||
                 Get(input, "source.particle.z"))
        Invalid("source.particle", "a/z apply only to nucleus");
      const auto dist =
          Get(fEffective, "source.distribution").as<std::string>();
      const double index = Numeric(Get(fEffective, "source.evolution_index"),
                                   "source.evolution_index");
      if (dist != "point" && dist != "continuous")
        Invalid("source.distribution", "choices: point|continuous");
      if (dist == "continuous" && index == 1000)
        Invalid("source.evolution_index",
                "1000 is reserved by the legacy source sampler");
      add("--m_z", dist == "point" ? "1000" : Number(index));
      const auto spec =
          Get(fEffective, "source.spectrum.type").as<std::string>();
      if (spec != "monoenergetic" && spec != "power_law")
        Invalid("source.spectrum.type", "choices: monoenergetic|power_law");
      const double cutoff =
          Numeric(Get(fEffective, "simulation.minimum_energy_eV"),
                  "simulation.minimum_energy_eV");
      if (spec == "monoenergetic") {
        for (const auto *key :
             {"index", "minimum_energy_eV", "maximum_energy_eV"})
          if (spectrum && spectrum[key])
            Invalid("source.spectrum", "power-law fields require power_law");
        const double e = Numeric(Get(fEffective, "source.spectrum.energy_eV"),
                                 "source.spectrum.energy_eV");
        if (e < cutoff)
          Invalid("source.spectrum.energy_eV", "must be >= propagation cutoff");
        add("--emax", Number(e));
      } else {
        if (spectrum["energy_eV"])
          Invalid("source.spectrum", "energy_eV applies only to monoenergetic");
        const double lo = Numeric(
                         Get(fEffective, "source.spectrum.minimum_energy_eV"),
                         "source.spectrum.minimum_energy_eV"),
                     hi = Numeric(
                         Get(fEffective, "source.spectrum.maximum_energy_eV"),
                         "source.spectrum.maximum_energy_eV"),
                     p = Numeric(Get(fEffective, "source.spectrum.index"),
                                 "source.spectrum.index");
        if (lo < cutoff || hi <= lo || p <= 0)
          Invalid("source.spectrum",
                  "require cutoff<=minimum<maximum and index>0");
        add("--emin_source", Number(lo));
        add("--emax", Number(hi));
        add("--power", Number(p));
      }
      for (const auto *key : {"seed", "threads", "batch_size"}) {
        auto n = Get(fEffective, std::string("simulation.") + key);
        if (n.IsScalar() && n.Scalar() == "auto")
          continue;
        const auto v = Unsigned(n, std::string("simulation.") + key);
        if (std::string(key) == "seed") {
          fHasSeed = true;
          fSeed = v;
          fAutomaticSeed = false;
        } else {
          if (!v ||
              v > static_cast<unsigned long>(std::numeric_limits<int>::max()))
            Invalid(key, "expected positive integer or auto");
          if (std::string(key) == "threads")
            Put(fEffective, "simulation.threads", v);
          else
            add("--batch", std::to_string(v));
        }
      }
      const auto logThread = Get(fEffective, "logging.thread");
      if (!logThread.IsNull()) {
        const auto v = Unsigned(logThread, "logging.thread");
        if (v > static_cast<unsigned long>(std::numeric_limits<int>::max()))
          Invalid("logging.thread", "out of range");
        add("--log-thread", std::to_string(v));
      }
      const auto prefix =
          std::string("backgrounds.ebl.high_redshift_extension.");
      const double plateau =
          Numeric(Get(fEffective, prefix + "constant_delta_z"),
                  prefix + "constant_delta_z");
      if (plateau < 0)
        Invalid(prefix + "constant_delta_z", "must be nonnegative");
      add("--eblDeltaZconst",
          Boolean(Get(fEffective, prefix + "enabled"), prefix + "enabled")
              ? Number(plateau)
              : "-1");
      if (version == 2) MagneticLegacyArguments(fEffective["magnetic_field"],fArguments);
      if (version == 1 && Numeric(Get(fEffective, "magnetic_field.strength_G"),
                  "magnetic_field.strength_G") > 0 &&
          Get(fEffective, "magnetic_field.model").as<std::string>() ==
              "kolmogorov" &&
          Numeric(Get(fEffective, "magnetic_field.minimum_scale_Mpc"),
                  "magnetic_field.minimum_scale_Mpc") >=
              Numeric(Get(fEffective, "magnetic_field.maximum_scale_Mpc"),
                      "magnetic_field.maximum_scale_Mpc"))
        Invalid("magnetic_field", "minimum scale must be below maximum");
    }
    for (auto &arg : fArguments)
      fArgv.push_back(arg.data());
    fArgv.push_back(nullptr);
  } catch (const std::exception &e) {
    Utils::Exception::Throw(std::string("Configuration: ") + e.what());
  }
}

std::string Configuration::NuclearTables() {
  return (fs::path(Utils::tables_dir) / "nuclear").string();
}
void Configuration::SetAdditionalTable(const std::string &path) {
  Put(fEffective, "backgrounds.additional.file", path);
}
void Configuration::ResolveExecution(unsigned long seed, int threads, int batch,
                                     const std::string &output, unsigned a,
                                     unsigned z) {
  Put(fEffective, "simulation.seed", seed);
  Put(fEffective, "simulation.threads", threads);
  Put(fEffective, "simulation.batch_size", batch);
  Put(fEffective, "output.directory",
      fs::absolute(output).lexically_normal().string());
  if (a) {
    Put(fEffective, "source.particle.a", a);
    Put(fEffective, "source.particle.z", z);
  } else {
    fEffective["source"]["particle"].remove("a");
    fEffective["source"]["particle"].remove("z");
  }
}
std::string Configuration::HashFile(const std::string &path) {
  return Sha256File(path);
}
void Configuration::VerifyBundle(const std::string &directory,
                                 bool requireHashes) {
  try {
    if (requireHashes) {
      const auto path = fs::path(directory) / "manifest.json";
      if (!fs::exists(path))
        Invalid("nuclear tables", "missing manifest: " + path.string());
      const auto manifest = YAML::LoadFile(path.string());
      if (!manifest["payload_sha256"].IsMap() ||
          manifest["payload_sha256"].size() == 0)
        Invalid("nuclear tables",
                "standard installed bundle requires payload hashes");
    }
    for (const auto *name :
         {"manifest.json", "decay_manifest.json", "photon_manifest.json"}) {
      const auto path = fs::path(directory) / name;
      if (!fs::exists(path))
        continue;
      const auto manifest = YAML::LoadFile(path.string());
      if (!manifest["payload_sha256"])
        continue;
      if (!manifest["payload_sha256"].IsMap())
        Invalid(name, "payload_sha256 must be a mapping");
      for (const auto &row : manifest["payload_sha256"]) {
        const fs::path relative = row.first.as<std::string>();
        if (relative.is_absolute())
          Invalid(name, "absolute manifest path");
        for (const auto &part : relative)
          if (part == "..")
            Invalid(name, "parent traversal in manifest");
        if (HashFile((fs::path(directory) / relative).string()) !=
            row.second.as<std::string>())
          Invalid(name, "checksum mismatch: " + relative.string());
      }
    }
  } catch (const std::exception &e) {
    Utils::Exception::Throw(e.what());
  }
}
void Configuration::WriteConfig(const std::string &directory) const {
  try {
    Emit(fs::path(directory) / "config.yaml", fEffective);
  } catch (const std::exception &e) {
    Utils::Exception::Throw(e.what());
  }
}
YAML::Node
Configuration::Provenance(const std::string &nuclearDirectory) const {
  YAML::Node out;
  out["schema_version"] = 1;
  out["status"] = "running";
  out["started_utc"] = Now();
  out["input"]["mode"] = IsYaml() ? "yaml" : "cli";
  out["execution"]["seed_origin"] = fAutomaticSeed ? "automatic" : "explicit";
  if (IsYaml()) {
    out["input"]["file"] = fInputFile;
    out["input"]["sha256"] = fInputHash;
  }
  out["input"]["executable"] = boost::dll::program_location().string();
  out["build"]["revision"] = MCRAY_REVISION;
  out["build"]["tracked_patch_sha256"] = MCRAY_PATCH_SHA256;
  out["build"]["application_sources_sha256"] = MCRAY_CONFIG_SHA256;
  out["build"]["compiler"] = MCRAY_COMPILER;
  out["build"]["type"] = MCRAY_BUILD_TYPE;
  out["build"]["executable_sha256"] =
      HashFile(boost::dll::program_location().string());
  out["tables"]["root"] =
      fs::absolute(Utils::tables_dir).lexically_normal().string();
  out["tables"]["inventory_scope"] =
      "all files installed under the table root; not a claim all were used";
  out["tables"]["installed_sha256"] = YAML::Node(YAML::NodeType::Map);
  Inventory(out["tables"]["installed_sha256"], Utils::tables_dir,
            Get(fEffective, "output.directory").as<std::string>());
  const auto extra = Get(fEffective, "backgrounds.additional.file");
  if (!extra.IsNull())
    out["tables"]["additional_sha256"] = HashFile(extra.as<std::string>());
  if (!nuclearDirectory.empty()) {
    out["tables"]["nuclear_lookup"] =
        std::find(fOriginal.begin(), fOriginal.end(), "--talys-data") !=
                    fOriginal.end() ||
                std::find(fOriginal.begin(), fOriginal.end(), "-talys-data") !=
                    fOriginal.end()
            ? "legacy_cli_override"
            : "standard_tables";
    out["tables"]["nuclear_directory"] =
        fs::absolute(nuclearDirectory).lexically_normal().string();
    out["tables"]["nuclear_sha256"] = YAML::Node(YAML::NodeType::Map);
    Inventory(out["tables"]["nuclear_sha256"], nuclearDirectory,
              Get(fEffective, "output.directory").as<std::string>());
    for (const auto *name :
         {"manifest.json", "decay_manifest.json", "photon_manifest.json"}) {
      const auto path = fs::path(nuclearDirectory) / name;
      if (fs::exists(path))
        out["tables"]["nuclear_manifests"][name] =
            YAML::LoadFile(path.string());
    }
    // Preserve explicit source-data exceptions without copying process
    // settings.
    for (const auto *name : {"terminals.tsv", "pd_photon_exceptions.tsv"}) {
      const auto path = fs::path(nuclearDirectory) / name;
      if (fs::exists(path)) {
        std::ifstream file(path);
        std::ostringstream text;
        text << file.rdbuf();
        if (!file || file.bad())
          throw std::runtime_error("Cannot read " + path.string());
        out["tables"]["nuclear_exceptions"][name] = text.str();
      }
    }
  }
  for (const auto *name :
       {"OMP_THREAD_LIMIT", "OMP_NUM_THREADS", "GSL_RNG_SEED"})
    if (const char *value = std::getenv(name))
      out["environment"][name] = value;
  return out;
}
RunRecord::RunRecord(const std::string &directory, YAML::Node metadata)
    : fDirectory(directory), fMetadata(metadata) {
  try {
    fMetadata["schema_version"] = 1;
    fMetadata["status"] = "running";
    fMetadata["started_utc"] = Now();
    fMetadata["config_sha256"] =
        Configuration::HashFile((fs::path(directory) / "config.yaml").string());
    Write();
  } catch (const std::exception &e) {
    Utils::Exception::Throw(e.what());
  }
}
void RunRecord::Write() { Emit(fs::path(fDirectory) / "run.yaml", fMetadata); }
void RunRecord::Update(YAML::Node metadata) {
  try {
    metadata["config_sha256"] = fMetadata["config_sha256"];
    fMetadata = metadata;
    Write();
  } catch (const std::exception &e) {
    Utils::Exception::Throw(e.what());
  }
}
void RunRecord::Finish(bool interrupted, unsigned long failures) {
  try {
    fMetadata["execution"]["discarded_cascades"] = failures;
    fMetadata["status"] = failures      ? "failed"
                          : interrupted ? "interrupted"
                                        : "completed";
    fMetadata["finished_utc"] = Now();
    Write();
    fFinished = true;
  } catch (const std::exception &e) {
    Utils::Exception::Throw(e.what());
  }
}
RunRecord::~RunRecord() {
  if (!fFinished)
    try {
      fMetadata["status"] = "failed";
      fMetadata["finished_utc"] = Now();
      Write();
    } catch (...) {
      std::cerr << "Cannot update failed run.yaml\n";
    }
}
