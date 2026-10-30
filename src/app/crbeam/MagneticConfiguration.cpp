#include "MagneticConfiguration.h"
#include "GridMagneticField.h"
#include "Sha256.h"
#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
[[noreturn]] void Invalid(const std::string &why) {
  throw std::runtime_error("magnetic_field: " + why);
}
std::string CliText(cors::cmdline::CmdLine &cli, const char *option) {
  auto parameter = cli.get_param(option);
  if (cli.has_param(option)) {
    if (!parameter.value || !*parameter.value)
      Invalid(std::string(option) + " requires a nonempty argument");
    return parameter.value;
  }
  return static_cast<const char *>(parameter);
}
std::string Text(const YAML::Node &node, const std::string &key) {
  if (!node.IsScalar())
    Invalid(key + " must be a string");
  const auto value = node.Scalar();
  std::string lower = value;
  std::transform(lower.begin(), lower.end(), lower.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (node.Tag() == "?" &&
      (lower == "true" || lower == "false" || lower == ".nan" ||
       lower == ".inf" || lower == "+.inf" || lower == "-.inf"))
    Invalid(key + " must be a string");
  if (node.Tag() == "?") {
    try {
      size_t end = 0;
      std::stod(value, &end);
      if (end == value.size())
        Invalid(key + " must be a string");
    } catch (const std::invalid_argument &) {
    } catch (const std::out_of_range &) {
    }
  }
  return value;
}
double Number(const YAML::Node &node, const std::string &key,
              bool zero = false) {
  if (!node.IsScalar() || node.Tag() == "!")
    Invalid(key + " must be numeric");
  size_t end = 0;
  double value;
  try {
    value = std::stod(node.Scalar(), &end);
  } catch (const std::exception &) {
    Invalid(key + " must be finite numeric data");
  }
  if (end != node.Scalar().size() || !std::isfinite(value) ||
      (zero ? value < 0 : value <= 0))
    Invalid(key + (zero ? " must be finite and nonnegative"
                        : " must be finite and positive"));
  return value;
}
int Integer(const YAML::Node &node, const std::string &key) {
  if (!node.IsScalar() || node.Tag() == "!")
    Invalid(key + " must be an integer");
  int value = 0;
  const auto &s = node.Scalar();
  const auto parsed = std::from_chars(s.data(), s.data() + s.size(), value);
  if (parsed.ec != std::errc() || parsed.ptr != s.data() + s.size() ||
      value <= 0)
    Invalid(key + " must be a positive supported integer");
  return value;
}
Interactions::SavedMF::FillingType Boundary(const std::string &value) {
  using Field = Interactions::SavedMF;
  if (value == "zero")
    return Field::FILLING_ZERO;
  if (value == "nearest")
    return Field::FILLING_NEAREST;
  if (value == "periodic")
    return Field::FILLING_REPEAT;
  if (value == "mirror")
    return Field::FILLING_MIRROR;
  Invalid("boundary choices: zero|nearest|periodic|mirror");
}
std::string GridKey(const YAML::Node &field) {
  return field["boundary"].as<std::string>() + ":" +
         field["file"].as<std::string>();
}
YAML::Node Normalize(const YAML::Node &requested, const YAML::Node &legacy,
                     const std::string &base, bool component) {
  if (!requested.IsMap())
    Invalid("expected a mapping");
  const auto model = requested["model"]
                         ? Text(requested["model"], "model")
                         : (component ? "" : legacy["model"].as<std::string>());
  if (model != "none" && model != "grid" && model != "sum" &&
      model != "kolmogorov" && model != "monochromatic")
    Invalid("model choices: none|monochromatic|kolmogorov|grid|sum");
  if (component && (model == "sum" || model == "none"))
    Invalid("sum components must be grid, monochromatic or kolmogorov fields");
  std::set<std::string> keys{"model"};
  if (!component) {
    keys.insert("realization");
    keys.insert("deflection_accuracy");
  }
  if (model == "grid") {
    keys.insert("file");
    keys.insert("boundary");
  }
  if (model == "sum")
    keys.insert("components");
  if (model == "monochromatic" || model == "kolmogorov") {
    keys.insert("strength_G");
    keys.insert("maximum_scale_Mpc");
    if (model == "kolmogorov") {
      keys.insert("minimum_scale_Mpc");
      keys.insert("modes");
    }
  }
  std::set<std::string> seen;
  for (const auto &entry : requested) {
    const auto key = Text(entry.first, "key");
    if (!seen.insert(key).second)
      Invalid("duplicate key " + key);
    if (!keys.count(key))
      Invalid("inapplicable or unknown key " + key + " for " + model);
  }
  YAML::Node result;
  result["model"] = model;
  auto value = [&](const char *name) {
    return requested[name] ? YAML::Clone(requested[name])
                           : YAML::Clone(legacy[name]);
  };
  if (!component) {
    result["deflection_accuracy"] =
        Integer(value("deflection_accuracy"), "deflection_accuracy");
    const auto realization = value("realization")
                                 ? Text(value("realization"), "realization")
                                 : "shared";
    if (realization != "shared" && realization != "per_primary")
      Invalid("realization choices: shared|per_primary");
    if ((model == "grid" || model == "none") && realization != "shared")
      Invalid(model + " requires shared realization");
    result["realization"] = realization;
  }
  if (model == "grid") {
    if (!requested["file"])
      Invalid("grid requires file");
    const auto file = Text(requested["file"], "file");
    if (file.empty())
      Invalid("file must not be empty");
    const std::filesystem::path path(file);
    result["file"] =
        std::filesystem::absolute(
            path.is_absolute() ? path : std::filesystem::path(base) / path)
            .lexically_normal()
            .string();
    const auto boundary = requested["boundary"]
                              ? Text(requested["boundary"], "boundary")
                              : "zero";
    Boundary(boundary);
    result["boundary"] = boundary;
  } else if (model == "sum") {
    const auto components = requested["components"];
    if (!components.IsSequence() || components.size() == 0)
      Invalid("sum requires a nonempty components sequence");
    result["components"] = YAML::Node(YAML::NodeType::Sequence);
    bool stochastic = false;
    for (const auto &child : components) {
      const auto normalized = Normalize(child, legacy, base, true);
      stochastic =
          stochastic || normalized["model"].as<std::string>() != "grid";
      result["components"].push_back(normalized);
    }
    if (!stochastic && result["realization"].as<std::string>() == "per_primary")
      Invalid("per_primary sum requires a stochastic component");
  } else if (model != "none") {
    result["strength_G"] = Number(value("strength_G"), "strength_G", true);
    result["maximum_scale_Mpc"] =
        Number(value("maximum_scale_Mpc"), "maximum_scale_Mpc");
    if (model == "kolmogorov") {
      result["minimum_scale_Mpc"] =
          Number(value("minimum_scale_Mpc"), "minimum_scale_Mpc");
      result["modes"] = Integer(value("modes"), "modes");
      if (result["modes"].as<int>() < 2)
        Invalid("kolmogorov requires at least two modes");
      if (result["minimum_scale_Mpc"].as<double>() >=
          result["maximum_scale_Mpc"].as<double>())
        Invalid("minimum scale must be below maximum");
    }
  }
  return result;
}
} // namespace

YAML::Node NormalizeMagneticField(const YAML::Node &requested,
                                  const YAML::Node &legacy,
                                  const std::string &base) {
  // If omitted entirely, preserve the legacy default and zero-strength
  // disabling.
  if (!requested) {
    YAML::Node defaults = YAML::Clone(legacy);
    if (defaults["strength_G"].as<double>() == 0) {
      YAML::Node none;
      none["model"] = "none";
      none["deflection_accuracy"] = defaults["deflection_accuracy"];
      return Normalize(none, legacy, base, false);
    }
    if (defaults["model"].as<std::string>() == "monochromatic") {
      defaults.remove("minimum_scale_Mpc");
      defaults.remove("modes");
    }
    return Normalize(defaults, legacy, base, false);
  }
  return Normalize(requested, legacy, base, false);
}

YAML::Node MagneticFieldFromCli(cors::cmdline::CmdLine &cli,
                                const YAML::Node &legacy) {
  const bool explicitModel = cli.has_param("--magnetic-field");
  if (!explicitModel && !cli.has_param("--magnetic-field-file") &&
      !cli.has_param("--magnetic-field-boundary"))
    return YAML::Node(YAML::NodeType::Undefined);
  if (!explicitModel)
    Invalid("--magnetic-field-file/boundary require --magnetic-field grid");
  const std::string model = CliText(cli, "--magnetic-field");
  YAML::Node requested;
  requested["model"] = model;
  requested["deflection_accuracy"] =
      cli.has_param("--accEGMF")
          ? Integer(YAML::Node(CliText(cli, "--accEGMF")),
                    "deflection_accuracy")
          : legacy["deflection_accuracy"].as<int>();
  requested["realization"] = legacy["realization"];
  if (model == "sum")
    Invalid("sum is available through YAML components only");
  if (cli.has_param("--magnetic-field-file"))
    requested["file"] = CliText(cli, "--magnetic-field-file");
  if (cli.has_param("--magnetic-field-boundary"))
    requested["boundary"] = CliText(cli, "--magnetic-field-boundary");
  if (model == "grid" || model == "none") {
    for (const char *option :
         {"--EGMF", "--lminEGMF", "--lmaxEGMF", "--nmodesEGMF",
          "--turbulentEGMF", "--randomEGMF"})
      if (cli.has_param(option))
        Invalid(std::string(option) + " is incompatible with " + model);
  } else if (model == "monochromatic" || model == "kolmogorov") {
    if (model == "monochromatic" && cli.has_param("--turbulentEGMF"))
      Invalid("--turbulentEGMF contradicts monochromatic model");
    if (model == "monochromatic" &&
        (cli.has_param("--lminEGMF") || cli.has_param("--nmodesEGMF")))
      Invalid("minimum scale/modes apply only to kolmogorov");
    auto numeric = [&](const char *option, const char *key, bool zero = false) {
      return cli.has_param(option)
                 ? Number(YAML::Node(CliText(cli, option)), key, zero)
                 : legacy[key].as<double>();
    };
    requested["strength_G"] = numeric("--EGMF", "strength_G", true);
    requested["maximum_scale_Mpc"] = numeric("--lmaxEGMF", "maximum_scale_Mpc");
    if (model == "kolmogorov") {
      requested["minimum_scale_Mpc"] =
          numeric("--lminEGMF", "minimum_scale_Mpc");
      requested["modes"] =
          cli.has_param("--nmodesEGMF")
              ? Integer(YAML::Node(CliText(cli, "--nmodesEGMF")), "modes")
              : legacy["modes"].as<int>();
    }
  }
  return NormalizeMagneticField(requested, legacy,
                                std::filesystem::current_path().string());
}

void MagneticLegacyArguments(const YAML::Node &field,
                             std::vector<std::string> &arguments) {
  auto add = [&](const char *key, const YAML::Node &value) {
    YAML::Emitter text;
    text.SetDoublePrecision(17);
    text << value;
    arguments.emplace_back(key);
    arguments.emplace_back(text.c_str());
  };
  add("--accEGMF", field["deflection_accuracy"]);
  const auto model = field["model"].as<std::string>();
  if (model == "grid" || model == "sum" || model == "none") {
    arguments.emplace_back("--EGMF");
    arguments.emplace_back("0");
    return;
  }
  add("--EGMF", field["strength_G"]);
  add("--lmaxEGMF", field["maximum_scale_Mpc"]);
  if (model == "kolmogorov") {
    arguments.emplace_back("--turbulentEGMF");
    add("--lminEGMF", field["minimum_scale_Mpc"]);
    add("--nmodesEGMF", field["modes"]);
  }
  if (field["realization"].as<std::string>() == "per_primary")
    arguments.emplace_back("--randomEGMF");
}

ConfiguredMagneticField::ConfiguredMagneticField(const YAML::Node &field)
    : field_(YAML::Clone(field)) {
  provenance_ = YAML::Node(YAML::NodeType::Map);
  LoadGrids(field_);
}
void ConfiguredMagneticField::LoadGrids(const YAML::Node &field) {
  const auto model = field["model"].as<std::string>();
  if (model == "sum") {
    for (const auto &child : field["components"])
      LoadGrids(child);
  }
  if (model != "grid")
    return;
  const auto file = field["file"].as<std::string>();
  const auto key = GridKey(field);
  if (grids_.count(key))
    return;
  std::ifstream input(file, std::ios::binary);
  if (!input)
    Invalid("Cannot hash/read grid file: " + file);
  const std::string bytes{std::istreambuf_iterator<char>(input), {}};
  if (input.bad())
    Invalid("Cannot read grid file: " + file);
  const auto hash = Sha256Bytes(bytes);
  const auto previous = static_cast<const YAML::Node &>(provenance_)[file];
  if (previous && previous.as<std::string>() != hash)
    Invalid("grid changed between component loads: " + file);
  std::istringstream frozen(bytes);
  grids_[key] = new Interactions::SavedMF(
      frozen, file, Boundary(field["boundary"].as<std::string>()));
  provenance_[file] = hash;
}
Utils::SmartPtr<Interactions::MagneticField>
ConfiguredMagneticField::Build(mcray::Randomizer &random) const {
  return BuildNode(field_, random);
}
Utils::SmartPtr<Interactions::MagneticField>
ConfiguredMagneticField::BuildNode(const YAML::Node &field,
                                   mcray::Randomizer &random) const {
  using namespace Interactions;
  const auto model = field["model"].as<std::string>();
  if (model == "grid")
    return grids_.at(GridKey(field));
  if (model == "sum") {
    std::vector<Utils::SmartPtr<MagneticField>> children;
    for (const auto &child : field["components"])
      children.push_back(BuildNode(child, random));
    return new ComboMF(std::move(children));
  }
  const auto strength = field["strength_G"].as<double>();
  const auto scale = field["maximum_scale_Mpc"].as<double>();
  if (model == "kolmogorov")
    return new TurbulentMF(random, field["minimum_scale_Mpc"].as<double>(),
                           scale, strength, field["modes"].as<int>());
  if (field_["realization"].as<std::string>() == "per_primary")
    return new MonochromaticMF(random, scale, strength);
  return new MonochromaticMF(scale, strength, 0.1, M_PI / 4, 0, 0);
}
