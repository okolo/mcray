#include "SpatialOutput.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

namespace mcray {
namespace {
coord_type Flight(const Particle &p) {
  const coord_type flight = p.Time.t() - p.SourceParticle->Time.t();
  if (!std::isfinite(flight) || flight < 0)
    Exception::ThrowFatal("Invalid primary flight time in spatial output");
  return flight;
}
void Precision(std::ostream &out) {
  out << std::defaultfloat
      << std::setprecision(std::numeric_limits<coord_type>::max_digits10);
}
} // namespace

SpatialOutput::SpatialOutput(std::string directory, bool overwrite,
                             bool saveSourceDirection, bool saveId,
                             std::shared_ptr<const SpatialObserver> observer)
    : RawOutput(std::move(directory), false, overwrite),
      saveSourceDirection_(saveSourceDirection), saveId_(saveId),
      observer_(std::move(observer)) {}

void SpatialOutput::WriteHeader(std::ostream &out, const Particle &first) {
  const auto flight = Flight(first);
  coord_type scale = flight;
  if (flight == 0) {
    scale = units.Mpc;
    for (const auto x : first.X)
      scale = std::max(scale, std::abs(x));
  }
  if (!std::isfinite(scale) || scale <= 0)
    Exception::ThrowFatal("Invalid spatial output normalization scale");
  scales_[first.Name()] = scale;
  Precision(out);
  out << "# " << first.Name() << "\tT/Mpc = ";
  // Old readers accept only digits and a decimal point in this field.
  // Preserve significant digits even when the scale is smaller than one Mpc.
  const auto value = scale / units.Mpc;
  if (!std::isfinite(value) || value <= 0)
    Exception::ThrowFatal("Spatial output scale is not representable in Mpc");
  const int decimals =
      std::max(0, std::numeric_limits<coord_type>::max_digits10 - 1 -
                      static_cast<int>(std::floor(std::log10(value))));
  out << std::fixed << std::setprecision(decimals) << value;
  Precision(out);
  out << "\tnormalization = "
      << (flight > 0 ? "first_particle_flight" : "fallback")
      << "\tspatial_format = 1";
  for (size_t i = 0; i < 3; ++i)
    out << "\tsource_" << "xyz"[i]
        << "_Mpc = " << first.SourceParticle->X[i] / units.Mpc;
  out << "\n#\n# "
         "E/eV\tweight\tx/T\ty/T\t1-(z/T)\tPx/P\tPy/P\t1-(Pz/P)\tz_production";
  if (observer_)
    out << "\tin_target";
  if (saveSourceDirection_)
    out << "\tsource_x\tsource_y\tsource_z";
  out << "\tz_src\tE_src/eV";
  if (saveId_)
    out << "\tId";
  out << "\tPathMpc\tFlightMpc";
}

void SpatialOutput::Write(std::ostream &out, const Particle &p) {
  const auto scale = scales_.at(p.Name());
  coord_type coordinates[3];
  for (size_t i = 0; i < 3; ++i) {
    coordinates[i] = p.X[i] / scale;
    if (!std::isfinite(coordinates[i]))
      Exception::ThrowFatal(
          "Spatial coordinate cannot be normalized with the file scale");
  }
  Precision(out);
  out << p.Energy / units.eV << '\t' << p.Weight << '\t' << coordinates[0]
      << '\t' << coordinates[1] << '\t' << 1 - coordinates[2] << '\t'
      << p.Pdir[0] << '\t' << p.Pdir[1] << '\t' << 1 - p.Pdir[2] << '\t'
      << p.fProductionTime.z();
  if (observer_)
    out << '\t' << observer_->Reached(p);
  if (saveSourceDirection_)
    for (const auto direction : p.SourceParticle->Pdir)
      out << '\t' << direction;
  out << '\t' << p.SourceParticle->Time.z() << '\t'
      << p.SourceParticle->Energy / units.eV;
  if (saveId_)
    out << '\t' << p.id;
  out << '\t' << (p.Time - p.fProductionTime) / units.Mpc << '\t'
      << Flight(p) / units.Mpc;
}
} // namespace mcray
