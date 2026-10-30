// Cartesian magnetic fields and immutable field composition.
// SPDX-License-Identifier: MIT
#ifndef MCRAY_GRID_MAGNETIC_FIELD_H
#define MCRAY_GRID_MAGNETIC_FIELD_H

#include "Deflection3D.h"
#include <array>
#include <string>
#include <vector>

namespace Interactions {

// Static field: input coordinates are comoving Mpc, components are Gauss.
// No redshift scaling or component sign change is applied by this class.
class SavedMF final : public MagneticField {
public:
    enum FillingType : int { FILLING_ZERO, FILLING_NEAREST, FILLING_REPEAT, FILLING_MIRROR };
    explicit SavedMF(const std::string& filename, FillingType filling = FILLING_ZERO);
    SavedMF(const SavedMF&) = delete;
    SavedMF& operator=(const SavedMF&) = delete;

    void GetValueGauss(const coord_type* x, const CosmoTime& time,
                       std::vector<double>& value) const override;
    double MinVariabilityScale(const CosmoTime& time) const override;

private:
    using Axis = std::vector<coord_type>;
    struct Cell { size_t index; coord_type fraction; };
    coord_type MapCoordinate(coord_type value, size_t axis) const;
    Cell Locate(coord_type value, size_t axis) const;
    size_t Index(size_t i, size_t j, size_t k) const;

    std::array<Axis, 3> axes_;
    std::vector<std::array<double, 3>> values_;
    FillingType filling_;
    double minimumScale_;
};

// Components are retained through the library's intrusive smart pointers.
// Assemble before worker threads start; component field queries must be const-safe.
class ComboMF final : public MagneticField {
public:
    explicit ComboMF(std::vector<SmartPtr<MagneticField>> fields);
    ComboMF(const ComboMF&) = delete;
    ComboMF& operator=(const ComboMF&) = delete;

    void GetValueGauss(const coord_type* x, const CosmoTime& time,
                       std::vector<double>& value) const override;
    double MinVariabilityScale(const CosmoTime& time) const override;

private:
    std::vector<SmartPtr<MagneticField>> fields_;
};

} // namespace Interactions
#endif
