// SPDX-License-Identifier: MIT
#include "GridMagneticField.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <locale>
#include <sstream>
#include <utility>

namespace Interactions {
namespace {
[[noreturn]] void GridError(const std::string& message) {
    // Input and field-contract errors cannot be repaired by retrying a history.
    Exception::ThrowFatal("Grid magnetic field: " + message);
    std::terminate();
}

void ValidateQuery(const coord_type* x, const std::vector<double>& value) {
    if (!x || value.size() != 3) GridError("query requires three coordinates and three output components");
    for (size_t a = 0; a < 3; ++a)
        if (!std::isfinite(x[a])) GridError("nonfinite query coordinate");
}

struct Row {
    std::array<coord_type, 3> position;
    std::array<double, 3> field;
};
} // namespace

SavedMF::SavedMF(const std::string& filename, FillingType filling) : filling_(filling) {
    std::ifstream input(filename);
    if (!input) GridError("cannot open '" + filename + "'");
    Load(input,filename);
}

SavedMF::SavedMF(std::istream& input, const std::string& source, FillingType filling) : filling_(filling) {
    Load(input,source);
}

void SavedMF::Load(std::istream& input, const std::string& filename) {
    if (filling_ < FILLING_ZERO || filling_ > FILLING_MIRROR) GridError("invalid boundary mode");
    input.imbue(std::locale::classic());
    std::vector<Row> rows;
    size_t lineNumber = 0;
    for (std::string line; std::getline(input, line);) {
        ++lineNumber;
        line.erase(line.find('#') == std::string::npos ? line.size() : line.find('#'));
        std::istringstream stream(line);
        stream.imbue(std::locale::classic());
        stream >> std::ws;
        if (stream.eof()) continue;
        Row row;
        if (!(stream >> row.position[0] >> row.position[1] >> row.position[2]
                     >> row.field[0] >> row.field[1] >> row.field[2]))
            GridError("'" + filename + "' line " + std::to_string(lineNumber) + ": expected six finite numbers");
        stream >> std::ws;
        if (!stream.eof()) GridError("'" + filename + "' line " + std::to_string(lineNumber) + ": extra column");
        for (size_t a = 0; a < 3; ++a) {
            if (!std::isfinite(row.position[a]) || !std::isfinite(row.field[a]))
                GridError("'" + filename + "' line " + std::to_string(lineNumber) + ": nonfinite value");
            axes_[a].push_back(row.position[a]);
        }
        rows.push_back(row);
    }
    if (input.bad()) GridError("read failure in '" + filename + "'");
    size_t count = 1;
    for (auto& axis : axes_) {
        std::sort(axis.begin(), axis.end());
        axis.erase(std::unique(axis.begin(), axis.end()), axis.end());
        if (axis.size() < 2) GridError("each axis requires at least two distinct coordinates");
        const coord_type span = axis.back() - axis.front();
        const coord_type step = span / (axis.size() - 1);
        // Permit decimal input rounding, but do not infer missing nodes or snap rows.
        const coord_type magnitude = std::max(std::abs(axis.front()), std::abs(axis.back()));
        const coord_type tolerance = 1e-10L * step +
            64 * std::numeric_limits<coord_type>::epsilon() * magnitude;
        if (!std::isfinite(span) || !std::isfinite(tolerance) || !(step > 0))
            GridError("axis extent is not representable");
        for (size_t i = 1; i < axis.size(); ++i) {
            const coord_type spacing = axis[i] - axis[i - 1];
            if (std::abs(spacing - step) > tolerance) GridError("axis is not regularly spaced");
        }
        if (count > rows.size() / axis.size()) GridError("incomplete Cartesian grid");
        count *= axis.size();
    }
    if (count != rows.size()) GridError("duplicate nodes or incomplete Cartesian grid");
    values_.resize(count);
    std::vector<bool> occupied(count, false);
    for (const Row& row : rows) {
        std::array<size_t, 3> index;
        for (size_t a = 0; a < 3; ++a)
            index[a] = std::lower_bound(axes_[a].begin(), axes_[a].end(), row.position[a]) - axes_[a].begin();
        const size_t flat = Index(index[0], index[1], index[2]);
        if (occupied[flat]) GridError("duplicate grid node");
        occupied[flat] = true;
        values_[flat] = row.field;
    }
    // Interpolate in Particle::X units. Keep a separate double-precision box
    // so that seam classification and periods use the fork's arithmetic.
    double minimumScale = std::numeric_limits<double>::max();
    for (size_t a = 0; a < axes_.size(); ++a) {
        auto& axis = axes_[a];
        bounds_[a] = {static_cast<double>(axis.front()) * units.Mpc,
                      static_cast<double>(axis.back()) * units.Mpc};
        const double legacyWidth = bounds_[a][1] - bounds_[a][0];
        if (!std::isfinite(legacyWidth) || !(legacyWidth > 0))
            GridError("box extent is not representable at boundary precision");
        // Match the fork's operation order: convert each bound to double
        // internal coordinates, subtract, then divide by node count. Rounding
        // a wider Mpc extent only at the end can change an integer step count.
        // Keep the true minimum; the fork's strict comparisons mishandle ties.
        minimumScale = std::min(minimumScale, legacyWidth / axis.size());
        for (auto& node : axis) {
            node *= units.Mpc;
            if (!std::isfinite(node)) GridError("coordinate is not representable in internal units");
        }
        for (size_t i = 1; i < axis.size(); ++i)
            if (!(axis[i] > axis[i - 1])) GridError("nodes collapse in internal units");
        const coord_type width = axis.back() - axis.front();
        if (!std::isfinite(width) ||
            (filling_ == FILLING_MIRROR && !std::isfinite(2 * legacyWidth)))
            GridError("boundary period is not representable in internal units");
    }
    minimumScale_ = minimumScale;
    if (!std::isfinite(minimumScale_) || !(minimumScale_ > 0))
        GridError("grid spacing is not representable in internal units");
}

size_t SavedMF::Index(size_t i, size_t j, size_t k) const {
    return k + axes_[2].size() * (j + axes_[1].size() * i);
}

coord_type SavedMF::MapCoordinate(coord_type value, size_t a) const {
    const auto& axis = axes_[a];
    const coord_type lo = bounds_[a][0], hi = bounds_[a][1];
    const double width = bounds_[a][1] - bounds_[a][0];
    const double period = filling_ == FILLING_MIRROR ? 2 * width : width;
    // Use the fork's subtraction and fmod order for representable differences.
    // Splitting the remainder is only needed for extreme finite coordinates.
    const coord_type difference = value - lo;
    coord_type offset;
    if (std::isfinite(difference)) {
        offset = std::fmod(difference, coord_type(period));
    } else {
        offset = std::fmod(value, coord_type(period)) -
                 std::fmod(lo, coord_type(period));
        offset = std::fmod(offset, coord_type(period));
    }
    if (offset < 0) offset += period;
    if (filling_ == FILLING_MIRROR && offset > width) offset = period - offset;
    coord_type mapped = std::clamp(lo + offset, lo, hi);
    // Match the fork's finite seam displacement. All axes are mapped once
    // any coordinate touches/leaves the box; this is part of its field model.
    const double inset = width / axis.size();
    if (mapped >= hi) mapped -= inset;
    else if (mapped <= lo) mapped += inset;
    return mapped;
}

SavedMF::Cell SavedMF::Locate(coord_type value, size_t a) const {
    const auto& axis = axes_[a];
    value = std::clamp(value, axis.front(), axis.back());
    auto right = std::upper_bound(axis.begin(), axis.end(), value);
    const size_t left = right == axis.end() ? axis.size() - 2 :
        static_cast<size_t>(right - axis.begin() - 1);
    return {left, std::clamp((value - axis[left]) / (axis[left + 1] - axis[left]),
                             coord_type(0), coord_type(1))};
}

void SavedMF::GetValueGauss(const coord_type* x, const CosmoTime&,
                            std::vector<double>& value) const {
    ValidateQuery(x, value);
    bool outside = false;
    for (size_t a = 0; a < 3; ++a)
        outside = outside || x[a] <= bounds_[a][0] || x[a] >= bounds_[a][1];
    std::fill(value.begin(), value.end(), 0.);
    if (outside && filling_ == FILLING_ZERO) return;
    if (outside && filling_ == FILLING_NEAREST) {
        // Outside means nearest Cartesian node on every axis, including axes
        // still inside the box. Clamping and interpolating changes this model.
        std::array<size_t, 3> node;
        for (size_t a = 0; a < 3; ++a) {
            const auto& axis = axes_[a];
            // Clamp before narrowing: distant finite coordinates may exceed
            // the double range even though their nearest endpoint is known.
            if (x[a] <= bounds_[a][0]) node[a] = 0;
            else if (x[a] >= bounds_[a][1]) node[a] = axis.size() - 1;
            else {
                const double coordinate = static_cast<double>(x[a]);
                const double index = std::round((axis.size() - 1) *
                    (coordinate - bounds_[a][0]) / (bounds_[a][1] - bounds_[a][0]));
                if (!std::isfinite(index)) GridError("nearest node index overflow");
                node[a] = static_cast<size_t>(index);
            }
        }
        const auto& field = values_[Index(node[0], node[1], node[2])];
        std::copy(field.begin(), field.end(), value.begin());
        return;
    }
    std::array<Cell, 3> cells;
    for (size_t a = 0; a < 3; ++a)
        cells[a] = Locate(outside ? MapCoordinate(x[a], a) : x[a], a);
    std::array<coord_type, 3> sum{};
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 2; ++j)
            for (size_t k = 0; k < 2; ++k) {
                const coord_type weight = (i ? cells[0].fraction : 1 - cells[0].fraction) *
                    (j ? cells[1].fraction : 1 - cells[1].fraction) *
                    (k ? cells[2].fraction : 1 - cells[2].fraction);
                const auto& field = values_[Index(cells[0].index + i, cells[1].index + j, cells[2].index + k)];
                for (size_t a = 0; a < 3; ++a) sum[a] += weight * field[a];
            }
    for (size_t a = 0; a < 3; ++a) {
        value[a] = static_cast<double>(sum[a]);
        if (!std::isfinite(value[a])) GridError("interpolated component overflow");
    }
}

double SavedMF::MinVariabilityScale(const CosmoTime&) const { return minimumScale_; }

ComboMF::ComboMF(std::vector<SmartPtr<MagneticField>> fields) : fields_(std::move(fields)) {
    if (fields_.empty()) GridError("composition requires at least one field");
    for (const auto& field : fields_)
        if (field == nullptr) GridError("null component in composition");
}

void ComboMF::GetValueGauss(const coord_type* x, const CosmoTime& time,
                           std::vector<double>& value) const {
    ValidateQuery(x, value);
    std::array<coord_type, 3> sum{};
    std::vector<double> component(3);
    for (const auto& field : fields_) {
        std::fill(component.begin(), component.end(), 0.);
        field->GetValueGauss(x, time, component);
        if (component.size() != 3) GridError("component returned an invalid vector size");
        for (size_t a = 0; a < 3; ++a) {
            if (!std::isfinite(component[a])) GridError("component returned a nonfinite field");
            sum[a] += component[a];
        }
    }
    for (size_t a = 0; a < 3; ++a) {
        value[a] = static_cast<double>(sum[a]);
        if (!std::isfinite(value[a])) GridError("field sum overflow");
    }
}

double ComboMF::MinVariabilityScale(const CosmoTime& time) const {
    double scale = std::numeric_limits<double>::max();
    for (const auto& field : fields_) {
        const double component = field->MinVariabilityScale(time);
        // The existing deflector divides by this scale, including in zero-field regions.
        if (!std::isfinite(component) || !(component > 0)) GridError("component scale must be finite and positive");
        scale = std::min(scale, component);
    }
    return scale;
}
} // namespace Interactions
