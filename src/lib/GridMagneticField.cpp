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
    if (filling < FILLING_ZERO || filling > FILLING_MIRROR) GridError("invalid boundary mode");
    std::ifstream input(filename);
    if (!input) GridError("cannot open '" + filename + "'");
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
    coord_type minimumSpacing = std::numeric_limits<coord_type>::max();
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
            minimumSpacing = std::min(minimumSpacing, spacing);
        }
        if (count > rows.size() / axis.size()) GridError("incomplete Cartesian grid");
        count *= axis.size();
    }
    if (count != rows.size()) GridError("duplicate nodes or incomplete Cartesian grid");
    minimumScale_ = static_cast<double>(minimumSpacing * units.Mpc);
    if (!std::isfinite(minimumScale_) || !(minimumScale_ > 0))
        GridError("grid spacing is not representable in internal units");
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
    // Store nodes in the same units as Particle::X, avoiding a divide-back
    // rounding error that could classify an exact boundary node as outside.
    for (auto& axis : axes_) {
        for (auto& node : axis) {
            node *= units.Mpc;
            if (!std::isfinite(node)) GridError("coordinate is not representable in internal units");
        }
        for (size_t i = 1; i < axis.size(); ++i)
            if (!(axis[i] > axis[i - 1])) GridError("nodes collapse in internal units");
        const coord_type width = axis.back() - axis.front();
        if (!std::isfinite(width) ||
            (filling_ == FILLING_MIRROR && !std::isfinite(2 * width)))
            GridError("boundary period is not representable in internal units");
    }
}

size_t SavedMF::Index(size_t i, size_t j, size_t k) const {
    return k + axes_[2].size() * (j + axes_[1].size() * i);
}

coord_type SavedMF::MapCoordinate(coord_type value, size_t a) const {
    const auto& axis = axes_[a];
    const coord_type lo = axis.front(), hi = axis.back();
    if (value >= lo && value <= hi) return value;
    if (filling_ == FILLING_NEAREST) return std::clamp(value, lo, hi);
    const coord_type width = hi - lo;
    const coord_type period = filling_ == FILLING_MIRROR ? 2 * width : width;
    // Separate remainders avoid overflow in value-lo for distant finite queries.
    coord_type v = std::fmod(value, period), origin = std::fmod(lo, period);
    if (v < 0) v += period;
    if (origin < 0) origin += period;
    coord_type offset = v - origin;
    if (offset < 0) offset += period;
    if (filling_ == FILLING_MIRROR && offset > width) offset = period - offset;
    return std::clamp(lo + offset, lo, hi);
}

SavedMF::Cell SavedMF::Locate(coord_type value, size_t a) const {
    const auto& axis = axes_[a];
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
        outside = outside || x[a] < axes_[a].front() || x[a] > axes_[a].back();
    std::fill(value.begin(), value.end(), 0.);
    if (outside && filling_ == FILLING_ZERO) return;
    std::array<Cell, 3> cells;
    for (size_t a = 0; a < 3; ++a) cells[a] = Locate(MapCoordinate(x[a], a), a);
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
