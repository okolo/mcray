#include "GridMagneticField.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <unistd.h>

namespace {
using namespace mcray;
using namespace Interactions;
using Utils::SmartPtr;

class GridFile {
public:
    explicit GridFile(const std::string& contents) {
        static std::atomic<unsigned> serial{0};
        path = std::filesystem::temp_directory_path() /
            ("mcray-grid-" + std::to_string(getpid()) + "-" + std::to_string(serial++));
        std::ofstream out(path); out << contents;
        if (!out) throw std::runtime_error("Cannot write test grid");
    }
    ~GridFile() { std::filesystem::remove(path); }
    std::filesystem::path path;
};

std::array<double, 3> Affine(coord_type x, coord_type y, coord_type z) {
    return {double(1 + 2*x - 3*y + z), double(4 - x + y + 2*z), double(-2 + x - y - z)};
}

std::string Grid(coord_type shift = 0, bool constant = false) {
    std::ostringstream out;
    out << std::setprecision(21) << "# x/Mpc y/Mpc z/Mpc Bx/G By/G Bz/G\n\n";
    // Reverse row order and use a translated rectangular grid with unequal spacing.
    for (int i = 2; i >= 0; --i)
        for (int j = 1; j >= 0; --j)
            for (int k = 2; k >= 0; --k) {
                const coord_type x = shift + i, y = shift + 2*j, z = shift + 3*k;
                const auto b = constant ? std::array<double,3>{0, 1e-9, 0} : Affine(x,y,z);
                out << x << ' ' << y << ' ' << z << ' ' << b[0] << ' ' << b[1] << ' ' << b[2] << " # node\n";
            }
    return out.str();
}

std::array<double, 3> Query(const MagneticField& field, coord_type x, coord_type y, coord_type z,
                           double redshift = 0) {
    const coord_type position[] = {x*units.Mpc, y*units.Mpc, z*units.Mpc};
    std::vector<double> b(3, -999);
    field.GetValueGauss(position, CosmoTime(redshift), b);
    return {b[0], b[1], b[2]};
}

void ExpectField(const MagneticField& field, coord_type x, coord_type y, coord_type z,
                 coord_type mx, coord_type my, coord_type mz) {
    const auto actual = Query(field,x,y,z), expected = Affine(mx,my,mz);
    for (size_t a=0; a<3; ++a) EXPECT_NEAR(actual[a],expected[a],2e-12) << "component " << a;
}

void ExpectFatal(const std::function<void()>& action, const std::string& text) {
    try { action(); FAIL() << "Expected fatal field error"; }
    catch (Utils::Exception* error) {
        std::unique_ptr<Utils::Exception> owner(error);
        EXPECT_FALSE(error->Retryable());
        EXPECT_NE(error->Message().find(text), std::string::npos) << error->Message();
    }
}

class GridFields : public ::testing::Test {
    void SetUp() override { if (!cosmology.IsInitialized()) cosmology.Init(); }
};

TEST_F(GridFields, AffineInteriorAndPositiveHistoricalScale) {
    GridFile input(Grid(10)); SavedMF field(input.path.string());
    for (coord_type x : {10.25L,11.L,11.9L})
        for (coord_type y : {10.2L,10.7L,11.9L})
            for (coord_type z : {10.1L,11.3L,13.L,15.9L})
                ExpectField(field,x,y,z,x,y,z);
    EXPECT_DOUBLE_EQ(field.MinVariabilityScale(CosmoTime(0)),2*units.Mpc/3);
    EXPECT_EQ(Query(field,11,11,11,0), Query(field,11,11,11,1));
    EXPECT_EQ(Query(field,9.999L,11,11),(std::array<double,3>{0,0,0}));
    EXPECT_EQ(Query(field,11,12.001L,11),(std::array<double,3>{0,0,0}));
    EXPECT_EQ(Query(field,11,11,16.001L),(std::array<double,3>{0,0,0}));
}

TEST_F(GridFields, NegativeAndFractionalDomainsUseTheirActualExtrema) {
    // Preserve the loader correction: the fork incorrectly initializes extrema
    // to zero, leaving uninitialized nodes for boxes not containing the origin.
    GridFile negative(Grid(-10)); SavedMF field(negative.path.string());
    ExpectField(field,-9.5L,-9.5L,-9,-9.5L,-9.5L,-9);
    GridFile fractional(Grid(0.3L)); SavedMF fractions(fractional.path.string());
    ExpectField(fractions,0.8L,0.8L,1.3L,0.8L,0.8L,1.3L);
}

TEST_F(GridFields, ZeroModeTreatsEveryFaceAsOutside) {
    GridFile input(Grid()); SavedMF field(input.path.string());
    const double maximum[] = {2,2,6};
    for (size_t axis=0;axis<3;++axis) for (bool upper : {false,true}) {
        coord_type position[] = {units.Mpc,units.Mpc,units.Mpc};
        // Boundary classification and periods retain the fork's double box.
        position[axis] = upper ? maximum[axis]*units.Mpc : 0;
        std::vector<double> value(3,99);
        field.GetValueGauss(position,CosmoTime(0),value);
        EXPECT_EQ(value,(std::vector<double>{0,0,0}));
    }
}

TEST_F(GridFields, NearestModeRoundsAllAxesWhenAnyAxisIsOutside) {
    GridFile input(Grid()); SavedMF field(input.path.string(),SavedMF::FILLING_NEAREST);
    ExpectField(field,-1,0.7L,1.3L,0,0,0);
    ExpectField(field,0,1.3L,2.1L,0,2,3);
    ExpectField(field,1.2L,0.7L,1.3L,1.2L,0.7L,1.3L); // Interior stays trilinear.
}

TEST_F(GridFields, NearestModeClampsBeforeNarrowingExtremeCoordinates) {
    GridFile input(Grid()); SavedMF field(input.path.string(),SavedMF::FILLING_NEAREST);
    // On Linux, coord_type can hold finite values outside the double range.
    // Query directly in internal units so conversion from Mpc cannot overflow.
    const coord_type far = std::numeric_limits<coord_type>::max()/2;
    for (bool positive : {false,true}) {
        const coord_type position[] = {positive ? far : -far, 0.7L*units.Mpc, 1.3L*units.Mpc};
        std::vector<double> actual(3);
        field.GetValueGauss(position,CosmoTime(0),actual);
        const auto expected = Affine(positive ? 2 : 0,0,0);
        for (size_t a=0;a<3;++a) EXPECT_DOUBLE_EQ(actual[a],expected[a]);
    }
}

TEST_F(GridFields, RepeatAndMirrorPreserveHistoricalSeamInsets) {
    GridFile input(Grid());
    SavedMF repeat(input.path.string(),SavedMF::FILLING_REPEAT);
    SavedMF mirror(input.path.string(),SavedMF::FILLING_MIRROR);
    ExpectField(repeat,-2,0,0,2.L/3,1,2);
    ExpectField(mirror,-2,0,0,4.L/3,1,2);
    // A box with translated limits exercises folding without seam coincidences.
    GridFile translated(Grid(10));
    SavedMF shiftedRepeat(translated.path.string(),SavedMF::FILLING_REPEAT);
    SavedMF shiftedMirror(translated.path.string(),SavedMF::FILLING_MIRROR);
    ExpectField(shiftedRepeat,8.25L,14.5L,22.5L,10.25L,10.5L,10.5L);
    ExpectField(shiftedMirror,12.25L,13.5L,16.5L,11.75L,10.5L,15.5L);
    ExpectField(shiftedMirror,9.75L,9.5L,9.5L,10.25L,10.5L,10.5L);
    const auto far = Query(repeat,1e200L,-1e200L,1e200L);
    for (double b : far) EXPECT_TRUE(std::isfinite(b));
}

TEST_F(GridFields, ScaleTiesTakeTheTrueMinimum) {
    // x and y both have extent/count=1; z has extent/count=2.
    // The fork's strict a<b && a<c comparisons incorrectly return 2 here.
    std::ostringstream data;
    for (int i=0;i<2;++i) for (int j=0;j<2;++j) for (int k=0;k<2;++k)
        data << 2*i << ' ' << 2*j << ' ' << 4*k << " 0 1 0\n";
    GridFile input(data.str()); SavedMF field(input.path.string());
    EXPECT_DOUBLE_EQ(field.MinVariabilityScale(CosmoTime(0)),units.Mpc);
}

TEST_F(GridFields, RejectsMalformedAndIncompleteInput) {
    GridFile valid(Grid());
    ExpectFatal([&] { SavedMF f(valid.path.string()+"-missing"); },"cannot open");
    ExpectFatal([&] { SavedMF f(valid.path.string(), static_cast<SavedMF::FillingType>(99)); },"boundary mode");
    for (const char* contents : {"", "# empty\n", "0 0 0 1 2\n", "0 0 0 1 2 3 extra\n",
                                      "0 0 0 nan 2 3\n", "inf 0 0 1 2 3\n", "0 0 0 1e999 2 3\n"}) {
        GridFile bad(contents);
        ExpectFatal([&] { SavedMF f(bad.path.string()); },"Grid magnetic field");
    }
    std::string one = Grid();
    one.erase(one.find_last_of('\n',one.size()-2)+1);
    GridFile missing(one);
    ExpectFatal([&] { SavedMF f(missing.path.string()); },"grid");
    GridFile duplicate(Grid()+"0 0 0 1 2 3\n");
    ExpectFatal([&] { SavedMF f(duplicate.path.string()); },"duplicate");
    std::string compensated = one + "2 2 6 1 2 3\n"; // Missing a node plus duplicate, unchanged row count.
    GridFile both(compensated);
    ExpectFatal([&] { SavedMF f(both.path.string()); },"duplicate");
    std::ostringstream nonuniform, flat;
    for (int i=0;i<3;++i) for (int j=0;j<2;++j) for (int k=0;k<2;++k)
        nonuniform << (i==2?3:i) << ' ' << j << ' ' << k << " 1 2 3\n";
    for (int j=0;j<2;++j) for (int k=0;k<2;++k) flat << "0 " << j << ' ' << k << " 1 2 3\n";
    GridFile irregular(nonuniform.str()), plane(flat.str());
    ExpectFatal([&] { SavedMF f(irregular.path.string()); },"regularly spaced");
    ExpectFatal([&] { SavedMF f(plane.path.string()); },"at least two");
}

TEST_F(GridFields, RejectsInvalidQueryWithoutResampling) {
    GridFile input(Grid()); SavedMF field(input.path.string());
    coord_type position[3] = {0,0,std::numeric_limits<coord_type>::infinity()};
    std::vector<double> b(3);
    ExpectFatal([&] { field.GetValueGauss(position,CosmoTime(0),b); },"nonfinite");
    position[2] = 0; b.resize(2);
    ExpectFatal([&] { field.GetValueGauss(position,CosmoTime(0),b); },"three output");
    b.resize(3);
    ExpectFatal([&] { field.GetValueGauss(nullptr,CosmoTime(0),b); },"three coordinates");
}

class ConstantField : public MagneticField {
public:
    ConstantField(std::array<double,3> b, double scale, int* deaths = nullptr)
        : b_(b), scale_(scale), deaths_(deaths) {}
    ~ConstantField() override { if (deaths_) ++*deaths_; }
    void GetValueGauss(const coord_type*,const CosmoTime&,std::vector<double>& out) const override {
        out.assign(b_.begin(),b_.end());
    }
    double MinVariabilityScale(const CosmoTime&) const override { return scale_; }
private:
    std::array<double,3> b_;
    double scale_;
    int* deaths_;
};

TEST_F(GridFields, CompositionRetainsComponentsAndUsesSmallestPositiveScale) {
    int deaths = 0;
    SmartPtr<MagneticField> combo;
    {
        SmartPtr<MagneticField> a = new ConstantField({1,2,3},5*units.Mpc,&deaths);
        SmartPtr<MagneticField> b = new ConstantField({-1,4,-3},2*units.Mpc,&deaths);
        combo = new ComboMF({a,b});
    }
    EXPECT_EQ(deaths,0);
    EXPECT_EQ(Query(*combo,0,0,0),(std::array<double,3>{0,6,0}));
    EXPECT_DOUBLE_EQ(combo->MinVariabilityScale(CosmoTime(0)),2*units.Mpc);
    std::unique_ptr<DeflectionInteraction> clone;
    {
        Deflection3D original(combo);
        clone.reset(original.Clone());
    }
    combo = nullptr;
    EXPECT_EQ(deaths,0); // The clone still owns the composed field and its children.
    Particle electron(Electron,0.01); electron.Energy=1e20*units.eV;
    EXPECT_TRUE(std::isfinite(clone->Rate(electron)));
    clone.reset();
    EXPECT_EQ(deaths,2);
    ExpectFatal([&] { ComboMF empty({}); },"at least one");
    ExpectFatal([&] { ComboMF null({SmartPtr<MagneticField>()}); },"null");
}

TEST_F(GridFields, CompositionRejectsInvalidScaleAndOverflow) {
    for (double scale : {0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        SmartPtr<MagneticField> child = new ConstantField({0,0,0},scale);
        ComboMF combo({child});
        ExpectFatal([&] { combo.MinVariabilityScale(CosmoTime(0)); },"finite and positive");
    }
    SmartPtr<MagneticField> huge = new ConstantField({std::numeric_limits<double>::max(),0,0},1);
    ComboMF overflow({huge,huge});
    ExpectFatal([&] { Query(overflow,0,0,0); },"overflow");
    SmartPtr<MagneticField> nan = new ConstantField({std::numeric_limits<double>::quiet_NaN(),0,0},1);
    ComboMF invalid({nan});
    ExpectFatal([&] { Query(invalid,0,0,0); },"nonfinite");
}

TEST_F(GridFields, ConcurrentQueriesShareImmutableField) {
    GridFile input(Grid());
    SmartPtr<MagneticField> grid = new SavedMF(input.path.string());
    ComboMF combo({grid,grid});
    std::array<std::array<double,3>,256> values;
    #pragma omp parallel for num_threads(4)
    for (int i=0;i<256;++i) values[i]=Query(combo,0.25L,0.5L,1.25L);
    const auto expected = Affine(0.25L,0.5L,1.25L);
    for (const auto& row : values)
        for (size_t a=0;a<3;++a) EXPECT_NEAR(row[a],2*expected[a],1e-13);
}

TEST_F(GridFields, DeflectorAndCloneMatchUniformFieldForElectronAndNucleus) {
    GridFile input(Grid(0,true));
    SmartPtr<MagneticField> grid = new SavedMF(input.path.string(),SavedMF::FILLING_NEAREST);
    SmartPtr<MagneticField> baseline = new ConstantField({0,1e-9,0},units.Mpc);
    SmartPtr<MagneticField> combo = new ComboMF({grid});
    Deflection3D deflector(combo,10), control(baseline,10);
    std::unique_ptr<DeflectionInteraction> clone(deflector.Clone());
    for (bool nucleus : {false,true}) {
        Particle p = nucleus ? Particle(56u,0.01) : Particle(Electron,0.01);
        p.Energy = 1e19*units.eV; p.Weight = 2.5;
        Particle expected=p, copied=p;
        Randomizer r(17), r2(17), r3(17);
        EXPECT_TRUE(control.Propagate(0.01L*units.Mpc,expected,r));
        EXPECT_TRUE(deflector.Propagate(0.01L*units.Mpc,p,r2));
        EXPECT_TRUE(clone->Propagate(0.01L*units.Mpc,copied,r3));
        for (size_t a=0;a<3;++a) {
            EXPECT_EQ(p.X[a],expected.X[a]); EXPECT_EQ(p.Pdir[a],expected.Pdir[a]);
            EXPECT_EQ(p.X[a],copied.X[a]); EXPECT_EQ(p.Pdir[a],copied.Pdir[a]);
        }
        EXPECT_DOUBLE_EQ(p.Energy,expected.Energy); EXPECT_DOUBLE_EQ(p.Weight,2.5);
        EXPECT_EQ(p.Time.t(),expected.Time.t());
        if (nucleus) { EXPECT_EQ(p.MassNumber(),56u); EXPECT_EQ(p.AtomicNumber(),26u); }
    }
}
} // namespace
