#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#include "ParticleStack.h"
#include "Utils.h"
#include "Cosmology.h"
#include "PropagationEngine.h"
#include <iostream>
#include <GZK.h>
#include <Inoue12IROSpectrum.h>
#include "PPP.h"
#include "ProtonPP.h"
#include "TableBackgrounds.h"
#include "ICS.h"
#include "Test.h"
#include "NeutronDecay.h"
#include "Deflection3D.h"
#include "SteckerEBL.h"
#include "Stecker16Background.h"
#include <omp.h>
#include <stdlib.h>
#include "MathUtils.h"
#include "Randomizer.h"
//#include "Jet.h"
#include <math.h>

using namespace std;
using namespace mcray;
using namespace Utils;
using namespace Backgrounds;
using namespace Interactions;

void testCPP(){
	std::cerr<<"Test from CPP success!"<<std::endl;
}

void testMCpropag(){
	double aOmegaVac=0.7;
	double aHubbleKm_sec_Mpc=70;
	if(!cosmology.IsInitialized())
		cosmology.Init(0.02 + 10, aOmegaVac, aHubbleKm_sec_Mpc);
	cerr<<"Cosmology init OK!"<<endl;
	cerr<<"units.Eunit = "<<units.Eunit<<endl;//MeV
	cerr<<"TEST OK"<<endl;
}

int add(int i, int j) {
	return i + j;
}

void cosmologyInit(double fZmax, double aOmegaVac, double aHubbleKm_sec_Mpc) {
	if(!cosmology.IsInitialized())
		cosmology.Init(fZmax, aOmegaVac, aHubbleKm_sec_Mpc);
	else
		cerr<<"Cosmology already initialized. Nothing do."<<endl;
}

class InteractionAdapter{ // Stub interface for addInteraction of PropagationEngine
	public:
		~InteractionAdapter(){
			if(destructorFlag)return;
			cerr<<"We in destructor!"<<endl;
			cerr<<"Please, do not create objects without purpose, it can raise memory leak!"<<endl;
			try{
				if(rndInt){delete rndInt;rndInt=nullptr;}
				if(celInt){delete celInt;celInt=nullptr;}
				if(dflInt){delete dflInt;dflInt=nullptr;}
			}catch(...){
				cerr<<"What The Failure error in InteractionAdapter?!"<<endl;
			}
		}
		RandomInteraction* rndInt=nullptr;
		CELInteraction* celInt=nullptr;
		DeflectionInteraction* dflInt=nullptr;
		int fInteractionType;
		std::string name="Unknown";
		bool destructorFlag=false;
};


/*class PyMagneticField{ // Magnetic Field adapter for avoid memory leak
	public:
		~PyMagneticField(){
			if(destructorFlag)return;
			try{
				if(pB){delete pB;pB=nullptr;}
			}catch(...){
				cerr<<"What The Failure error in MagneticField adapter?!"<<endl;
			}
		}
		MagneticField* pB=nullptr;
		std::string name="Unknown";
		bool destructorFlag=false;
} //*/

//  PYBIND

namespace py = pybind11;

class CustomMF : public MagneticField{
	public:
		CustomMF(py::function &aGetB,py::function &aGetScale){
			fGetB=aGetB;
			fGetScale=aGetScale;
		};
		virtual void GetValueGauss(const coord_type *x, const CosmoTime &aTime,
								std::vector<double> &outValue) const {
			py::list xv;
			std::vector<double> result_py = fGetB(x[0],x[1],x[2],aTime).cast<std::vector<double>>();
			outValue[0] = result_py[0];
			outValue[1] = result_py[1];
			outValue[2] = result_py[2];
		};
		virtual double MinVariabilityScale(const CosmoTime &aTime) const {
			return fGetScale(aTime).cast<double>();
		};
	private:
		py::function fGetB;
		py::function fGetScale;
};

PYBIND11_MODULE(pybeam, m) {
	m.doc() = R"pbdoc(
		Pybind11 plugin for CRbeam
		-----------------------

		Tests:
		   add
		   testCPP
		   testMCpropag
		   subtract
	)pbdoc";

//  MCPROPAG UNITS

	py::module pyunits = m.def_submodule("units");
	pyunits.attr("eV")=units.eV;
	pyunits.attr("phTemperature_mult")=Units::phTemperature_mult;
	pyunits.attr("Eunit")=units.Eunit;
	pyunits.attr("cm3")=units.cm3;
	pyunits.attr("Vunit")=units.Vunit;
	pyunits.attr("Mpc")=units.Mpc;
	pyunits.attr("kpc")=units.kpc;

	py::module pyParticleTypes = m.def_submodule("particle_types");
	pyParticleTypes.attr("Electron")=(int)Electron;
	pyParticleTypes.attr("Positron")=(int)Positron;
	pyParticleTypes.attr("Photon")=(int)Photon;
	pyParticleTypes.attr("NeutrinoE")=(int)NeutrinoE;
	pyParticleTypes.attr("NeutrinoM")=(int)NeutrinoM;
	pyParticleTypes.attr("NeutrinoT")=(int)NeutrinoT;
	pyParticleTypes.attr("NeutrinoAE")=(int)NeutrinoAE;
	pyParticleTypes.attr("NeutrinoAM")=(int)NeutrinoAM;
	pyParticleTypes.attr("NeutrinoAT")=(int)NeutrinoAT;
	pyParticleTypes.attr("Neutron")=(int)Neutron;
	pyParticleTypes.attr("Proton")=(int)Proton;

//  MCPROPAG BACKGROUND UTILS

	py::module PyBackgroundUtils = m.def_submodule("BackgroundUtils");
	PyBackgroundUtils.def("CalcIntegralDensity",&BackgroundUtils::CalcIntegralDensity, R"pbdoc(
		Example: pybeam.BackgroundUtils.CalcIntegralDensity(backgr)*pybeam.units.cm3
	)pbdoc",py::arg("Background"), py::arg("aRelError") = 1e-3);
	PyBackgroundUtils.def("Print",[](IBackground * b, double z, std::string s, int n) {
		std::ofstream backgrOut;
		backgrOut.open(s.c_str(),std::ios::out);
		BackgroundUtils::Print(b, z, backgrOut, n);
	}, R"pbdoc(
		Example: pybeam.BackgroundUtils.Print(backgr,0,'output_filename',100)
	)pbdoc",py::arg("Background"), py::arg("z"), py::arg("filename"), py::arg("nIntervals")=100);

//  MCPROPAG COSMOLOGY INIT

	m.def("cosmologyInit", &cosmologyInit, R"pbdoc(
		Init cosmology, call it only once after import
	)pbdoc",py::arg("z") = 0.02+10, py::arg("aOmegaVac") = 0.7, py::arg("aHubbleKm_sec_Mpc") = 70);

//  MCPROPAG OUTPUT 3D

	py::class_<IOutput>(m, "IOutput"); // Do not call it
	py::class_<RawOutput3D,IOutput>(m, "RawOutput3D")
		.def(py::init<std::string,bool,bool,bool,bool,bool,bool>())
		.def("SetOutputDir", &RawOutput3D::SetOutputDir,py::arg("dirname"), py::arg("OverwriteExisting")=false);

//  MCPROPAG BACKGROUND CLASSES

	py::class_<IBackground,std::unique_ptr<IBackground, py::nodelete>>(m, "IBackground")
		.def("n", &IBackground::n,R"pbdoc(
				n(double aE, double aZ)
			)pbdoc",py::arg("aE"), py::arg("aZ"))
		.def("Name", &IBackground::Name)
		.def("n", &IBackground::n,R"pbdoc(
				n(double aE, double aZ)
			)pbdoc",py::arg("aE"), py::arg("aZ"))
		.def("MinE", &IBackground::MinE)
		.def("MaxE", &IBackground::MaxE)
		.def("MinZ", &IBackground::MinZ)
		.def("MaxZ", &IBackground::MaxZ)
		.def("__repr__",
		[](const IBackground &a) {
			return "<pybeam.IBackground '" + a.Name() + "'>";
		}
		);

	py::class_<CompoundBackground,IBackground,std::unique_ptr<CompoundBackground, py::nodelete>>(m, "CompoundBackground")
		.def(py::init<>())
		.def("AddComponent", &CompoundBackground::AddComponent, R"pbdoc(
				CompoundBackground.AddComponent(component, weight = 1)
				Add background component
			)pbdoc",py::arg("Background"), py::arg("weight") = 1)
		.def("__repr__",
		[](const CompoundBackground &a) {
			return "<pybeam.CompoundBackground '" + a.Name() + "'>";
		}
		);

	py::class_<Inoue12BaselineIROSpectrum,IBackground>(m, "Inoue12BaselineIROSpectrum")
		.def(py::init<>());
	py::class_<MatrixBackground,IBackground,std::unique_ptr<MatrixBackground, py::nodelete>>(m, "MatrixBackground")
		.def(py::init<std::string, std::string, bool, bool>(),R"pbdoc(
				MatrixBackground from file in \"tables\" directory
				Example:
				b = pybeam.MatrixBackground("Inoue12Baseline", "EBL_Inoue/EBL_proper_baseline.dat", false, true)
			)pbdoc",py::arg("Name"), py::arg("TableFile"), py::arg("IsComoving")=false, py::arg("ExtendToZero")=true);

	py::class_<PlankBackground,IBackground,std::unique_ptr<PlankBackground, py::nodelete>>(m, "PlankBackground")
		.def(py::init<double,double,double,double,double,bool>(),R"pbdoc(
				PlankBackground is CMB background
				Example:

				cmbTemp = 2.73/pybeam.units.phTemperature_mult/pybeam.units.Eunit
				b = pybeam.PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp, 0, fZmax + 1.)
			)pbdoc",py::arg("aT0"), py::arg("aEmin"), py::arg("aEmax"), py::arg("aMinZ")=0, py::arg("aMaxZ")=10000, py::arg("aFixedTemperature")=false);

	py::class_<CuttedBackground,IBackground,std::unique_ptr<CuttedBackground, py::nodelete>>(m, "CuttedBackground")
		.def(py::init<IBackground*, double, double, double, double>(),R"pbdoc(
				Cut background
				Example:
				b = pybeam.CuttedBackground(backgr, 0, 100)
			)pbdoc",py::arg("Background"), py::arg("MinE")=0, py::arg("MaxE")=1e300, py::arg("MinZ")=0, py::arg("MaxZ")=1e300);

	py::class_<HighRedshiftBackgrExtension,IBackground,std::unique_ptr<HighRedshiftBackgrExtension, py::nodelete>>(m, "HighRedshiftBackgrExtension")
		.def(py::init<IBackground*, double, double, double, double>(),R"pbdoc(
				HighRedshiftBackgrExtension
				Example:
				b = pybeam.HighRedshiftBackgrExtension(backgr, fExtDeltaZconst, fExtPowerLow, fExtDeltaZexp)
			)pbdoc",py::arg("Background"), py::arg("DeltaZconst"), py::arg("PowerLow"), py::arg("DeltaZexp"), py::arg("Zmax")=1000);

	py::class_<Function>(m, "Function"); // Do not call it
	py::class_<ConstFunction,Function,std::unique_ptr<ConstFunction, py::nodelete>>(m, "ConstFunction")
		.def(py::init<double,double,double>(),R"pbdoc(
				ConstFunction stub for interface
				Example:
				k1 = pybeam.ConstFunction(centralE1)
			)pbdoc",py::arg("aValue")=0, py::arg("Xmin")=-std::numeric_limits<double>::max(), py::arg("Xmax")=std::numeric_limits<double>::max() )
		.def("__call__",&ConstFunction::f,R"pbdoc(
				Call const function is a stub and always return init value
			)pbdoc",py::arg("x")=0) // Const function is a stub
		.def("__float__",[](const ConstFunction &a){return a.f(1);}); // Const function is a stub

	py::class_<BackgroundIntegral>(m, "BackgroundIntegral"); // Do not call it
	py::class_<MonochromaticBackgroundIntegral,BackgroundIntegral,std::unique_ptr<MonochromaticBackgroundIntegral, py::nodelete>>(m, "MonochromaticBackgroundIntegral")
		.def(py::init<Function&,Function&,double,double>(),R"pbdoc(
				MonochromaticBackgroundIntegral with Functional Concentration and Energy
				Example:
				backgrI = pybeam.MonochromaticBackgroundIntegral(c1, k1, logStepK, epsRel)
			)pbdoc",py::arg("Concentration"), py::arg("Energy"), py::arg("LogStepE"), py::arg("EpsRel")) //*/
		.def(py::init([](double centralE1,double n1,double logStepK,double epsRel){
			ConstFunction k1(centralE1);
			ConstFunction c1(n1);
			return new MonochromaticBackgroundIntegral(c1, k1, logStepK, epsRel);
		}),R"pbdoc(
				MonochromaticBackgroundIntegral with constant float Concentration and Energy
				Same result but ConstFunction is baked in.
				Example:
				backgrI = pybeam.MonochromaticBackgroundIntegral(c1, k1, logStepK, epsRel)
			)pbdoc",py::arg("Concentration"), py::arg("Energy"), py::arg("LogStepE"), py::arg("EpsRel"));

	py::class_<ContinuousBackgroundIntegral,BackgroundIntegral,std::unique_ptr<ContinuousBackgroundIntegral, py::nodelete>>(m, "ContinuousBackgroundIntegral")
		.def(py::init<IBackground&,double,double,double,double>(),R"pbdoc(
				ContinuousBackgroundIntegral
				Example:
				backgrI = pybeam.ContinuousBackgroundIntegral(backgr, stepZ, logStepK, fZmax, epsRel)
			)pbdoc",py::arg("Background"), py::arg("StepZ"), py::arg("LogStepE"), py::arg("ZmaxLimit"), py::arg("EpsRel") = 0);

//  MCPROPAG RANDOMIZER
	py::class_<Randomizer,std::unique_ptr<Randomizer, py::nodelete>>(m, "Randomizer")
		.def(py::init<unsigned long int>(),R"pbdoc(
				Randomizer
				Example:
				rand = pybeam.Randomizer(seed)
			)pbdoc",py::arg("seed"))
		.def(py::init<>([](){
			struct timeval tv;
			gettimeofday(&tv, NULL);
			unsigned long int seed = tv.tv_sec * 1000 + tv.tv_usec / 1000;
			return new Randomizer(seed);
		}),R"pbdoc(
				Randomizer with automatic seed from current time
				Example:
				rand = pybeam.Randomizer()
			)pbdoc")
		.def("Rand", &Randomizer::Rand)
		.def("RandZero", &Randomizer::RandZero)
		.def("RandGauss", &Randomizer::RandGauss)
		.def("CreateIndependent", &Randomizer::CreateIndependent);

//  MCPROPAG ENGINE HELPERS CLASSES: FILTERS AND ETC
	py::class_<IFilter>(m, "IFilter"); // Do not call it
	py::class_<EnergyBasedFilter,IFilter,std::unique_ptr<EnergyBasedFilter, py::nodelete>>(m, "EnergyBasedFilter")
		.def(py::init<double,double>(),R"pbdoc(
				EnergyBasedFilter
				Example:
				resultFilter = pybeam.EnergyBasedFilter(Emin,np.pi)
			)pbdoc",py::arg("Emin")=0, py::arg("Xmin")=M_PI);
	py::class_<ParticleTypeFilter,IFilter,std::unique_ptr<ParticleTypeFilter, py::nodelete>>(m, "ParticleTypeFilter")
		.def(py::init([](const std::vector<int> &v){
			ParticleType* arr = new ParticleType[v.size()+1];
			for(int i=0;i<v.size();i++){
				if(v[i]<Electron || v[i]>Proton){
					std::cerr<<v[i]<<" -- bad particle type"<<std::endl;
					throw pybind11::value_error(std::to_string(v[i])+" is bad particle type");
				}else arr[i] =(ParticleType)v[i];
			}
			arr[v.size()] = ParticleTypeEOF;
			return new ParticleTypeFilter(arr);
		}),R"pbdoc(
				ParticleTypeFilter
				Example:
				discardedParticles = pybeam.ParticleTypeFilter([0,1,2]) # Electron, Positron, Photon
			)pbdoc")
		.def("__repr__",
		[](const ParticleTypeFilter &a) {
			std::string repr="<pybeam.ParticleTypeFilter([";
			bool flag = false;
			for(int i=0;i<ParticleTypeEOF;i++){
				Particle particle( (ParticleType)i, 0.1);
				if( !a.Pass(particle) ){
					repr+=std::to_string(i)+",";
					flag = true;
				}
			}
			if(flag)return repr.substr(0,repr.length()-1) + "])";
			return repr + "])";
		}
		);

	py::class_<IResult>(m, "IResult")
		.def("ContinuePropagation", &IResult::ContinuePropagation)
		.def("GetMinTravelTimeLeft", &IResult::GetMinTravelTimeLeft)
		.def("AbortRequested", &IResult::AbortRequested); // Do not call it
	py::class_<Result,IResult,std::unique_ptr<Result, py::nodelete>>(m, "Result")
		.def(py::init<IFilter*,bool>(),R"pbdoc(
				Result
				Example:
				result = pybeam.Result(resultFilter,True)
			)pbdoc")
		.def("SetAutoFlushInterval", &Result::SetAutoFlushInterval)
		.def("EnableFlushOnCtrlC", &Result::EnableFlushOnCtrlC)
		.def("AddOutput", &Result::AddOutput)
		.def("Flush", &Result::Flush,py::arg("final")=true)
		.def("SetEndTime", &Result::SetEndTime);

	py::class_<ParticleStack>(m, "ParticleStack")
		.def(py::init<int>(),R"pbdoc(
				ParticleStack
				Example:
				particles = pybeam.ParticleStack()
			)pbdoc",py::arg("aDebugOutputPeriod")=-1)
		.def("AddPrimary", &ParticleStack::AddPrimary);

	py::class_<PropagationEngine>(m, "PropagationEngine") // PROPAGATION ENGINE
		.def(py::init<ParticleStack&,IResult&,unsigned long int>(),R"pbdoc(
				PropagationEngine
				Example:
				pe = pybeam.PropagationEngine(particles, result, rand.CreateIndependent())
			)pbdoc")
		.def("SetThinning", &PropagationEngine::SetThinning)
		.def("SetPropagationFilter", &PropagationEngine::SetPropagationFilter)
		.def("Run", &PropagationEngine::Run)
		.def("RunMultithread", &PropagationEngine::RunMultithread)
		.def("RunMultithreadReleaseOnly", &PropagationEngine::RunMultithreadReleaseOnly)
		.def("AddInteraction", [](PropagationEngine &a,InteractionAdapter* b){
			if(b->fInteractionType==0)a.AddInteraction(b->rndInt);
			else if(b->fInteractionType==1)a.AddInteraction(b->celInt);
			else if(b->fInteractionType==2)a.AddInteraction(b->dflInt);
			else throw pybind11::value_error("Not implemented interaction type: fInteractionType="+std::to_string(b->fInteractionType));
			b->destructorFlag=true; // Do not call destructor from Python
		});

	py::class_<IThinning>(m, "IThinning"); // Do not call it
//	py::class_<PhotonThinning,IThinning,std::unique_ptr<PhotonThinning, py::nodelete>>(m, "PhotonThinning") // TODO fix it
//		.def(py::init<double,double>(),R"pbdoc(
//				PropagationEngine
//				Example:
//				thinning = pybeam.PhotonThinning(alphaThinning, EThinning)
//			)pbdoc",py::arg("Alpha"),py::arg("Threshold"));
//
	// Для всяких взаимодействий строю костыльный огород
	py::class_<InteractionAdapter>(m, "InteractionAdapter") // Do not call it
		.def("__repr__", [](InteractionAdapter &a){
			if(a.fInteractionType==0)return "<pybeam RandomInteraction -> "+a.name+" >";
			else if(a.fInteractionType==1)return "<pybeam CELInteraction -> "+a.name+" >";
			else if(a.fInteractionType==2)return "<pybeam DeflectionInteraction -> "+a.name+" >";
			else throw pybind11::value_error("Not implemented interaction type");
		})
		.def("Propagate", [](InteractionAdapter &a,cosmo_time deltaT, Particle &aParticle, Randomizer &aRandomizer){
			if(a.fInteractionType==2)return a.dflInt->Propagate(deltaT, aParticle, aRandomizer);
			else throw pybind11::value_error("Methot Propagate implemented only for DeflectionInteraction");
		});

	m.def("ICSInteraction", [](BackgroundIntegral* backgrI,double Emin){
			InteractionAdapter* tmp=new InteractionAdapter();
			tmp->fInteractionType = 0;
			tmp->rndInt = new Interactions::ICS(backgrI, Emin);
			tmp->name = "Interactions::ICS with Emin="+std::to_string(Emin);
			return tmp;
		}, R"pbdoc(
				ICS interaction, returns InteractionAdapter instance
				Example:
				ics = pybeam.ICSInteraction(backgrI, EThinning)
			)pbdoc");
	m.def("GZK", [](BackgroundIntegral* backgrI,double seed){
			InteractionAdapter* tmp=new InteractionAdapter();
			tmp->fInteractionType = 0;
			tmp->rndInt = new GZK(backgrI, seed);
			tmp->name = "GZK with seed="+std::to_string(seed);
			return tmp;
		},R"pbdoc(
				GZK interaction, returns InteractionAdapter instance
				Example:
				gzk = pybeam.GZK(backgrI)
		)pbdoc",py::arg("Background"),py::arg("RandomSeed")=0);
	m.def("NeutronDecay", [](){
			InteractionAdapter* tmp=new InteractionAdapter();
			tmp->fInteractionType = 0;
			tmp->rndInt = new NeutronDecay();
			tmp->name = "Neutron Decay";
			return tmp;
		},R"pbdoc(
				NeutronDecay interaction, returns InteractionAdapter instance
				Example:
				ndecay = pybeam.NeutronDecay()
		)pbdoc");
	m.def("ProtonPPcel", [](BackgroundIntegral* backgrI){
			InteractionAdapter* tmp=new InteractionAdapter();
			tmp->fInteractionType = 1;
			tmp->celInt = new ProtonPPcel(backgrI);
			tmp->name = "ProtonPPcel";
			return tmp;
		},R"pbdoc(
				ProtonPPcel interaction, returns InteractionAdapter instance
				Example:
				ppcel = pybeam.ProtonPPcel(backgrI)
		)pbdoc",py::arg("Background"));
	m.def("Deflection3D", [](MagneticField * aMagnFieldGauss, uint aAccuracy){
			InteractionAdapter* tmp=new InteractionAdapter();
			tmp->fInteractionType = 2;
			tmp->dflInt = new Deflection3D(aMagnFieldGauss,aAccuracy);
			tmp->name = "Deflection3D with accuracy="+std::to_string(aAccuracy);
			return tmp;
		},R"pbdoc(
				Deflection3D interaction, returns InteractionAdapter instance
				Example:
				defl = pybeam.Deflection3D(MagnField,10)
		)pbdoc",py::arg("MagnField"),py::arg("Accuracy")=10);

	py::class_<MagneticField>(m, "MagneticField") // Do not call it
		.def("GetValueGauss", [](MagneticField &a,const std::vector<coord_type> &x,cosmo_time t){
			if(x.size()!=3)throw pybind11::value_error("The first argument must be with length 3 (x,y,z)");
			std::vector<double> B(3);
			a.GetValueGauss(&x[0],t,B); // See https://stackoverflow.com/questions/2923272/how-to-convert-vector-to-array
			return B;
		})
		.def("savetxt", [](MagneticField &a,std::string &fname,double maxMpc,double stepMpc){
			CosmoTime t(0);
			if(stepMpc<1e-15)stepMpc=a.MinVariabilityScale(t);
			std::ofstream mfout(fname.c_str());
			a.Print(maxMpc/3., mfout, stepMpc); // maxMpc TODO: Not implemented
			mfout.close();
		},R"pbdoc(
				Dump magnetic to file
				Example:
				mf.savetxt("/tmp/magneticField.out",10)
				if stepMpc not passed, it will be get via mf.MinVariabilityScale
		)pbdoc",py::arg("fname"),py::arg("maxMpc"),py::arg("stepMpc")=-1)
		.def("MinVariabilityScale", &MagneticField::MinVariabilityScale);
//	py::class_<TurbulentMF,MagneticField,std::unique_ptr<TurbulentMF, py::nodelete>>(m, "TurbulentMF")
//		.def(py::init<Randomizer&,double,double,double, double, double>(),R"pbdoc(
//				TurbulentMF with Kolmogorov spectrum
//				Example:
//				mf = pybeam.TurbulentMF(rand, 1, 1e-9, 1.01, 0.1)
//				Randomizer -- pybeam.Randomizer instance
//				Lcor_Mpc -- Extragalactic magnetic field correlation length in Mpc at z=0
//				Variance_Gauss -- Amplitude of magnetic field in Gauss
//				MultK -- Ratio of two neighbourhood modes in EGMF spectrum
//				MinVarScale -- Minimal variability scale of the MF as fraction of Lc
//				PMFCut_p=-1 -- reset magnetic field to 0 after *this param* corr lengths, debug feature, not need for physics
//			)pbdoc",py::arg("Randomizer"), py::arg("Lcor_Mpc"), py::arg("Variance_Gauss"), py::arg("MultK")=1.05, py::arg("MinVarScale")=0.1,py::arg("PMFCut_p")=-1); // TODO: Fix constructor args

	py::class_<SavedMF,MagneticField,std::unique_ptr<SavedMF, py::nodelete>>(m, "SavedMF")
		.def(py::init<const char *>(),R"pbdoc(
				SavedMF
				Example:
				mf = pybeam.SavedMF("/tmp/magneticField.out")
				File format:
				#x,Mpc	y,Mpc	z,Mpc	Bx	By	Bz
			)pbdoc",py::arg("fname"))
		.def("GetValueGauss_NEAREST", [](SavedMF &a,const std::vector<coord_type> &x,cosmo_time t){
			if(x.size()!=3)throw pybind11::value_error("The first argument must be with length 3 (x,y,z)");
			std::vector<double> B(3);
			a.GetValueGauss_NEAREST(&x[0],t,B); // See https://stackoverflow.com/questions/2923272/how-to-convert-vector-to-array
			return B;
		})
		.def("GetValueGauss_LINEAR", [](SavedMF &a,const std::vector<coord_type> &x,cosmo_time t){
			if(x.size()!=3)throw pybind11::value_error("The first argument must be with length 3 (x,y,z)");
			std::vector<double> B(3);
			a.GetValueGauss_LINEAR(&x[0],t,B); // See https://stackoverflow.com/questions/2923272/how-to-convert-vector-to-array
			return B;
		});

	py::class_<CustomMF,MagneticField,std::unique_ptr<CustomMF, py::nodelete>>(m, "CustomMF")
		.def(py::init<py::function&,py::function&>(),R"pbdoc(
				Custom MF implemented on Python side
				Example:
				def getValueGauss(x,y,z,t):
					return [0,0,0]
				def getScale(t):
					return 0.1*pybeam.units.Mpc
				mf = pybeam.CustomMF(getValueGauss,getScale)
				getValueGauss must return Python list with length 3
			)pbdoc",py::arg("getValueGauss"),py::arg("getScale"))
		.def("__repr__", [](CustomMF &a){
			return "TODO";
		});

	py::class_<CosmoTime>(m, "CosmoTime")
		.def(py::init<cosmo_time>(),R"pbdoc(
				CosmoTime class connects redshift and time
				Example:
				t = pybeam.CosmoTime(0.03) # By Redshift
				t0 = pybeam.CosmoTime() # z=0

				Moreover you can use this class for converting redshift to distance:
				Example:
				pybeam.cosmologyInit() # Make sure before work with CosmoTime
				D = -pybeam.CosmoTime(0.0024).t()/pybeam.units.kpc # D = 10260.737417323564
				or just
				D = -pybeam.CosmoTime(0.0024)/pybeam.units.kpc # D = 10260.737417323564
				CosmoTime->float casting returns CosmoTime.t()
				some arithmetic operations are casts CosmoTime->float automatically
			)pbdoc",py::arg("Redshift")=0)
		.def(py::init<CosmoTime>(),R"pbdoc(
				CosmoTime copy constructor
				Example:
				t = pybeam.CosmoTime(t) # By another CosmoTime
			)pbdoc")
		.def("__float__", [](const CosmoTime &a){return a.t();})
		.def("__neg__", [](const CosmoTime &a){return -a.t();})
		.def("__lt__", [](const CosmoTime &a,const CosmoTime &b){return a<b;})
		.def("__gt__", [](const CosmoTime &a,const CosmoTime &b){return !(a<b) && a.t()!=b.t();})
		.def("__eq__", [](const CosmoTime &a,const CosmoTime &b){return a.t()==b.t();})
		.def("__ne__", [](const CosmoTime &a,const CosmoTime &b){return a.t()!=b.t();})
		.def("__sub__", [](const CosmoTime &a,const CosmoTime &b){return a.t()-b.t();})
		.def("__iadd__", [](CosmoTime &a,cosmo_time b){return a+=b;})
		.def("__truediv__", [](const CosmoTime &a,cosmo_time b){return a.t()/b;})
		.def("t", &CosmoTime::t)
		.def("z", &CosmoTime::z)
		.def("setZ", &CosmoTime::setZ)
		.def("setT", &CosmoTime::setT)
		.def("__repr__", [](const CosmoTime &a){
			return "<pybeam.CosmoTime object: "+a.ToString()+">";
		});

	py::class_<Particle,std::unique_ptr<Particle, py::nodelete>>(m, "Particle", py::dynamic_attr())
		.def(py::init([](int a,cosmo_time b){
			return new Particle((ParticleType)a,b);
		}),R"pbdoc(
				Particle class
				Example:
				particle = pybeam.Particle(pybeam.particle_types.Proton, 0.0024)
			)pbdoc",py::arg("type"),py::arg("Zsource"))
		.def("copy",[](const Particle &a){
			Particle* p=new Particle(a.Type,a.fProductionTime.z());
			(*p)=a;
//			memcpy(p, &a, sizeof(Particle));
			return *p;
		})
		.def("Mass",[](const Particle &a){return a.Mass();})
		.def("ElectricCharge",[](const Particle &a){return a.ElectricCharge();})
		.def("beta", &Particle::beta)
		.def("SetParent", &Particle::SetParent)
		.def("PropagateFreely", &Particle::PropagateFreely)
		.def("JetOpenningAngle", &Particle::JetOpenningAngle)
		.def("LastDeflectionTime", &Particle::LastDeflectionTime)
		.def("ObservationAngle", &Particle::ObservationAngle)
		.def("DeflectionAngle", &Particle::DeflectionAngle)
		.def_readwrite("ParticleType",&Particle::Type)
		.def_readwrite("Weight",&Particle::Weight)
		.def_readwrite("Time",&Particle::Time)
		.def_readwrite("Energy",&Particle::Energy)
		.def_readwrite("Ninteractions",&Particle::Ninteractions)
		.def("getX",[](const Particle &a){
			std::vector<double> X(3);
			X[0]=a.X[0];
			X[1]=a.X[1];
			X[2]=a.X[2];
			return X;
		})
		.def("setX",[](Particle &a,std::vector<double> X){
			if(X.size()!=3)throw pybind11::value_error("The argument must be with length 3 (x,y,z)");
			a.X[0]=X[0];
			a.X[1]=X[1];
			a.X[2]=X[2];
		})
		.def("getP",[](const Particle &a){
			std::vector<double> X(3);
			X[0]=a.Pdir[0];
			X[1]=a.Pdir[1];
			X[2]=a.Pdir[2];
			return X;
		})
		.def("getStartP",[](const Particle &a){
			std::vector<double> X(3);
			X[0]=a.PdirStart[0];
			X[1]=a.PdirStart[1];
			X[2]=a.PdirStart[2];
			return X;
		})
		.def("setP",[](Particle &a,std::vector<double> X,bool asStart){
			if(X.size()!=3)throw pybind11::value_error("The argument must be with length 3 (Px,Py,Pz)");
			double norm = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
			a.Pdir[0]=X[0]/norm;
			a.Pdir[1]=X[1]/norm;
			a.Pdir[2]=X[2]/norm;
			if(asStart){
				a.PdirStart[0]=X[0]/norm;
				a.PdirStart[1]=X[1]/norm;
				a.PdirStart[2]=X[2]/norm;
			}
		},R"pbdoc(
				Set P direction. If |P|!=1 it will be automatically normalized
				asStart -- set this P as initial PdirStart, default is false
				Example:
				particle.setP([0.238445,-0.455763,-0.857569],True)
			)pbdoc",py::arg("X"),py::arg("asStart")=false)
//		.def("setTracer",[](Particle &a,const char * fname){ // TODO: implement it
//			if(a.Tracer){
//				cerr<<"Warning: previous tracer will be closed"<<endl;
//				if(a.Tracer->is_open())a.Tracer->close();
//				delete a.Tracer;
//			}
//			a.Tracer=new std::ofstream(fname);
//		})
		.def_readwrite("id",&Particle::id)
		.def_readwrite("fPrevId",&Particle::fPrevId)
//		.def_readwrite("isInTarget",&Particle::isInTarget) // TODO: not implemented
		.def("__repr__", [](Particle &a){
			return "<pybeam.Particle object: "+a.ToString()+">";
		});


	m.def("omp_get_max_threads", &omp_get_max_threads);

	m.def("testCPP", &testCPP, R"pbdoc(
		Test C++ call.
	)pbdoc");

	m.def("testMCpropag", &testMCpropag, R"pbdoc(
		Test C++ call with MCpropag features.
	)pbdoc");

#ifdef VERSION_INFO
	m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
	m.attr("__version__") = "dev";
#endif
}
