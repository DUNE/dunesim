// WhiteChannelNoiseService.cxx

#include "dune/DetSim/Service/WhiteChannelNoiseService.h"
#include "larcore/Geometry/Geometry.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "TH1F.h"

using std::ostream;
using std::cout;
using std::endl;
using std::string;
using rndm::NuRandomService;
using CLHEP::HepJamesRandom;

//**********************************************************************

WhiteChannelNoiseService::
WhiteChannelNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_RandomSeed(0), m_LogLevel(1),
  fNoiseHist(nullptr),
  m_pran(nullptr) {
  const string myname = "WhiteChannelNoiseService::ctor: ";
  // Fetch parameters.
  bool haveSeed = pset.get_if_present<int>("RandomSeed", m_RandomSeed);
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  if ( m_RandomSeed == 0 ) haveSeed = false;
  // Create random number engine.
  if ( haveSeed ) {
    if ( m_LogLevel > 0 ) cout << myname << "WARNING: Using hardwired seed." << endl;
    m_pran = new HepJamesRandom(m_RandomSeed);
  } else {
    string rname = "WhiteChannelNoiseService";
    if ( m_LogLevel > 0 ) cout << myname << "Using NuRandomService." << endl;
    art::ServiceHandle<NuRandomService> seedSvc;
    m_pran = new HepJamesRandom;
    if ( m_LogLevel > 0 ) cout << myname << "    Initial seed: " << m_pran->getSeed() << endl;
    seedSvc->registerEngine(NuRandomService::CLHEPengineSeeder(m_pran), rname);
  }
  if ( m_LogLevel > 0 ) cout << myname << "  Registered seed: " << m_pran->getSeed() << endl;
  art::ServiceHandle<art::TFileService> tfs;
  fNoiseHist     = tfs->make<TH1F>("Noise", ";Noise  (ADC);", 1000,   -10., 10.);
}

//**********************************************************************

WhiteChannelNoiseService::~WhiteChannelNoiseService() {
  const string myname = "WhiteChannelNoiseService::dtor: ";
  if ( m_LogLevel > 0 ) {
    cout << myname << "Deleting random engine with seed " << m_pran->getSeed() << endl;
  }
  delete m_pran;
}

//**********************************************************************

int WhiteChannelNoiseService::addNoise(Channel chan, AdcSignalVector& sigs) const {
  art::ServiceHandle<util::SignalShapingServiceDUNE> sss;
  float fASICGain      = sss->GetASICGain(chan);
  double fShapingTime   = sss->GetShapingTime(chan);
  std::map< double, int > fShapingTimeOrder;
  fShapingTimeOrder = { {0.5, 0}, {1.0, 1}, {2.0, 2}, {3.0, 3} };
  DoubleVec fNoiseFactVec;
  auto tempNoiseVec = sss->GetNoiseFactVec();
  if ( fShapingTimeOrder.find(fShapingTime) != fShapingTimeOrder.end() ) {
    size_t i = 0;
    fNoiseFactVec.resize(2);
    for ( auto& item : tempNoiseVec ) {
      fNoiseFactVec[i] = item.at(fShapingTimeOrder.find( fShapingTime )->second);
      fNoiseFactVec[i] *= fASICGain/4.7;
      ++i;
    }
  } else {
    throw cet::exception("WhiteChannelNoiseService")
      << "\033[93m"
      << "Shaping Time received from signalservices_dune.fcl is not one of allowed values"
      << std::endl
      << "Allowed values: 0.5, 1.0, 2.0, 3.0 usec"
      << "\033[00m"
      << std::endl;
  }
  art::ServiceHandle<geo::Geometry> geo;
  const geo::View_t view = geo->View(chan);
#ifdef UseSeedService
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine& engine = rng->getEngine("WhiteChannelNoiseService");
#else
  CLHEP::HepRandomEngine& engine = *m_pran;
#endif
  CLHEP::RandGaussQ rGauss_Ind(engine, 0.0, fNoiseFactVec[0]);
  CLHEP::RandGaussQ rGauss_Col(engine, 0.0, fNoiseFactVec[1]);
  for ( AdcSignal& sig : sigs ) {
    double tnoise = 0.0;
    if      ( view==geo::kU ) tnoise = rGauss_Ind.fire();
    else if ( view==geo::kV ) tnoise = rGauss_Ind.fire();
    else                      tnoise = rGauss_Col.fire();
    sig += tnoise;
  }
  return 0;
}

//**********************************************************************

ostream& WhiteChannelNoiseService::print(ostream& out, string prefix) const {
  string myprefix = prefix + "  ";
  out << myprefix << "      Random seed: " << m_pran->getSeed() << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(WhiteChannelNoiseService, ChannelNoiseService)

//**********************************************************************
