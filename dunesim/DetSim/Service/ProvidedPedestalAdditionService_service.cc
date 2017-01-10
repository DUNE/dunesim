// ProvidedPedestalAdditionService.cxx

#include "dune/DetSim/Service/ProvidedPedestalAdditionService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "art/Framework/Services/Optional/TFileService.h"
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

ProvidedPedestalAdditionService::
ProvidedPedestalAdditionService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_RandomSeed(0), m_LogLevel(1),
  m_PedNoiseHist(nullptr),
  m_pran(nullptr),
  m_PedestalProvider(art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider()) {
  const string myname = "ProvidedPedestalAdditionService::ctor: ";
  m_NoiseScale = pset.get<float>("NoiseScale");
  bool haveSeed = pset.get_if_present<int>("RandomSeed", m_RandomSeed);
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  if ( m_RandomSeed == 0 ) haveSeed = false;
  if ( haveSeed ) {
    if ( m_LogLevel > 0 ) cout << myname << "WARNING: Using hardwired seed." << endl;
    m_pran = new HepJamesRandom(m_RandomSeed);
  } else {
    string rname = "ProvidedPedestalAdditionService";
    if ( m_LogLevel > 0 ) cout << myname << "Using NuRandomService." << endl;
    art::ServiceHandle<NuRandomService> seedSvc;
    m_pran = new HepJamesRandom;
    if ( m_LogLevel > 0 ) cout << myname << "    Initial seed: " << m_pran->getSeed() << endl;
    seedSvc->registerEngine(NuRandomService::CLHEPengineSeeder(m_pran), rname);
  }
  if ( m_LogLevel > 0 ) cout << myname << "  Registered seed: " << m_pran->getSeed() << endl;
  art::ServiceHandle<art::TFileService> tfs;
  m_PedNoiseHist  = tfs->make<TH1F>("PedNoise", ";Pedestal noise  (ADC);", 1000,  -10., 10.);
}

//**********************************************************************

ProvidedPedestalAdditionService::~ProvidedPedestalAdditionService() {
  const string myname = "ProvidedPedestalAdditionService::dtor: ";
  if ( m_LogLevel > 0 ) {
    cout << myname << "Deleting random engine with seed " << m_pran->getSeed() << endl;
  }
  delete m_pran;
}

//**********************************************************************

int ProvidedPedestalAdditionService::
addPedestal(Channel chan, AdcSignalVector& sigs, float& ped, float& pedrms) const {
  // Fetch the pedestal for this channel.
  float ped_mean = m_PedestalProvider.PedMean(chan);
  float ped_rms =  m_PedestalProvider.PedRms(chan);
  for ( unsigned int itck=0; itck<sigs.size(); ++itck ) {
    sigs[itck] += ped_mean;
  }
  if ( ped_rms > 0 && m_NoiseScale > 0 ) {
    CLHEP::RandGaussQ rGauss_Ped(*m_pran, 0.0, m_NoiseScale*ped_rms);
    for ( unsigned int itck=0; itck<sigs.size(); ++itck ) {
      double ped_variation = rGauss_Ped.fire();
      m_PedNoiseHist->Fill(ped_variation);
      sigs[itck] += ped_variation;
    }
  }
  ped = ped_mean;
  pedrms = ped_rms;
  return 0;
}

//**********************************************************************

ostream& ProvidedPedestalAdditionService::print(ostream& out, string prefix) const {
  out << prefix << "ProvidedPedestalAdditionService:" << endl;
  out << prefix << "  Noise scale: " << m_NoiseScale;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(ProvidedPedestalAdditionService, PedestalAdditionService)

//**********************************************************************
