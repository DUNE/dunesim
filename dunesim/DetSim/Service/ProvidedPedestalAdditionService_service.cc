// ProvidedPedestalAdditionService.cxx

#include "dune/DetSim/Service/ProvidedPedestalAdditionService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "artextensions/SeedService/SeedService.hh"
#include "art/Framework/Core/EngineCreator.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "TH1F.h"

using std::ostream;
using std::endl;
using std::string;

#undef UseSeedService

//**********************************************************************

ProvidedPedestalAdditionService::
ProvidedPedestalAdditionService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_PedNoiseHist(nullptr),
  m_PedestalProvider(art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider()) {
  m_NoiseScale = pset.get<float>("NoiseScale");
#ifdef UseSeedService
  art::ServiceHandle<artext::SeedService> seedSvc;
  int seed = seedSvc->getSeed("ProvidedPedestalAdditionService");
#else
  int seed = 1007;
#endif
  art::ServiceHandle<art::TFileService> tfs;
  m_PedNoiseHist  = tfs->make<TH1F>("PedNoise", ";Pedestal noise  (ADC);", 1000,  -10., 10.);
  art::EngineCreator ecr;
  m_pran = &ecr.createEngine(seed, "HepJamesRandom", "ProvidedPedestalAdditionService");
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
