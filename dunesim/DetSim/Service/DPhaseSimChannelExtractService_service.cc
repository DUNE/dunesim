// DPhaseSimChannelExtractService.cxx

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "DPhaseSimChannelExtractService.h"
#include <string>
#include "Geometry/Geometry.h"
#include "Simulation/SimChannel.h"

//#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
//#include "artextensions/SeedService/SeedService.hh"
//#include "CLHEP/Random/RandGaussQ.h"

//#undef UseSeedService

using std::string;

//**********************************************************************

DPhaseSimChannelExtractService::
DPhaseSimChannelExtractService(const fhicl::ParameterSet& pset, art::ActivityRegistry&)
: m_ntick(m_pfft->FFTSize()) 
{
  
  if( !(m_ntick > 0 && !(m_ntick & (m_ntick-1))) )
    {
      throw cet::exception("DPhaseSimChannelExtractService")
	<< "FFT size is not a power of 2. ";
    }
      
  // get LEM gain per view
  fDPGainPerView  = pset.get<float> ("DPGainPerView");
  //fRedENC         = pset.get<float> ("RedENC");
  mf::LogInfo("DPhaseSimChannelExtractService")<<" Gain per view "<<fDPGainPerView;

  /*
  // for adding noise fluctuations
#ifdef UseSeedService
  art::ServiceHandle<artext::SeedService> seedSvc;
  int seed = seedSvc->getSeed("DPhaseSimChannelExtractService");
#else
  int seed = 26012016;
#endif
  art::EngineCreator ecr;
  m_pran = &ecr.createEngine(seed, "HepJamesRandom", "DPhaseSimChannelExtractService");
  */
}

//**********************************************************************

int DPhaseSimChannelExtractService::
extract(const sim::SimChannel* psc, AdcSignalVector& sigs) const {

  // clear and resize temporary ADC buffer
  sigs.clear();
  sigs.resize(m_ntick, 0.0);

  if ( psc == nullptr ) return 0;

  // get the channel number 
  unsigned int chan = psc->Channel();
  
  //CLHEP::RandGaussQ rGauss(*m_pran, 0.0, fRedENC);
  
  for ( size_t itck=0; itck<sigs.size(); ++itck ) 
    {
      sigs[itck] = fDPGainPerView * psc->Charge(itck);
      //if(fRedENC > 0)
      //sigs[itck] += rGauss.fire();
    }
  
  // perform convolution
  m_psss->Convolute(chan, sigs);
  
  
  return 0;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DPhaseSimChannelExtractService, SimChannelExtractService)

//**********************************************************************
