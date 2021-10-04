// DPhaseSimChannelExtractService.cxx

#include <iostream>
#include "dune/ArtSupport/DuneToolManager.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "DPhaseSimChannelExtractService.h"

using std::cout;
using std::endl;
using std::string;

//**********************************************************************

DPhaseSimChannelExtractService::
DPhaseSimChannelExtractService(const fhicl::ParameterSet& pset, art::ActivityRegistry&)
: m_ntick(m_pfft->FFTSize())
{
  const string myname = "DPhaseSimChannelExtractService::ctor: ";
  if( !(m_ntick > 0 && !(m_ntick & (m_ntick-1))) )
    {
      throw cet::exception("DPhaseSimChannelExtractService")
	<< "FFT size is not a power of 2. ";
    }
  
  /*
  m_CrpGainToolName = pset.get<string>("CrpGainToolName");  
  DuneToolManager* ptm = DuneToolManager::instance("");
  if ( ptm == nullptr ) 
    {
      throw cet::exception("DPhaseSimChannelExtractService")
	<< "Unable to retrieve tool manaager.";
    } 
  
  m_CrpGainTool = ptm->getPrivate<CrpGainSimTool>(m_CrpGainToolName);
  
  if( !m_CrpGainTool )
    {
      throw cet::exception("DPhaseSimChannelExtractService")
	<< "Unable to retrieve "<<m_CrpGainToolName<<".";
    }
  */
  //fDPGainPerView  = pset.get<float> ("DPGainPerView");
  //mf::LogInfo("DPhaseSimChannelExtractService")<<" Gain per view "<<fDPGainPerView;
}

//**********************************************************************

int DPhaseSimChannelExtractService::
extract(detinfo::DetectorClocksData const& clockData,
        const sim::SimChannel* psc, AdcSignalVector& sigs) const {

  if ( psc == nullptr ) return 0;

  // get the channel number
  unsigned int chan = psc->Channel();
  
  // clear and resize temporary ADC buffer
  sigs.clear();
  sigs.resize(m_ntick, 0.0);

  for ( size_t itck=0; itck<sigs.size(); ++itck ){
    sigs[itck] = m_crpgain->viewCharge( psc, itck );
  }
  
  // perform convolution
  m_psss->Convolute(clockData, chan, sigs);
  
  return 0;
}

//**********************************************************************

std::ostream& DPhaseSimChannelExtractService::print(std::ostream& out, std::string prefix) const {
  out << prefix << "DPhaseSimChannelExtractService";
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DPhaseSimChannelExtractService, SimChannelExtractService)

//**********************************************************************
