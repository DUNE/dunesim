// GenericSimChannelExtractService.cxx

#include "GenericSimChannelExtractService.h"
#include <string>
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimChannel.h"

using std::string;

//**********************************************************************

GenericSimChannelExtractService::
GenericSimChannelExtractService(const fhicl::ParameterSet&, art::ActivityRegistry&)
: m_ntick(m_pfft->FFTSize()) {
  if ( m_ntick%2 != 0 )
    throw cet::exception("GenericSimChannelExtractService")
          << "FFT size is not a power of 2.";
}

//**********************************************************************

int GenericSimChannelExtractService::
extract(const sim::SimChannel* psc, AdcSignalVector& sigs) const {
  sigs.clear();
  sigs.resize(m_ntick, 0.0);
  if ( psc == nullptr ) return 0;
  for ( size_t itck=0; itck<sigs.size(); ++itck ) {
    sigs[itck] = psc->Charge(itck);
  }
  unsigned int chan = psc->Channel();
  m_psss->Convolute(chan, sigs);
  return 0;
}

//**********************************************************************

std::ostream& GenericSimChannelExtractService::
print(std::ostream& out, std::string prefix) const {
  out << prefix << "GenericSimChannelExtractService";
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(GenericSimChannelExtractService, SimChannelExtractService)

//**********************************************************************
