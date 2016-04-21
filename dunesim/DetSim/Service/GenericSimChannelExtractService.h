// GenericSimChannelExtractService.h

// David Adams
// December 2015
//
// Interface for a service that extracts charge from
// a SimChannel object and assigns it to ticks.
//
// The charge is distributed over two arrays: sig and xsig.
// The first is for normal collection/induction. The second
// is for charge that is collected on the wire even if it is
// an induction plane.

#ifndef GenericSimChannelExtractService_H
#define GenericSimChannelExtractService_H

#include <vector>
#include "dune/DuneInterface/SimChannelExtractService.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/LArFFT.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"

namespace sim {
class SimChannel;
}

class GenericSimChannelExtractService : public SimChannelExtractService {

public:

  GenericSimChannelExtractService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int extract(const sim::SimChannel* psc, AdcSignalVector& sig) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  art::ServiceHandle<util::LArFFT> m_pfft;
  art::ServiceHandle<util::SignalShapingServiceDUNE> m_psss;
  unsigned int m_ntick;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(GenericSimChannelExtractService, SimChannelExtractService, LEGACY)

#endif

