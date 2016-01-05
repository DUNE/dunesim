// ReplaceCompressService.h
//
// David Adams
// November 2015
//
// Utility/service to compress an ADC vector by replacing filtered out values
// with a configurable zero value plus a passed offset (pedestal).
//
// FCL parameters
//   Zero [0] - replacement value

#ifndef ReplaceCompressService_H
#define ReplaceCompressService_H

#include "dune/DuneInterface/AdcCompressService.h"

namespace fhicl {
class ParameterSet;
}
namespace art {
class ActivityRegistry;
}

class ReplaceCompressService : public AdcCompressService {

public:

  // Ctor from parameters that characterize the algorithm.
  ReplaceCompressService(AdcCount azero =0);

  // Ctor from fcl.
  ReplaceCompressService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Compress a vector of signals.
  // Suppressed signals are replaced with the value of offset + (FCL parameter) Zero.
  int compress(AdcCountVector& sigs,
               const AdcFilterVector& keep,
               AdcCount offset,
               raw::Compress_t& comp) const;

  // Return the value assigned to suppressed channels.
  AdcCount zero() const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="  ") const;

private:

  AdcCount m_zero;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(ReplaceCompressService, AdcCompressService, LEGACY)

#endif
