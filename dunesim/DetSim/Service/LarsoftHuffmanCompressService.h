// LarsoftHuffmanCompressService.h
//
// David Adams
// November 2015
//
// Utility/service to compress an ADC vector using the LArSoft
// encoding in larsoft/lardata/RawData/raw.cxx.
// Uncmpress in that same utility and be used to uncompress.
//
// No FCL parameters

#ifndef LarsoftHuffmanCompressService_H
#define LarsoftHuffmanCompressService_H

#include "dune/DuneInterface/AdcCompressService.h"

namespace fhicl {
class ParameterSet;
}
namespace art {
class ActivityRegistry;
}

class LarsoftHuffmanCompressService : public AdcCompressService {

public:

  // Ctor from parameters that characterize the algorithm.
  LarsoftHuffmanCompressService();

  // Ctor from fcl.
  LarsoftHuffmanCompressService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Compress a vector of signals.
  // Suppressed signals are replaced with the value of offset + (FCL parameter) Zero.
  int compress(AdcCountVector& sigs,
               const AdcFilterVector& keep,
               AdcCount offset,
               raw::Compress_t& comp) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="  ") const;

private:

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(LarsoftHuffmanCompressService, AdcCompressService, LEGACY)

#endif
