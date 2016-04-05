// LarsoftHuffmanCompressService.h
//
// David Adams
// February 2016
//
// Service to compress an ADC vector using the LArSoft
// conventions in larsoft/lardata/RawData/raw.cxx.
// Uncompress in that same utility can be used to uncompress.
// Data may be put in block format, Huffman encoded or both or neither.
// In the latter case, any zero-supressed values are replaced with zero.
//
// FCL parameters:
//    UseBlock - Put the data in larsoft block format
//    UseHuffman - Do Huffman encoding
//    LogLevel - Log messaging level
//               0 - no messages
//               1 - Initialization only
//               2 - Short message for each event
//               3 - Long message for each event
//               4 - Very long message for each event

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
  LarsoftHuffmanCompressService(bool useBlock, bool useHuffman, int logLevel);

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

  // Convert to block format.
  void block(const AdcCountVector& oldsigs, const AdcFilterVector& keep, AdcCountVector& newsigs) const;

private:

  bool m_UseBlock;
  bool m_UseHuffman;
  int m_LogLevel;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(LarsoftHuffmanCompressService, AdcCompressService, LEGACY)

#endif
