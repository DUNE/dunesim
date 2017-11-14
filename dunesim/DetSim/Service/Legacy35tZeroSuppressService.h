// Legacy35tZeroSuppressService.h
//
// David Adams
// November 2015
//
// Service wrapper for legacy 35-ton zero suppression.
// Code taken from run.cxx in dunetpc v04.29.02.
// Parameters:
//   AdcThreshold - threshold
//   TickRange - range
//   SuppressStickyBits - treat sticky bits as below threshold
// If a tick I has pedestal-corrected |ADC| > AdcThreshold, ticks
// [i-TickRange, i+TickRange] are retained.
// If tick I has sticky bits and has pedestal-corrected ADC < 64,
// it is treated as though its pedstal-corrected ADC is zero.
// If the gap between unsuppressed blocks is less that MinTickGap,
// then the bins in that gap are not suppressed. Edges (bins
// outside the tick range are treated as unsupressed.

#ifndef Legacy35tZeroSuppressService_H
#define Legacy35tZeroSuppressService_H

#include "dune/DuneInterface/AdcSuppressService.h"

#include <memory>
#include <string>
#include <iostream>

namespace fhicl {
class ParameterSet;
}
namespace art {
class ActivityRegistry;
}

class Legacy35tZeroSuppressService : public AdcSuppressService {

public:

  // Ctor from fcl.
  Legacy35tZeroSuppressService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Ctor from direct params.
  Legacy35tZeroSuppressService(float aAdcThreshold,
                               unsigned int aTickRange,
                               unsigned int aMinTickGap,
                               bool aSuppressStickyBits);

  // Filter an array of signals. Result is written to keep.
  int filter(const AdcCountVector& sigs, Channel chan, AdcPedestal ped, AdcFilterVector& keep) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;


private:

  // Parameters.
  float         m_AdcThreshold;
  unsigned int  m_TickRange;
  unsigned int  m_MinTickGap;
  bool          m_SuppressStickyBits;


};

DECLARE_ART_SERVICE_INTERFACE_IMPL(Legacy35tZeroSuppressService, AdcSuppressService, LEGACY)

#endif
