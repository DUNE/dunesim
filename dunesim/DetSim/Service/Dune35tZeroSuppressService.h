// Dune35tZeroSuppressService.h
//
// David Adams
// November 2015
//
// Utility/service wrapper for 35-ton zero suppression.
//
// The algorithm is that described by James Russell at the September 2015
// DUNE collaboration meeting:
//   https://indico.fnal.gov/conferenceOtherViews.py?view=standard&confId=10100
// and is described with later modifications at
//   https://cdcvs.fnal.gov/redmine/projects/35ton/wiki/Data_compression_and_zero_suppression

#ifndef Dune35tZeroSuppressService_H
#define Dune35tZeroSuppressService_H

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

class Dune35tZeroSuppressService : public AdcSuppressService {

public:

  typedef unsigned int Index;

  // Ctor from parameters.
  Dune35tZeroSuppressService(AdcCount ts, AdcCount tl, AdcCount td, Index ns, Index nl, Index nd, Index nt);

  // Ctor from FCL.
  Dune35tZeroSuppressService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Filter an array of signals. Result is written to keep.
  int filter(const AdcCountVector& sigs, Channel chan, AdcPedestal ped, AdcFilterVector& keep) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

  // Set the debug flag.
  void setDebug(int dbg);

private:

  // Parameters.
  AdcCount m_ts;
  AdcCount m_tl;
  AdcCount m_td;
  Index m_ns;
  Index m_nl;
  Index m_nd;
  Index m_nt;
  int m_LogLevel;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(Dune35tZeroSuppressService, AdcSuppressService, LEGACY)

#endif
