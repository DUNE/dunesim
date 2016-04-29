// Dune35tSimChannelExtractService.h

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

#ifndef Dune35tSimChannelExtractService_H
#define Dune35tSimChannelExtractService_H

#include <vector>
#include "dune/DuneInterface/SimChannelExtractService.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/LArFFT.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"

namespace sim {
class SimChannel;
}

class Dune35tSimChannelExtractService : public SimChannelExtractService {

public:

  Dune35tSimChannelExtractService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int extract(const sim::SimChannel* psc, AdcSignalVector& sig) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  typedef enum {
    NONACTIVE, UCOMB, VCOMB, ACTIVE, HORIZGAP, VERTGAP
  } GapType_t;

  void init();
  GapType_t combtest35t(double x, double y, double z) const;

  int GapHasDeflector(double x, double y, double z) const;

  art::ServiceHandle<util::LArFFT> m_pfft;
  art::ServiceHandle<util::SignalShapingServiceDUNE> m_psss;
  unsigned int m_ntick;

  unsigned int fFirstCollectionChannel;  // 1st coll channel is used for shaping extra charge

  bool m_init = false;

  std::vector<float> fFractUUCollect;    // fraction of charge that collects on U (non-transparency) when charge drifts over the comb holding U wires
  std::vector<float> fFractUVCollect;    // fraction of charge that collects on U (non-transparency) when charge drifts over the comb holding V wires
  std::vector<float> fFractVUCollect;    // fraction of charge that collects on V (non-transparency) when charge drifts over the comb holding U wires
  std::vector<float> fFractVVCollect;    // fraction of charge that collects on V (non-transparency) when charge drifts over the comb holding V wires
  std::vector<float> fFractUUMiss;       // fraction of charge that gets missed on U when charge drifts over the comb holding U
  std::vector<float> fFractUVMiss;       // fraction of charge that gets missed on U when charge drifts over the comb holding V
  std::vector<float> fFractVUMiss;       // fraction of charge that gets missed on V when charge drifts over the comb holding U
  std::vector<float> fFractVVMiss;       // fraction of charge that gets missed on V when charge drifts over the comb holding V
  std::vector<float> fFractZUMiss;       // fraction of charge that gets missed on Z (collection)  when charge drifts over the comb holding U
  std::vector<float> fFractZVMiss;       // fraction of charge that gets missed on Z (collection)  when charge drifts over the comb holding V
  std::vector<float> fFractHorizGapUMiss;     // fraction of charge in the horizontal gap that is missing on U (and not collected)
  std::vector<float> fFractVertGapUMiss;     // fraction of charge in the horizontal gaps that is missing on U
  std::vector<float> fFractHorizGapVMiss;     // fraction of charge in the horizontal gap that is missing on V
  std::vector<float> fFractVertGapVMiss;     // fraction of charge in the horizontal gaps that is missing on V
  std::vector<float> fFractHorizGapZMiss;     // fraction of charge in the horizontal gap that is missing on Z (collection)
  std::vector<float> fFractVertGapZMiss;     // fraction of charge in the horizontal gaps that is missing on Z (collection
  std::vector<float> fFractHorizGapUCollect;     // fraction of charge in the horizontal gap that collects on U
  std::vector<float> fFractVertGapUCollect;     // fraction of charge in the horizontal gaps that collects on U
  std::vector<float> fFractHorizGapVCollect;     // fraction of charge in the horizontal gap that collects on V
  std::vector<float> fFractVertGapVCollect;     // fraction of charge in the horizontal gaps that collects on V

  // boundaries of the combs -- cached here for speed
  double zcomb1,zcomb2,zcomb3,zcomb4,zcomb5,zcomb6;
  double zcomb7,zcomb8,zcomb9,zcomb10,zcomb11,zcomb12;
  double zcomb13,zcomb14,zcomb15,zcomb16,zcomb17,zcomb18;
  double ycomb1,ycomb2,ycomb3,ycomb4,ycomb5,ycomb6;
  double ycomb7,ycomb8,ycomb9,ycomb10,ycomb11,ycomb12;
  double ycomb13,ycomb14,ycomb15,ycomb16,ycomb17,ycomb18;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(Dune35tSimChannelExtractService, SimChannelExtractService, LEGACY)

#endif

