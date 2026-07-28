////////////////////////////////////////////////////////////////////////
// \file VDMaracasPositionDistorterService.h
//
// \brief PositionDistorterService serving a VDMaracasDistorter provider.
////////////////////////////////////////////////////////////////////////
#ifndef DUNESIM_MARACAS_VDMARACASPOSITIONDISTORTERSERVICE_H
#define DUNESIM_MARACAS_VDMARACASPOSITIONDISTORTERSERVICE_H

#include "dunesim/MARACAS/VDMaracasDistorter.h"
#include "lardata/DetectorInfoServices/PositionDistorterService.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

namespace art {
  class ActivityRegistry;
}
namespace fhicl {
  class ParameterSet;
}

namespace maracas {
  class VDMaracasPositionDistorterService : public detinfo::PositionDistorterService {
  public:
    VDMaracasPositionDistorterService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
      : fProp{pset}
    {}

  private:
    detinfo::IPositionDistorter const* provider() const override { return &fProp; }

    maracas::VDMaracasDistorter fProp;
  };
} // namespace maracas

DECLARE_ART_SERVICE_INTERFACE_IMPL(maracas::VDMaracasPositionDistorterService,
                                   detinfo::PositionDistorterService,
                                   SHARED)

#endif // MARACAS_VDMARACASPOSITIONDISTORTERSERVICE_H
