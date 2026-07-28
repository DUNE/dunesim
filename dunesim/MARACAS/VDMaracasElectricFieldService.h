////////////////////////////////////////////////////////////////////////
// \file VDMaracasElectricFieldService.h
//
// \brief ElectricFieldService serving a VDMaracasElectricField provider.
////////////////////////////////////////////////////////////////////////
#ifndef MARACAS_VDMARACASELECTRICFIELDSERVICE_H
#define MARACAS_VDMARACASELECTRICFIELDSERVICE_H

#include "dunesim/MARACAS/VDMaracasElectricField.h"
#include "lardata/DetectorInfoServices/ElectricFieldService.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

namespace art {
  class ActivityRegistry;
}
namespace fhicl {
  class ParameterSet;
}

namespace maracas {
  class VDMaracasElectricFieldService : public detinfo::ElectricFieldService {
  public:
    VDMaracasElectricFieldService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
      : fProp{pset}
    {}

  private:
    detinfo::IElectricFieldProvider const* provider() const override { return &fProp; }

    maracas::VDMaracasElectricField fProp;
  };
} // namespace maracas

DECLARE_ART_SERVICE_INTERFACE_IMPL(maracas::VDMaracasElectricFieldService,
                                   detinfo::ElectricFieldService,
                                   SHARED)

#endif // MARACAS_VDMARACASELECTRICFIELDSERVICE_H
