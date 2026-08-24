////////////////////////////////////////////////////////////////////////
// \file CSUEFieldServiceProtoDUNE.h
//
// \brief ElectricFieldService serving a spacecharge::CSUEFieldProtoDUNE
//        provider (CSU ProtoDUNE space-charge E-field maps).
////////////////////////////////////////////////////////////////////////
#ifndef DUNESIM_SPACECHARGESERVICES_CSUEFIELDSERVICEPROTODUNE_H
#define DUNESIM_SPACECHARGESERVICES_CSUEFIELDSERVICEPROTODUNE_H

#include "dunesim/SpaceCharge/CSUEFieldProtoDUNE.h"
#include "lardata/DetectorInfoServices/ElectricFieldService.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

namespace art {
  class ActivityRegistry;
}
namespace fhicl {
  class ParameterSet;
}

namespace spacecharge {
  class CSUEFieldServiceProtoDUNE : public detinfo::ElectricFieldService {
  public:
    CSUEFieldServiceProtoDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
      : fProp{pset}
    {}

  private:
    detinfo::IElectricFieldProvider const* provider() const override { return &fProp; }

    spacecharge::CSUEFieldProtoDUNE fProp;
  };
} // namespace spacecharge

DECLARE_ART_SERVICE_INTERFACE_IMPL(spacecharge::CSUEFieldServiceProtoDUNE,
                                   detinfo::ElectricFieldService,
                                   SHARED)

#endif // DUNESIM_SPACECHARGESERVICES_CSUEFIELDSERVICEPROTODUNE_H
