////////////////////////////////////////////////////////////////////////
// \file SpaceChargeServiceProtoDUNEdp.h
//
// \brief header of service for storing/accessing space charge distortions for ProtoDUNE
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef SPACECHARGESERVICEPROTODUNEDP_H
#define SPACECHARGESERVICEPROTODUNEDP_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "dune/SpaceCharge/SpaceChargeProtoDUNEdp.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


namespace spacecharge{
  class SpaceChargeServiceProtoDUNEdp : public SpaceChargeService {
    public:
      
      // this enables art to print the configuration help:
      //using Parameters = art::ServiceTable<spacecharge::SpaceChargeProtoDUNE::ConfigurationParameters_t>;
      
      SpaceChargeServiceProtoDUNEdp(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

      virtual void   reconfigure(fhicl::ParameterSet const& pset);
      void   preBeginRun(const art::Run& run);


      virtual const  provider_type* provider() const override { return fProp.get();}

    private:

      std::unique_ptr<spacecharge::SpaceChargeProtoDUNEdp> fProp;

    }; // class SpaceChargeServiceProtoDUNE
} //namespace spacecharge
DECLARE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceProtoDUNEdp, spacecharge::SpaceChargeService, LEGACY)
#endif // SPACECHARGESERVICEPROTODUNE_H
