////////////////////////////////////////////////////////////////////////
// \file SpaceChargeServiceProtoDUNEhd.h
//
// \brief header of service for storing/accessing space charge distortions for ProtoDUNEhd
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef SPACECHARGESERVICEPROTODUNE_H
#define SPACECHARGESERVICEPROTODUNE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "dunesim/SpaceCharge/SpaceChargeProtoDUNEhd.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


namespace spacecharge{
  class SpaceChargeServiceProtoDUNEhd : public SpaceChargeService {
    public:
      
      // this enables art to print the configuration help:
      //using Parameters = art::ServiceTable<spacecharge::SpaceChargeProtoDUNEhd::ConfigurationParameters_t>;
      
      SpaceChargeServiceProtoDUNEhd(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

      void   reconfigure(fhicl::ParameterSet const& pset);
      void   preBeginRun(const art::Run& run);


      virtual const  provider_type* provider() const override { return fProp.get();}

    private:

      std::unique_ptr<spacecharge::SpaceChargeProtoDUNEhd> fProp;

    }; // class SpaceChargeServiceProtoDUNEhd
} //namespace spacecharge
DECLARE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceProtoDUNEhd, spacecharge::SpaceChargeService, LEGACY)
#endif // SPACECHARGESERVICEPROTODUNE_H
