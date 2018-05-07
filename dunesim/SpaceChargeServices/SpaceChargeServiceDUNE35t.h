////////////////////////////////////////////////////////////////////////
// \file SpaceChargeServiceDUNE35t.h
//
// \brief header of service for storing/accessing space charge distortions for DUNE 35ton
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef SPACECHARGESERVICEDUNE35T_H
#define SPACECHARGESERVICEDUNE35T_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "dune/SpaceCharge/SpaceChargeDUNE35t.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


namespace spacecharge{
  class SpaceChargeServiceDUNE35t : public SpaceChargeService {
    public:
      
      // this enables art to print the configuration help:
      //using Parameters = art::ServiceTable<spacecharge::SpaceChargeDUNE35t::ConfigurationParameters_t>;
      
      SpaceChargeServiceDUNE35t(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

      virtual void   reconfigure(fhicl::ParameterSet const& pset) override;
      void   preBeginRun(const art::Run& run);


      virtual const  provider_type* provider() const override { return fProp.get();}

    private:

      std::unique_ptr<spacecharge::SpaceChargeDUNE35t> fProp;

    }; // class SpaceChargeServiceDUNE35t
} //namespace spacecharge
DECLARE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceDUNE35t, spacecharge::SpaceChargeService, LEGACY)
#endif // SPACECHARGESERVICEDUNE35T_H
