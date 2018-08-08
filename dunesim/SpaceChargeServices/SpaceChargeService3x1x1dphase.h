////////////////////////////////////////////////////////////////////////
// \file SpaceChargeServiceProtoDUNE.h
//
// \brief header of service for storing/accessing space charge distortions for the 3x1x1 detector.
// \Adapted from SpaceChargeServiceProtoDUNE.h
//
// \author kevin.fusshoeller@cern.ch
// 
////////////////////////////////////////////////////////////////////////
#ifndef SPACECHARGESERVICE3X1X1DPHASE_H
#define SPACECHARGESERVICE3X1X1DPHASE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "dune/SpaceCharge/SpaceCharge3x1x1dphase.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


namespace spacecharge{
  class SpaceChargeService3x1x1dphase : public SpaceChargeService {
    public:
      
      // this enables art to print the configuration help:
      //using Parameters = art::ServiceTable<spacecharge::SpaceCharge3x1x1dphase::ConfigurationParameters_t>;
      
      SpaceChargeService3x1x1dphase(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

      virtual void   reconfigure(fhicl::ParameterSet const& pset)  override;
      void   preBeginRun(const art::Run& run);


      virtual const  provider_type* provider() const override { return fProp.get();}

    private:

      std::unique_ptr<spacecharge::SpaceCharge3x1x1dphase> fProp;

    }; // class SpaceChargeService3x1x1dphase
} //namespace spacecharge
DECLARE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeService3x1x1dphase, spacecharge::SpaceChargeService, LEGACY)
#endif // SPACECHARGESERVICE3X1X1DPHASE_H
