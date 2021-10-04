////////////////////////////////////////////////////////////////////////
// \file SpaceChargeProtoDUNE.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for ProtoDUNE
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "dune/SpaceChargeServices/SpaceChargeServiceProtoDUNEdp.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"


// ROOT includes
#include "TMath.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
spacecharge::SpaceChargeServiceProtoDUNEdp::SpaceChargeServiceProtoDUNEdp(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new spacecharge::SpaceChargeProtoDUNEdp(pset));
  
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  fProp->Configure(pset,detProp);

  reg.sPreBeginRun.watch(this, &SpaceChargeServiceProtoDUNEdp::preBeginRun);
}

//----------------------------------------------
void spacecharge::SpaceChargeServiceProtoDUNEdp::preBeginRun(const art::Run& run)
{
  fProp->Update(run.id().run());
}

//------------------------------------------------
void spacecharge::SpaceChargeServiceProtoDUNEdp::reconfigure(fhicl::ParameterSet const& pset)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  fProp->Configure(pset,detProp);
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceProtoDUNEdp, spacecharge::SpaceChargeService)
