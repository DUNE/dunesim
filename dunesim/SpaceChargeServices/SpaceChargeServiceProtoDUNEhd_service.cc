////////////////////////////////////////////////////////////////////////
// \file SpaceChargeProtoDUNEhd.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for ProtoDUNEhd
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "dunesim/SpaceChargeServices/SpaceChargeServiceProtoDUNEhd.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ROOT includes
#include "TMath.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
spacecharge::SpaceChargeServiceProtoDUNEhd::SpaceChargeServiceProtoDUNEhd(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new spacecharge::SpaceChargeProtoDUNEhd(pset));
  
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  fProp->Configure(pset,detProp);

  reg.sPreBeginRun.watch(this, &SpaceChargeServiceProtoDUNEhd::preBeginRun);
}

//----------------------------------------------
void spacecharge::SpaceChargeServiceProtoDUNEhd::preBeginRun(const art::Run& run)
{
  fProp->Update(run.id().run());
}

//------------------------------------------------
void spacecharge::SpaceChargeServiceProtoDUNEhd::reconfigure(fhicl::ParameterSet const& pset)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  fProp->Configure(pset,detProp);
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceProtoDUNEhd, spacecharge::SpaceChargeService)
