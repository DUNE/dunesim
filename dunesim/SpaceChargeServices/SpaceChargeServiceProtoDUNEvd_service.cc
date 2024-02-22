////////////////////////////////////////////////////////////////////////
// \file SpaceChargeProtoDUNEvd.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for ProtoDUNEvd
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "dunesim/SpaceChargeServices/SpaceChargeServiceProtoDUNEvd.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ROOT includes
#include "TMath.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
spacecharge::SpaceChargeServiceProtoDUNEvd::SpaceChargeServiceProtoDUNEvd(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new spacecharge::SpaceChargeProtoDUNEvd(pset));
  
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  fProp->Configure(pset,detProp);

  reg.sPreBeginRun.watch(this, &SpaceChargeServiceProtoDUNEvd::preBeginRun);
}

//----------------------------------------------
void spacecharge::SpaceChargeServiceProtoDUNEvd::preBeginRun(const art::Run& run)
{
  fProp->Update(run.id().run());
}

//------------------------------------------------
void spacecharge::SpaceChargeServiceProtoDUNEvd::reconfigure(fhicl::ParameterSet const& pset)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  fProp->Configure(pset,detProp);
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceProtoDUNEvd, spacecharge::SpaceChargeService)
