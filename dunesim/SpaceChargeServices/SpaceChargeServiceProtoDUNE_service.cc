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
#include "dune/SpaceChargeServices/SpaceChargeServiceProtoDUNE.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ROOT includes
#include "TMath.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
spacecharge::SpaceChargeServiceProtoDUNE::SpaceChargeServiceProtoDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new spacecharge::SpaceChargeProtoDUNE(pset));
  
  auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fProp->Configure(pset,detprop);

  reg.sPreBeginRun.watch(this, &SpaceChargeServiceProtoDUNE::preBeginRun);
}

//----------------------------------------------
void spacecharge::SpaceChargeServiceProtoDUNE::preBeginRun(const art::Run& run)
{
  fProp->Update(run.id().run());
}

//------------------------------------------------
void spacecharge::SpaceChargeServiceProtoDUNE::reconfigure(fhicl::ParameterSet const& pset)
{
  auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();  
  fProp->Configure(pset,detprop);  
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceProtoDUNE, spacecharge::SpaceChargeService)
