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

// ROOT includes
#include "TMath.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
spacecharge::SpaceChargeServiceProtoDUNE::SpaceChargeServiceProtoDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new spacecharge::SpaceChargeProtoDUNE(pset));

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
  fProp->Configure(pset);  
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceProtoDUNE, spacecharge::SpaceChargeService)
