////////////////////////////////////////////////////////////////////////
// \file SpaceCharge3x1x1dphase.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for the 3x1x1 detector.
// \Adapted from SpaceChargeServiceProtoDUNE_service.cc
//
// \author kevin.fusshoeller@cern.ch
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

// LArSoft includes
#include "dune/SpaceChargeServices/SpaceChargeService3x1x1dphase.h"

// ROOT includes
#include "TMath.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
spacecharge::SpaceChargeService3x1x1dphase::SpaceChargeService3x1x1dphase(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new spacecharge::SpaceCharge3x1x1dphase(pset));

  reg.sPreBeginRun.watch(this, &SpaceChargeService3x1x1dphase::preBeginRun);
}

//----------------------------------------------
void spacecharge::SpaceChargeService3x1x1dphase::preBeginRun(const art::Run& run)
{
  fProp->Update(run.id().run());
}

//------------------------------------------------
void spacecharge::SpaceChargeService3x1x1dphase::reconfigure(fhicl::ParameterSet const& pset)
{
  fProp->Configure(pset);  
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeService3x1x1dphase, spacecharge::SpaceChargeService)
