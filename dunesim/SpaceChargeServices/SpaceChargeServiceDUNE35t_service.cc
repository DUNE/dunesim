////////////////////////////////////////////////////////////////////////
// \file SpaceChargeDUNE35t.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for DUNE 35ton
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

// LArSoft includes
#include "dune/SpaceChargeServices/SpaceChargeServiceDUNE35t.h"

// ROOT includes
#include "TMath.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
spacecharge::SpaceChargeServiceDUNE35t::SpaceChargeServiceDUNE35t(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new spacecharge::SpaceChargeDUNE35t(pset));

  reg.sPreBeginRun.watch(this, &SpaceChargeServiceDUNE35t::preBeginRun);
}

//----------------------------------------------
void spacecharge::SpaceChargeServiceDUNE35t::preBeginRun(const art::Run& run)
{
  fProp->Update(run.id().run());
}

//------------------------------------------------
void spacecharge::SpaceChargeServiceDUNE35t::reconfigure(fhicl::ParameterSet const& pset)
{
  fProp->Configure(pset);  
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceDUNE35t, spacecharge::SpaceChargeService)
