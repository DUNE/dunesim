// DPhaseFixedDetPedestalService_service.cc

#include "dune/DetSim/Service/DPhaseFixedDetPedestalService.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"

//**********************************************************************

DPhaseFixedDetPedestalService::
DPhaseFixedDetPedestalService(const fhicl::ParameterSet& pset, art::ActivityRegistry&)
: DPhaseFixedDetPedestalService(pset.get<fhicl::ParameterSet>("FixedDetPedestal")) { }

//**********************************************************************

DPhaseFixedDetPedestalService::
DPhaseFixedDetPedestalService(const fhicl::ParameterSet& pset) {
  m_PedMeanY    = pset.get<float>("PedMeanY");
  m_PedMeanZ    = pset.get<float>("PedMeanZ");
  m_PedRmsY     = pset.get<float>("PedRmsY");
  m_PedRmsZ     = pset.get<float>("PedRmsZ");
  m_PedMeanErrY = pset.get<float>("PedMeanErrY");
  m_PedMeanErrZ = pset.get<float>("PedMeanErrZ");
  m_PedRmsErrY  = pset.get<float>("PedRmsErrY");
  m_PedRmsErrZ  = pset.get<float>("PedRmsErrZ");
}

//**********************************************************************

float DPhaseFixedDetPedestalService::PedMean(raw::ChannelID_t ch) const {
  if      ( m_hgeo->View(ch) == geo::kZ ) return m_PedMeanZ;
  else if ( m_hgeo->View(ch) == geo::kY ) return m_PedMeanY;
  return -999.0;
}

//**********************************************************************

float DPhaseFixedDetPedestalService::PedRms(raw::ChannelID_t ch) const {
  return m_PedRmsZ;
}

//**********************************************************************

float DPhaseFixedDetPedestalService::PedMeanErr(raw::ChannelID_t ch) const {
  return m_PedMeanErrZ;
}

//**********************************************************************

float DPhaseFixedDetPedestalService::PedRmsErr(raw::ChannelID_t ch) const {
  return m_PedRmsErrZ;
}

//**********************************************************************

const lariov::DetPedestalProvider&
DPhaseFixedDetPedestalService::DoGetPedestalProvider() const {
  return *this;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DPhaseFixedDetPedestalService, lariov::DetPedestalService)

//**********************************************************************
