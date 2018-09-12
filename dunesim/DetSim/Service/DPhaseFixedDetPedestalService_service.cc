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
  m_PedMeanX    = pset.get<float>("PedMeanX");
  m_PedMeanY    = pset.get<float>("PedMeanY");
  m_PedMeanZ    = pset.get<float>("PedMeanZ");
  m_PedRmsX     = pset.get<float>("PedRmsX");
  m_PedRmsY     = pset.get<float>("PedRmsY");
  m_PedRmsZ     = pset.get<float>("PedRmsZ");
  m_PedMeanErrX = pset.get<float>("PedMeanErrX");
  m_PedMeanErrY = pset.get<float>("PedMeanErrY");
  m_PedMeanErrZ = pset.get<float>("PedMeanErrZ");
  m_PedRmsErrX  = pset.get<float>("PedRmsErrX");
  m_PedRmsErrY  = pset.get<float>("PedRmsErrY");
  m_PedRmsErrZ  = pset.get<float>("PedRmsErrZ");
}

//**********************************************************************

float DPhaseFixedDetPedestalService::PedMean(raw::ChannelID_t ch) const {
  if      ( m_hgeo->View(ch) == geo::kX ) return m_PedMeanX;
  else if ( m_hgeo->View(ch) == geo::kY ) return m_PedMeanY;
  else if ( m_hgeo->View(ch) == geo::kZ ) return m_PedMeanZ;
  return -999.0;
}

//**********************************************************************

float DPhaseFixedDetPedestalService::PedRms(raw::ChannelID_t ch) const {
  return m_PedRmsZ;
}

//**********************************************************************

float DPhaseFixedDetPedestalService::PedMeanErr(raw::ChannelID_t ch) const {
  if      ( m_hgeo->View(ch) == geo::kX ) return m_PedMeanErrX;
  else if ( m_hgeo->View(ch) == geo::kY ) return m_PedMeanErrY;
  else if ( m_hgeo->View(ch) == geo::kZ ) return m_PedMeanErrZ;
  return -999.0;
}

//**********************************************************************

float DPhaseFixedDetPedestalService::PedRmsErr(raw::ChannelID_t ch) const {
  if      ( m_hgeo->View(ch) == geo::kX ) return m_PedRmsErrX;
  else if ( m_hgeo->View(ch) == geo::kY ) return m_PedRmsErrY;
  else if ( m_hgeo->View(ch) == geo::kZ ) return m_PedRmsErrZ;
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
