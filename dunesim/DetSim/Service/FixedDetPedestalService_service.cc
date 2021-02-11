// FixedDetPedestalService_service.cc

#include "dune/DetSim/Service/FixedDetPedestalService.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"

//**********************************************************************

FixedDetPedestalService::
FixedDetPedestalService(const fhicl::ParameterSet& pset, art::ActivityRegistry&)
: FixedDetPedestalService(pset.get<fhicl::ParameterSet>("FixedDetPedestal")) { }

//**********************************************************************

FixedDetPedestalService::
FixedDetPedestalService(const fhicl::ParameterSet& pset) {
  m_PedMeanU    = pset.get<float>("PedMeanU");
  m_PedMeanV    = pset.get<float>("PedMeanV");
  m_PedMeanX    = pset.get<float>("PedMeanX");
  m_PedMeanY    = pset.get<float>("PedMeanY");
  m_PedMeanZ    = pset.get<float>("PedMeanZ");
  m_PedRmsU     = pset.get<float>("PedRmsU");
  m_PedRmsV     = pset.get<float>("PedRmsV");
  m_PedRmsX     = pset.get<float>("PedRmsX");
  m_PedRmsY     = pset.get<float>("PedRmsY");
  m_PedRmsZ     = pset.get<float>("PedRmsZ");
  m_PedMeanErrU = pset.get<float>("PedMeanErrU");
  m_PedMeanErrV = pset.get<float>("PedMeanErrV");
  m_PedMeanErrX = pset.get<float>("PedMeanErrX");
  m_PedMeanErrY = pset.get<float>("PedMeanErrY");
  m_PedMeanErrZ = pset.get<float>("PedMeanErrZ");
  m_PedRmsErrU  = pset.get<float>("PedRmsErrU");
  m_PedRmsErrV  = pset.get<float>("PedRmsErrV");
  m_PedRmsErrX  = pset.get<float>("PedRmsErrX");
  m_PedRmsErrY  = pset.get<float>("PedRmsErrY");
  m_PedRmsErrZ  = pset.get<float>("PedRmsErrZ");
}

//**********************************************************************

float FixedDetPedestalService::PedMean(raw::ChannelID_t ch) const {
  if      ( m_hgeo->View(ch) == geo::kX ) return m_PedMeanX;
  else if ( m_hgeo->View(ch) == geo::kY ) return m_PedMeanY;
  else if ( m_hgeo->View(ch) == geo::kZ ) return m_PedMeanZ;
  else if ( m_hgeo->View(ch) == geo::kU ) return m_PedMeanU;
  else if ( m_hgeo->View(ch) == geo::kV ) return m_PedMeanV;
  return -999.0;
}

//**********************************************************************

float FixedDetPedestalService::PedRms(raw::ChannelID_t ch) const {
  return m_PedRmsZ;
}

//**********************************************************************

float FixedDetPedestalService::PedMeanErr(raw::ChannelID_t ch) const {
  return m_PedMeanErrZ;
}

//**********************************************************************

float FixedDetPedestalService::PedRmsErr(raw::ChannelID_t ch) const {
  return m_PedRmsErrZ;
}

//**********************************************************************

const lariov::DetPedestalProvider&
FixedDetPedestalService::DoGetPedestalProvider() const {
  return *this;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(FixedDetPedestalService, lariov::DetPedestalService)

//**********************************************************************
