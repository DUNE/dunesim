// FixedDetPedestalService.h
//
// Implementation of lariov::DetPedestalService and its provider lariov::DetPedestal
// that returns fixed pedestals depending on wire orientation.

#ifndef FixedDetPedestalService_H
#define FixedDetPedestalService_H

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larcore/Geometry/Geometry.h"

namespace fhicl {
class ParameterSet;
}
namespace art {
class ActivityRegistry;
}

class FixedDetPedestalService : public lariov::DetPedestalService, public lariov::DetPedestalProvider {

public:

  // Service ctor.
  FixedDetPedestalService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);

  // Provider ctor.
  FixedDetPedestalService(const fhicl::ParameterSet& pset);

  // Retrieve pedestal information (provider interface).    
  float PedMean(raw::ChannelID_t ch) const;
  float PedRms(raw::ChannelID_t ch) const;
  float PedMeanErr(raw::ChannelID_t ch) const;
  float PedRmsErr(raw::ChannelID_t ch) const;

private:

  // Return provider (service interface).
  const lariov::DetPedestalProvider& DoGetPedestalProvider() const;

  // We use the geometry service to obtain the orientation for each channel.
  art::ServiceHandle<geo::Geometry> m_hgeo;

  // Pedestal value.
  float m_PedMeanU;
  float m_PedMeanV;
  float m_PedMeanZ;
  float m_PedRmsU;
  float m_PedRmsV;
  float m_PedRmsZ;
  float m_PedMeanErrU;
  float m_PedMeanErrV;
  float m_PedMeanErrZ;
  float m_PedRmsErrU;
  float m_PedRmsErrV;
  float m_PedRmsErrZ;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(FixedDetPedestalService, lariov::DetPedestalService, LEGACY)

#endif
