// DPhaseFixedDetPedestalService.h
//
// R.Sulej: adopted to kY and kZ views
//
// Implementation of lariov::DetPedestalService and its provider lariov::DetPedestal
// that returns fixed pedestals depending on the readout strip orientation.

#ifndef DPhaseFixedDetPedestalService_H
#define DPhaseFixedDetPedestalService_H

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

class DPhaseFixedDetPedestalService : public lariov::DetPedestalService, public lariov::DetPedestalProvider {

public:

  // Service ctor.
  DPhaseFixedDetPedestalService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);

  // Provider ctor.
  DPhaseFixedDetPedestalService(const fhicl::ParameterSet& pset);

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
  float m_PedMeanX;
  float m_PedMeanY;
  float m_PedMeanZ;
  float m_PedRmsX;
  float m_PedRmsY;
  float m_PedRmsZ;
  float m_PedMeanErrX;
  float m_PedMeanErrY;
  float m_PedMeanErrZ;
  float m_PedRmsErrX;
  float m_PedRmsErrY;
  float m_PedRmsErrZ;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(DPhaseFixedDetPedestalService, lariov::DetPedestalService, LEGACY)

#endif
