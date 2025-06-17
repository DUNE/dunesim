#ifndef GAPCHARGETRANSPORT_HH
#define GAPCHARGETRANSPORT_HH

#include "GapInfo.hh"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/IonizationScintillation/IGapChargeTransport.hh"

namespace gap {
  class GapChargeTransport : public gap::IGapChargeTransport {
  public:
    explicit GapChargeTransport(fhicl::ParameterSet const& pset);

    std::pair<geo::Point_t, int> GetOffset(double x, double y, double z, int n) const override;

    bool EnableGapChargeTransport() const override { return fEnableGapChargeTransport; }

  private:
    bool fEnableGapChargeTransport;
    art::ServiceHandle<geo::Geometry> geom;
    std::string functionType;
    double p1;
    double p2;
    double p3;
    double activeVolumeSizeX, activeVolumeSizeY, activeVolumeSizeZ;
    double average_gap_X, average_gap_Y, average_gap_Z;
    double minZ, minY, minX, maxX, maxY, maxZ;
    int nTPCsX, nTPCsY, nTPCsZ, totalTPCs;
    double detector_size_X, detector_size_Y, detector_size_Z;
    double DistanceToNearestActiveVolume(GapInfo gap, double x) const;
    double ComputeShiftProbability(GapInfo gap, double x, double y, double z) const;
    double LinearShiftProbability(GapInfo gap, double x, double y, double z) const;
    GapInfo Gap(double x, double y, double z) const;
  };
}
#endif
