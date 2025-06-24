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

    std::string fp4_volume;
    std::string Volume() const override;

  private:
    enum class FunctionType { Linear };
    FunctionType ParseFunctionType(const std::string& s) const;
    art::ServiceHandle<geo::Geometry> geom;
    FunctionType fFunctionType;
    double fp1_dist, fp2_prob, fp3_shift;
    double fActiveVolumeSizeX, fActiveVolumeSizeY, fActiveVolumeSizeZ;
    double fAverage_gap_X, fAverage_gap_Y, fAverage_gap_Z;
    double fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ;
    int fnTPCsX, fnTPCsY, fnTPCsZ, fTotalTPCs;
    double fDetector_size_X, fDetector_size_Y, fDetector_size_Z;
    double DistanceToNearestActiveVolume(GapInfo gap, double x) const;
    double ComputeShiftProbability(GapInfo gap, double x, double y, double z) const;
    double LinearShiftProbability(GapInfo gap, double x, double y, double z) const;
    GapInfo Gap(double x, double y, double z) const;
  };
}
#endif
