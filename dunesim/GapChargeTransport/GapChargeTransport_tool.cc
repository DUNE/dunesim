#include "art/Utilities/ToolMacros.h"

#include "GapChargeTransport.hh"
#include "GapInfo.hh"
#include "larcore/Geometry/Geometry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gap {

  GapChargeTransport::FunctionType GapChargeTransport::ParseFunctionType(const std::string& s) const
  {
    if (s == "linear") return GapChargeTransport::FunctionType::Linear;
    throw cet::exception("GapChargeTransport") << "Unknown FunctionType";
  }

  //constructor
  GapChargeTransport::GapChargeTransport(const fhicl::ParameterSet& pset)
  {
    auto const* geom = lar::providerFrom<geo::Geometry>();

    fFunctionType = ParseFunctionType(pset.get<std::string>("FunctionType"));
    fp1_dist = pset.get<double>("DistanceZero");
    fp2_prob = pset.get<double>("ProbabilityBorder");
    fp3_shift = pset.get<double>("Shift");
    fp4_volume = pset.get<std::string>("Volume");

    const geo::TPCGeo& one_tpc = geom->TPC(geo::TPCID(0, 0));
    geo::BoxBoundedGeo activeVolume = one_tpc.ActiveBoundingBox();

    fActiveVolumeSizeX = activeVolume.MaxX() - activeVolume.MinX();
    fActiveVolumeSizeY = activeVolume.MaxY() - activeVolume.MinY();
    fActiveVolumeSizeZ = activeVolume.MaxZ() - activeVolume.MinZ();

    fTotalTPCs = geom->NTPC();

    fMinX = std::numeric_limits<double>::max();
    fMaxX = std::numeric_limits<double>::lowest();
    fMinY = std::numeric_limits<double>::max();
    fMaxY = std::numeric_limits<double>::lowest();
    fMinZ = std::numeric_limits<double>::max();
    fMaxZ = std::numeric_limits<double>::lowest();

    for (int i = 0; i < fTotalTPCs; ++i) {
      const geo::TPCGeo& tpc = geom->TPC(geo::TPCID(0, i));
      geo::BoxBoundedGeo bbox = tpc.ActiveBoundingBox();
      fMinX = std::min(fMinX, bbox.MinX());
      fMaxX = std::max(fMaxX, bbox.MaxX());
      fMinY = std::min(fMinY, bbox.MinY());
      fMaxY = std::max(fMaxY, bbox.MaxY());
      fMinZ = std::min(fMinZ, bbox.MinZ());
      fMaxZ = std::max(fMaxZ, bbox.MaxZ());
    }

    fDetector_size_X = fMaxX - fMinX;
    fDetector_size_Y = fMaxY - fMinY;
    fDetector_size_Z = fMaxZ - fMinZ;
    fnTPCsX = static_cast<int>(fDetector_size_X / fActiveVolumeSizeX + 0.5);
    fnTPCsY = static_cast<int>(fDetector_size_Y / fActiveVolumeSizeY + 0.5);
    fnTPCsZ = static_cast<int>(fDetector_size_Z / fActiveVolumeSizeZ + 0.5);
    fAverage_gap_Z =
      (fDetector_size_Z - (double)fnTPCsZ * fActiveVolumeSizeZ) / (double)(fnTPCsZ - 1);
    fAverage_gap_Y =
      (fDetector_size_Y - (double)fnTPCsZ * fActiveVolumeSizeY) / (double)(fnTPCsY - 1);
    fAverage_gap_X =
      (fDetector_size_X - (double)fnTPCsX * fActiveVolumeSizeX) / (double)(fnTPCsX - 1);
    if (fAverage_gap_Z > 2 || fAverage_gap_Y > 2 || fAverage_gap_X > 2) {
      mf::LogWarning("GapChargeTransport")
        << "WARNING: The gaps are weirdly huge: y = " << fAverage_gap_Y
        << ", z = " << fAverage_gap_Z << ", x = " << fAverage_gap_X;
    }
  }

  //here the parameters of the gap are calculated
  GapInfo GapChargeTransport::Gap(double x, double y, double z) const
  {
    int id_z = static_cast<int>((z - fMinZ) / (fActiveVolumeSizeZ));
    int id_y = static_cast<int>((y - fMinY) / (fActiveVolumeSizeY));
    int id_x = static_cast<int>((x - fMinX) / (fActiveVolumeSizeX));
    mf::LogInfo("GapChargeTransport")
      << "The bottom left point (min): (" << fMinX << ", " << fMinY << ", " << fMinZ << ")";
    mf::LogInfo("GapChargeTransport")
      << "The active volume dimensions: x = " << fActiveVolumeSizeX
      << ", y = " << fActiveVolumeSizeY << ", z = " << fActiveVolumeSizeZ;
    mf::LogInfo("GapChargeTransport")
      << "id_x = " << id_x << ", id_y = " << id_y << ", id_z = " << id_z;
    int id_total = (id_x) * (fnTPCsY) + (id_y) + (id_z) * (fnTPCsX * fnTPCsY);
    GapInfo::GapType type = GapInfo::GapType::NotAGap;
    GapInfo::MovingDirection direction = GapInfo::MovingDirection::LEFT;
    double left = 0;
    double right = 0;
    if (z < fMinZ || z > fMaxZ || y < fMinY || y > fMaxY || x < fMinX || x > fMaxX ||
        id_total >= fTotalTPCs) {
      type = GapInfo::GapType::NotInDetectorVolume;
    }
    else {
      const geo::TPCGeo& tpc = geom->TPC(geo::TPCID(0, id_total));
      geo::BoxBoundedGeo bbox = tpc.ActiveBoundingBox();
      mf::LogInfo("GapChargeTransport") << "id_total = " << id_total << ", MaxX = " << bbox.MaxX()
                                        << ", MaxY = " << bbox.MaxY() << ", MaxZ = " << bbox.MaxZ();

      if (z < bbox.MinZ()) { //not in the "right" active volume
        left = geom->TPC(geo::TPCID(0, id_total - fnTPCsX * fnTPCsY)).ActiveBoundingBox().MaxZ();
        if (left < z) {
          right = bbox.MinZ();
          type = GapInfo::GapType::Z;
          if (z > (left + right) / 2.) direction = GapInfo::MovingDirection::RIGHT;
        }
      }
      else {
        if (y < bbox.MinY()) {
          left = geom->TPC(geo::TPCID(0, id_total - 1)).ActiveBoundingBox().MaxY();
          if (left < y) {
            type = GapInfo::GapType::Y;
            right = bbox.MinY();
            if (y > (left + right) / 2.) direction = GapInfo::MovingDirection::RIGHT;
          }
        }
        else {
          if (x > bbox.MaxX()) {
            left = geom->TPC(geo::TPCID(0, id_total - fnTPCsY)).ActiveBoundingBox().MaxX();
            if (left < x) {
              right = bbox.MinX();
              type = GapInfo::GapType::X;
              if (x > (left + right) / 2.) direction = GapInfo::MovingDirection::RIGHT;
            }
          }
        }
      }
    }
    mf::LogInfo("GapChargeTransport")
      << "Type: " << GapInfo::GapToString(type) << ", left = " << left << ", right = " << right
      << ", direction = " << GapInfo::GapToString(direction);
    double size = right - left;
    GapInfo the_gap = GapInfo(type, size, left, right, direction);
    return the_gap;
  }

  //in case there will be several modeles, not only linear
  double GapChargeTransport::ComputeShiftProbability(GapInfo the_gap,
                                                     double x,
                                                     double y,
                                                     double z) const
  {
    if (fFunctionType == GapChargeTransport::FunctionType::Linear) {
      return LinearShiftProbability(the_gap, x, y, z);
    }
    else {
      mf::LogWarning("GapChargeTransport") << "Unknown function type, defaulting to linear.";
      return LinearShiftProbability(the_gap, x, y, z);
    }
  }

  //self-explanatory
  double GapChargeTransport::DistanceToNearestActiveVolume(GapInfo the_gap, double x) const
  {
    double distance;
    if (the_gap.GetDirection() == GapInfo::MovingDirection::LEFT) {
      distance = x - the_gap.GetLeftBorder();
    }
    else {
      distance = the_gap.GetRightBorder() - x;
    }
    mf::LogInfo("GapChargeTransport")
      << "left" << the_gap.GetLeftBorder() << " " << x << " " << the_gap.GetRightBorder();
    return distance;
  }

  //computing the probability to move the charge
  double GapChargeTransport::LinearShiftProbability(GapInfo the_gap,
                                                    double x,
                                                    double y,
                                                    double z) const
  {

    double coordinate;
    switch (the_gap.GetType()) {
    case GapInfo::GapType::X: coordinate = x; break;
    case GapInfo::GapType::Y: coordinate = y; break;
    case GapInfo::GapType::Z: coordinate = z; break;
    default: coordinate = std::numeric_limits<double>::lowest();
    }
    mf::LogInfo("GapChargeTransport")
      << GapInfo::GapToString(the_gap.GetType()) << ", coordinate = " << coordinate;
    double gapWidth = the_gap.GetSize();
    double distance = DistanceToNearestActiveVolume(the_gap, coordinate);
    mf::LogInfo("GapChargeTransport")
      << "GapWidth = " << gapWidth << ", distance = " << distance << ", p2 (prob) = " << fp2_prob
      << ", p1 (dist) =  " << fp1_dist;
    double probability = fp2_prob - distance / fp1_dist * fp2_prob;
    probability = std::clamp(probability, 0.0, 1.0);

    mf::LogInfo("GapChargeTransport") << "Computed probability at (x = " << x << ", y = " << y
                                      << ", z = " << z << ") is " << probability;
    if (the_gap.GetType() == GapInfo::GapType::NotAGap ||
        the_gap.GetType() == GapInfo::GapType::NotInDetectorVolume)
      probability = -1;
    return probability;
  }

  // the main function - here the new postion and the new charge are defined, this is called inside the IonAndScint module
  std::pair<geo::Point_t, int> GapChargeTransport::GetOffset(double x,
                                                             double y,
                                                             double z,
                                                             int n_el) const
  {
    mf::LogInfo("GapChargeTransport")
      << "Calculation Initialized, x =" << x << ", y = " << y << ", z = ";
    mf::LogInfo("GapChargeTransport") << "nelectrons = " << n_el;
    GapInfo the_gap = Gap(x, y, z);
    double probability = ComputeShiftProbability(the_gap, x, y, z);
    int n = round(probability * n_el);
    double shifted_x, shifted_y, shifted_z;
    shifted_x = x;
    shifted_y = y;
    shifted_z = z;
    while (the_gap.GetType() == GapInfo::GapType::X || the_gap.GetType() == GapInfo::GapType::Y ||
           the_gap.GetType() == GapInfo::GapType::Z) {
      double dir = -1.;
      if (the_gap.GetDirection() == GapInfo::MovingDirection::RIGHT) dir = 1;
      mf::LogInfo("GapChargeTransport") << dir;
      mf::LogInfo("GapChargeTransport") << GapInfo::GapToString(the_gap.GetType());
      switch (the_gap.GetType()) {
      case GapInfo::GapType::X:
        shifted_x = shifted_x + dir * (DistanceToNearestActiveVolume(the_gap, x) + fp3_shift);
      case GapInfo::GapType::Y:
        shifted_y = shifted_y + dir * (DistanceToNearestActiveVolume(the_gap, y) + fp3_shift);
      case GapInfo::GapType::Z:
        shifted_z = shifted_z + dir * (DistanceToNearestActiveVolume(the_gap, z) + fp3_shift);
      default: break;
      }
      mf::LogInfo("GapChargeTransport") << "MOVING THE CHARGE TO THE ACTIVE VOLUME " << shifted_x
                                        << " " << shifted_y << " " << shifted_z;
      the_gap = Gap(shifted_x, shifted_y, shifted_z);
    }
    return std::make_pair(geo::Point_t(shifted_x, shifted_y, shifted_z), n);
  }

  //return volume name and throw exception in case there is a typo in the name (in .fcl)
  std::string GapChargeTransport::Volume() const
  {
    if (fp4_volume == "LArG4DetectorServicevolEnclosureTPC") return fp4_volume;
    throw cet::exception("GapChargeTransport")
      << "Incorrect name of the volume containing the gaps";
  }
}

DEFINE_ART_CLASS_TOOL(gap::GapChargeTransport)
