#include "art/Utilities/ToolMacros.h"

#include "GapChargeTransport.hh"
#include "GapInfo.hh"
#include "larcore/Geometry/Geometry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gap {
  //constructor
  GapChargeTransport::GapChargeTransport(const fhicl::ParameterSet& pset)
  {
    auto const* geom = lar::providerFrom<geo::Geometry>();

    functionType = pset.get<std::string>("functionType", "linear");
    p1 = pset.get<double>("DistanceZero", 0.5);
    p2 = pset.get<double>("ProbabilityBorder", 0.9);
    p3 = pset.get<double>("Shift", 0.01);
    const geo::TPCGeo& one_tpc = geom->TPC(geo::TPCID(0, 0));
    geo::BoxBoundedGeo activeVolume = one_tpc.ActiveBoundingBox();

    activeVolumeSizeX = activeVolume.MaxX() - activeVolume.MinX();
    activeVolumeSizeY = activeVolume.MaxY() - activeVolume.MinY();
    activeVolumeSizeZ = activeVolume.MaxZ() - activeVolume.MinZ();

    totalTPCs = geom->NTPC();

    minX = std::numeric_limits<double>::max();
    maxX = std::numeric_limits<double>::lowest();
    minY = std::numeric_limits<double>::max();
    maxY = std::numeric_limits<double>::lowest();
    minZ = std::numeric_limits<double>::max();
    maxZ = std::numeric_limits<double>::lowest();

    for (int i = 0; i < totalTPCs; ++i) {
      const geo::TPCGeo& tpc = geom->TPC(geo::TPCID(0, i));
      geo::BoxBoundedGeo bbox = tpc.ActiveBoundingBox();
      minX = std::min(minX, bbox.MinX());
      maxX = std::max(maxX, bbox.MaxX());
      minY = std::min(minY, bbox.MinY());
      maxY = std::max(maxY, bbox.MaxY());
      minZ = std::min(minZ, bbox.MinZ());
      maxZ = std::max(maxZ, bbox.MaxZ());
    }
    std::cout << "The bottom left point (min): (" << minX << ", " << minY << ", " << minZ << ")"
              << std::endl;
    std::cout << "The active volume dimensions: x = " << activeVolumeSizeX
              << ", y = " << activeVolumeSizeY << "z = " << activeVolumeSizeZ << std::endl;
    std::cout.flush();
    detector_size_X = maxX - minX;
    detector_size_Y = maxY - minY;
    detector_size_Z = maxZ - minZ;
    nTPCsX = static_cast<int>(detector_size_X / activeVolumeSizeX + 0.5);
    nTPCsY = static_cast<int>(detector_size_Y / activeVolumeSizeY + 0.5);
    nTPCsZ = static_cast<int>(detector_size_Z / activeVolumeSizeZ + 0.5);
    average_gap_Z = (detector_size_Z - (double)nTPCsZ * activeVolumeSizeZ) / (double)(nTPCsZ - 1);
    average_gap_Y = (detector_size_Y - (double)nTPCsZ * activeVolumeSizeY) / (double)(nTPCsY - 1);
    average_gap_X = (detector_size_X - (double)nTPCsX * activeVolumeSizeX) / (double)(nTPCsX - 1);
    if (average_gap_Z > 2 || average_gap_Y > 2 || average_gap_X > 2) {
      std::cout << "WARNING: The gaps are weirdly huge: y = " << average_gap_Y
                << ", z = " << average_gap_Z << ", x = " << average_gap_X << std::endl;
      std::cout.flush();
    }
  }

  //here the parameters of the gap are calculated
  GapInfo GapChargeTransport::Gap(double x, double y, double z) const
  {
    int id_z = static_cast<int>((z - minZ) / (activeVolumeSizeZ));
    int id_y = static_cast<int>((y - minY) / (activeVolumeSizeY));
    int id_x = static_cast<int>((x - minX) / (activeVolumeSizeX));
    std::cout << "The bottom left point (min): (" << minX << ", " << minY << ", " << minZ << ")"
              << std::endl;
    std::cout << "The active volume dimensions: x = " << activeVolumeSizeX
              << ", y = " << activeVolumeSizeY << ", z = " << activeVolumeSizeZ << std::endl;
    std::cout << "id_x = " << id_x << ", id_y = " << id_y << ", id_z = " << id_z << std::endl;
    std::cout.flush();
    int id_total = (id_x) * (nTPCsY) + (id_y) + (id_z) * (nTPCsX * nTPCsY);
    GapInfo::GapType type = GapInfo::GapType::NotAGap;
    GapInfo::MovingDirection direction = GapInfo::MovingDirection::LEFT;
    double left = 0;
    double right = 0;
    if (z < minZ || z > maxZ || y < minY || y > maxY || x < minX || x > maxX ||
        id_total >= totalTPCs) {
      type = GapInfo::GapType::NotInDetectorVolume;
    }
    else {
      const geo::TPCGeo& tpc = geom->TPC(geo::TPCID(0, id_total));
      geo::BoxBoundedGeo bbox = tpc.ActiveBoundingBox();
      std::cout << "id_total = " << id_total << ", MaxX = " << bbox.MaxX()
                << ", MaxY = " << bbox.MaxY() << ", MaxZ = " << bbox.MaxZ() << std::endl;
      std::cout.flush();

      if (z < bbox.MinZ()) { //not in the "right" active volume
        if (id_z > 0) {      //there is active volume on the left
          left = geom->TPC(geo::TPCID(0, id_total - nTPCsX * nTPCsY)).ActiveBoundingBox().MaxZ();
          if (left < z) {
            right = bbox.MinZ();
            type = GapInfo::GapType::Z;
            if (z > (left + right) / 2.) direction = GapInfo::MovingDirection::RIGHT;
          }
        }
      }
      else {
        if (y < bbox.MinY()) {
          if (id_y > 0) {
            left = geom->TPC(geo::TPCID(0, id_total - 1)).ActiveBoundingBox().MaxY();
            if (left < y) {
              type = GapInfo::GapType::Y;
              right = bbox.MinY();
              if (y > (left + right) / 2.) direction = GapInfo::MovingDirection::RIGHT;
            }
          }
        }
        else {
          if (x > bbox.MaxX()) {
            if (id_x > 0) {
              left = geom->TPC(geo::TPCID(0, id_total - nTPCsY)).ActiveBoundingBox().MaxX();
              if (left < x) {
                right = bbox.MinX();
                type = GapInfo::GapType::X;
                if (x > (left + right) / 2.) direction = GapInfo::MovingDirection::RIGHT;
              }
            }
          }
        }
      }
    }
    std::cout << "Type: " << GapInfo::GapToString(type) << ", left = " << left
              << ", right = " << right << ", direction = " << GapInfo::GapToString(direction)
              << std::endl;
    std::cout.flush();
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
    if (functionType == "linear") { return LinearShiftProbability(the_gap, x, y, z); }
    else {
      std::cout << "Unknown function type: " << functionType << ". Defaulting to linear."
                << std::endl;
      std::cout.flush();
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
    std::cout << "left" << the_gap.GetLeftBorder() << " " << x << " " << the_gap.GetRightBorder()
              << std::endl;
    std::cout.flush();
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
    std::cout << GapInfo::GapToString(the_gap.GetType()) << ", coordinate = " << coordinate
              << std::endl;
    std::cout.flush();
    double gapWidth = the_gap.GetSize();
    double distance = DistanceToNearestActiveVolume(the_gap, coordinate);
    std::cout << "GapWidth = " << gapWidth << ", distance = " << distance << ", p2 (prob) = " << p2
              << ", p1 (dist) =  " << p1 << std::endl;
    std::cout.flush();
    double probability = p2 - distance / p1 * p2;
    probability = std::clamp(probability, 0.0, 1.0);

    std::cout << "Computed probability at (x = " << x << ", y = " << y << ", z = " << z << ") is "
              << probability << std::endl;
    std::cout.flush();
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
    std::cout << "Calculation Initialized, x =" << x << ", y = " << y << ", z = " << z << std::endl;
    std::cout << "nelectrons = " << n_el << std::endl;
    std::cout.flush();
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
      std::cout << dir << std::endl;
      std::cout << GapInfo::GapToString(the_gap.GetType()) << std::endl;
      std::cout.flush();
      switch (the_gap.GetType()) {
      case GapInfo::GapType::X:
        shifted_x = shifted_x + dir * (DistanceToNearestActiveVolume(the_gap, x) + p3);
      case GapInfo::GapType::Y:
        shifted_y = shifted_y + dir * (DistanceToNearestActiveVolume(the_gap, y) + p3);
      case GapInfo::GapType::Z:
        shifted_z = shifted_z + dir * (DistanceToNearestActiveVolume(the_gap, z) + p3);
      default: break;
      }
      std::cout << "MOVING THE CHARGE TO THE ACTIVE VOLUME " << shifted_x << " " << shifted_y << " "
                << shifted_z << std::endl;
      std::cout.flush();
      the_gap = Gap(shifted_x, shifted_y, shifted_z);
    }
    return std::make_pair(geo::Point_t(shifted_x, shifted_y, shifted_z), n);
  }
}

DEFINE_ART_CLASS_TOOL(gap::GapChargeTransport)
