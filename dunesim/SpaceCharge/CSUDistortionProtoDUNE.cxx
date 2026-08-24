////////////////////////////////////////////////////////////////////////
// \file CSUDistortionProtoDUNE.cxx
////////////////////////////////////////////////////////////////////////
#include "dunesim/SpaceCharge/CSUDistortionProtoDUNE.h"

#include "fhiclcpp/ParameterSet.h"

namespace spacecharge {

  //--------------------------------------------------------------------
  CSUDistortionProtoDUNE::CSUDistortionProtoDUNE(fhicl::ParameterSet const& pset)
    : fMaps{pset} // backend validates config (tolerates tool_type)
  {}

  //--------------------------------------------------------------------
  geo::Point_t CSUDistortionProtoDUNE::Distort(geo::Point_t const& point) const
  {
    // Simulation (forward) spatial offsets. Reproduce the legacy IonAndScint
    // sign convention EXACTLY: X is SUBTRACTED, Y and Z are ADDED. (GetPosOffsets
    // already carries the internal X negation, so the net X displacement is
    // +Dx.) Diverter handled by a separate chain link.
    geo::Vector_t const off = fMaps.GetPosOffsets(point);
    geo::Point_t const distorted{
      point.X() - off.X(), point.Y() + off.Y(), point.Z() + off.Z()};
    return CallNextDistort(distorted);
  }

  //--------------------------------------------------------------------
  geo::Point_t CSUDistortionProtoDUNE::Correct(geo::Point_t const& point) const
  {
    // Calibration (backward) spatial offsets. Drift side from position sign,
    // replacing the legacy hardcoded ProtoDUNE-SP TPCid selection. Same
    // per-component sign convention as Distort (X subtracted, Y/Z added) --
    // VERIFY against a calibration consumer (e.g. Calorimetry) before relying
    // on the correction direction.
    int const side = (point.X() > 0.) ? 1 : -1;
    geo::Vector_t const off = fMaps.GetCalPosOffsets(point, side);
    geo::Point_t const corrected{
      point.X() - off.X(), point.Y() + off.Y(), point.Z() + off.Z()};
    return CallNextCorrect(corrected);
  }

} // namespace spacecharge
