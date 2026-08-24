////////////////////////////////////////////////////////////////////////
// \file CSUEFieldProtoDUNE.cxx
////////////////////////////////////////////////////////////////////////
#include "dunesim/SpaceCharge/CSUEFieldProtoDUNE.h"

#include "fhiclcpp/ParameterSet.h"

namespace spacecharge {

  //--------------------------------------------------------------------
  CSUEFieldProtoDUNE::CSUEFieldProtoDUNE(fhicl::ParameterSet const& pset)
    : fMaps{pset}                                 // backend validates the map keys
    , fE0{pset.get<double>("NominalEField")}      // kV/cm; must match DetectorProperties.Efield
  {}

  //--------------------------------------------------------------------
  geo::Vector_t CSUEFieldProtoDUNE::Efield(geo::Point_t const& point) const
  {
    // GetEfieldOffsets returns the fractional dE/E offsets (legacy CSU sign
    // convention). Build the absolute field from the nominal E0.
    //
    // LOAD-BEARING: the per-drift-volume sign and this combination must be
    // pinned by the golden IonAndScint comparison against the pre-refactor
    // DetectorProperties::Efield. If the golden test diverges, this is the
    // first place to revisit.
    // Outside the map's active volume the field is zero: with no field,
    // IonAndScint's ISCalc draws from a zero-mean distribution and the deposit
    // gets 0 electrons (the desired, clean result; the pre-refactor path yielded
    // -1 there). This mirrors the zero-outside-volume behavior of the other
    // E-field providers (Box/DriftVol) and the DetectorProperties fallback.
    if (!fMaps.IsInsideActiveVolume(point)) return {0., 0., 0.};

    geo::Vector_t const f = fMaps.GetEfieldOffsets(point);
    double const s = (point.X() > 0.) ? 1. : -1.; // drift-volume field direction
    return {s * fE0 * (1. + f.X()), fE0 * f.Y(), fE0 * f.Z()};
  }

} // namespace spacecharge
