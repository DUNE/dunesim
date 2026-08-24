////////////////////////////////////////////////////////////////////////
// \file CSUSpaceChargeMaps.h
//
// \brief Art-free backend holding the CSU ProtoDUNE space-charge map engine
//        (TH3 / Splines representations only), shared by the CSU pluggable
//        E-field provider (CSUEFieldProtoDUNE) and distortion link
//        (CSUDistortionProtoDUNE).
//
// This is a faithful lift of the "Voxelized_TH3" / "Splines_TH3" paths of
// spacecharge::SpaceChargeProtoDUNE (SP/hd/vd, ~identical), with the only
// per-detector differences (the map-boundary box) externalized to fhicl.
// The legacy "Voxelized"/"Parametric" paths, the Transform{X,Y,Z} helpers,
// the Build_TH3 tree-based construction, and the electron-diverter add-on
// are intentionally NOT ported (the diverter is a separate chain link).
//
// Constructed from a fhicl pset only (no DetectorPropertiesData): the file
// load that SpaceChargeProtoDUNE did in Configure() happens in the ctor. The
// nominal E-field is NOT needed here (only the E-field provider uses it).
////////////////////////////////////////////////////////////////////////
#ifndef DUNESIM_SPACECHARGE_CSUSPACECHARGEMAPS_H
#define DUNESIM_SPACECHARGE_CSUSPACECHARGEMAPS_H

#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "TH3.h"
#include "TSpline.h"

#include <string>
#include <vector>

namespace fhicl {
  class ParameterSet;
}

namespace spacecharge {

  class CSUSpaceChargeMaps {
  public:
    // Map-boundary box for one detector (TH3/Splines path). All differences
    // between SpaceChargeProtoDUNE{,hd,vd} live here. The |X| lower bound is
    // always 0 (implicit); the "too far" |X| bound equals XMax.
    struct BoundaryConfig {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<double> XMax{Name("XMax"), Comment("|X| upper bound [cm]")};
      fhicl::Sequence<double, 2> YRange{
        Name("YRange"), Comment("inside-boundary Y {low,high} [cm]")};
      fhicl::Sequence<double, 2> ZRange{
        Name("ZRange"), Comment("inside-boundary Z {low,high} [cm]")};
      fhicl::Sequence<double, 2> YFar{
        Name("YFar"), Comment("too-far Y {low,high} [cm]")};
      fhicl::Sequence<double, 2> ZFar{
        Name("ZFar"), Comment("too-far Z {low,high} [cm]")};
    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool> EnableSimSpatialSCE{Name("EnableSimSpatialSCE"), Comment("")};
      fhicl::Atom<bool> EnableSimEfieldSCE{Name("EnableSimEfieldSCE"), Comment("")};
      fhicl::Atom<bool> EnableCalSpatialSCE{Name("EnableCalSpatialSCE"), Comment("")};
      fhicl::Atom<bool> EnableCalEfieldSCE{Name("EnableCalEfieldSCE"), Comment("")};
      fhicl::Atom<std::string> RepresentationType{
        Name("RepresentationType"), Comment("\"Voxelized_TH3\" or \"Splines_TH3\"")};
      fhicl::Atom<std::string> InputFilename{Name("InputFilename"), Comment("sim map ROOT file")};
      fhicl::Atom<std::string> CalibrationInputFilename{
        Name("CalibrationInputFilename"), Comment("calibration map ROOT file")};
      fhicl::Table<BoundaryConfig> Boundary{Name("Boundary"), Comment("map boundary box")};
    };

    explicit CSUSpaceChargeMaps(fhicl::ParameterSet const& pset);

    bool EnableSimSpatialSCE() const { return fEnableSimSpatialSCE; }
    bool EnableSimEfieldSCE() const { return fEnableSimEfieldSCE; }
    bool EnableCalSpatialSCE() const { return fEnableCalSpatialSCE; }
    bool EnableCalEfieldSCE() const { return fEnableCalEfieldSCE; }

    // driftSide = +1 (positive-drift maps) or -1 (negative-drift maps),
    // replacing the legacy hardcoded ProtoDUNE-SP TPCid selection.
    geo::Vector_t GetPosOffsets(geo::Point_t const& point) const;
    geo::Vector_t GetEfieldOffsets(geo::Point_t const& point) const;
    geo::Vector_t GetCalPosOffsets(geo::Point_t const& point, int driftSide) const;
    geo::Vector_t GetCalEfieldOffsets(geo::Point_t const& point, int driftSide) const;

    /// True if the point is within the map's nominal active region (the
    /// "inside boundaries" box, NOT the looser "too far" box that includes the
    /// PretendAtBoundary extrapolation margin). Outside this, the CSU E-field
    /// provider returns a zero field so the deposit gets a zero-mean RNG draw
    /// (0 electrons) -- approximating the pre-refactor geometry active volume.
    bool IsInsideActiveVolume(geo::Point_t const& point) const
    {
      return IsInsideBoundaries(point);
    }

  private:
    std::vector<double> GetOffsetsVoxel(
      geo::Point_t const& point, TH3F* hX, TH3F* hY, TH3F* hZ, int maptype, int driftvol) const;
    TSpline3* MakeSpline(
      TH3F* spline_hist, int dim1, int dim2_bin, int dim3_bin, int maptype, int driftvol) const;
    double InterpolateSplines(
      TH3F* interp_hist, double xVal, double yVal, double zVal, int dim, int maptype, int driftvol) const;

    bool IsInsideBoundaries(geo::Point_t const& point) const;
    bool IsTooFarFromBoundaries(geo::Point_t const& point) const;
    geo::Point_t PretendAtBoundary(geo::Point_t const& point) const;

    bool fEnableSimSpatialSCE;
    bool fEnableSimEfieldSCE;
    bool fEnableCalSpatialSCE;
    bool fEnableCalEfieldSCE;

    std::string fRepresentationType;
    std::string fInputFilename;
    std::string fCalInputFilename;

    // Boundary box (from Config).
    double fXMax;
    double fYLow, fYHigh, fZLow, fZHigh;
    double fYFarLow, fYFarHigh, fZFarLow, fZFarHigh;

    // Histograms are Dx,Dy,Dz,dEx/E0,dEy/E0,dEz/E0 (positive; repeat for negative).
    std::vector<TH3F*> SCEhistograms = std::vector<TH3F*>(12);
    std::vector<TH3F*> CalSCEhistograms = std::vector<TH3F*>(12);

    std::vector<std::vector<TSpline3*>> spline_dx_fwd_neg;
    std::vector<std::vector<TSpline3*>> spline_dy_fwd_neg;
    std::vector<std::vector<TSpline3*>> spline_dz_fwd_neg;
    std::vector<std::vector<TSpline3*>> spline_dx_bkwd_neg;
    std::vector<std::vector<TSpline3*>> spline_dy_bkwd_neg;
    std::vector<std::vector<TSpline3*>> spline_dz_bkwd_neg;
    std::vector<std::vector<TSpline3*>> spline_dEx_neg;
    std::vector<std::vector<TSpline3*>> spline_dEy_neg;
    std::vector<std::vector<TSpline3*>> spline_dEz_neg;
    std::vector<std::vector<TSpline3*>> spline_dx_fwd_pos;
    std::vector<std::vector<TSpline3*>> spline_dy_fwd_pos;
    std::vector<std::vector<TSpline3*>> spline_dz_fwd_pos;
    std::vector<std::vector<TSpline3*>> spline_dx_bkwd_pos;
    std::vector<std::vector<TSpline3*>> spline_dy_bkwd_pos;
    std::vector<std::vector<TSpline3*>> spline_dz_bkwd_pos;
    std::vector<std::vector<TSpline3*>> spline_dEx_pos;
    std::vector<std::vector<TSpline3*>> spline_dEy_pos;
    std::vector<std::vector<TSpline3*>> spline_dEz_pos;

  }; // class CSUSpaceChargeMaps

} // namespace spacecharge

#endif // DUNESIM_SPACECHARGE_CSUSPACECHARGEMAPS_H
