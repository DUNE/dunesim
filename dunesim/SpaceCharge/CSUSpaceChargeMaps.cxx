////////////////////////////////////////////////////////////////////////
// \file CSUSpaceChargeMaps.cxx
//
// Faithful lift of the Voxelized_TH3 / Splines_TH3 paths of
// spacecharge::SpaceChargeProtoDUNE. See CSUSpaceChargeMaps.h.
////////////////////////////////////////////////////////////////////////
#include "dunesim/SpaceCharge/CSUSpaceChargeMaps.h"

#include "lardataalg/DetectorInfo/IgnorableToolConfigKeys.h"

#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"

#include "TFile.h"
#include "TMath.h"
#include "TString.h"

#include <memory>
#include <string>
#include <vector>

namespace spacecharge {

  //--------------------------------------------------------------------
  CSUSpaceChargeMaps::CSUSpaceChargeMaps(fhicl::ParameterSet const& pset)
  {
    // Tolerate the framework service_provider / tool_type keys so the same
    // backend config works inside both the E-field service and the tool, plus
    // the E-field provider's NominalEField, which rides in the same flat pset
    // when this backend is constructed inside CSUEFieldProtoDUNE.
    auto ignorable = detinfo::IgnorableToolConfigKeys();
    ignorable.insert("NominalEField");
    fhicl::Table<Config> const config{pset, ignorable};

    fEnableSimSpatialSCE = config().EnableSimSpatialSCE();
    fEnableSimEfieldSCE = config().EnableSimEfieldSCE();
    fEnableCalSpatialSCE = config().EnableCalSpatialSCE();
    fEnableCalEfieldSCE = config().EnableCalEfieldSCE();

    fRepresentationType = config().RepresentationType();
    fInputFilename = config().InputFilename();
    fCalInputFilename = config().CalibrationInputFilename();

    fXMax = config().Boundary().XMax();
    auto const yr = config().Boundary().YRange();
    fYLow = yr[0];
    fYHigh = yr[1];
    auto const zr = config().Boundary().ZRange();
    fZLow = zr[0];
    fZHigh = zr[1];
    auto const yf = config().Boundary().YFar();
    fYFarLow = yf[0];
    fYFarHigh = yf[1];
    auto const zf = config().Boundary().ZFar();
    fZFarLow = zf[0];
    fZFarHigh = zf[1];

    bool created_efield_splines = false;

    // ---- Simulation (forward) maps -------------------------------------
    if ((fEnableSimSpatialSCE == true) || (fEnableSimEfieldSCE == true)) {
      std::string fname;
      cet::search_path sp("FW_SEARCH_PATH");
      sp.find_file(fInputFilename, fname);

      std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
      if (!infile->IsOpen())
        throw cet::exception("CSUSpaceChargeMaps")
          << "Could not find the space charge effect file '" << fname << "'!\n";

      if ((fRepresentationType == "Voxelized_TH3") || (fRepresentationType == "Splines_TH3")) {

        TH3F* hDx_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Pos");
        TH3F* hDy_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Pos");
        TH3F* hDz_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Pos");
        TH3F* hEx_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_X_Pos");
        TH3F* hEy_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Pos");
        TH3F* hEz_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Pos");

        TH3F* hDx_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Neg");
        TH3F* hDy_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Neg");
        TH3F* hDz_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Neg");
        TH3F* hEx_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_X_Neg");
        TH3F* hEy_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Neg");
        TH3F* hEz_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Neg");

        TH3F* hDx_sim_pos = (TH3F*)hDx_sim_pos_orig->Clone("hDx_pos");
        TH3F* hDy_sim_pos = (TH3F*)hDy_sim_pos_orig->Clone("hDy_pos");
        TH3F* hDz_sim_pos = (TH3F*)hDz_sim_pos_orig->Clone("hDz_pos");
        TH3F* hEx_sim_pos = (TH3F*)hEx_sim_pos_orig->Clone("hEx_pos");
        TH3F* hEy_sim_pos = (TH3F*)hEy_sim_pos_orig->Clone("hEy_pos");
        TH3F* hEz_sim_pos = (TH3F*)hEz_sim_pos_orig->Clone("hEz_pos");

        TH3F* hDx_sim_neg = (TH3F*)hDx_sim_neg_orig->Clone("hDx_neg");
        TH3F* hDy_sim_neg = (TH3F*)hDy_sim_neg_orig->Clone("hDy_neg");
        TH3F* hDz_sim_neg = (TH3F*)hDz_sim_neg_orig->Clone("hDz_neg");
        TH3F* hEx_sim_neg = (TH3F*)hEx_sim_neg_orig->Clone("hEx_neg");
        TH3F* hEy_sim_neg = (TH3F*)hEy_sim_neg_orig->Clone("hEy_neg");
        TH3F* hEz_sim_neg = (TH3F*)hEz_sim_neg_orig->Clone("hEz_neg");

        hDx_sim_pos->SetDirectory(0);
        hDy_sim_pos->SetDirectory(0);
        hDz_sim_pos->SetDirectory(0);
        hEx_sim_pos->SetDirectory(0);
        hEy_sim_pos->SetDirectory(0);
        hEz_sim_pos->SetDirectory(0);

        hDx_sim_neg->SetDirectory(0);
        hDy_sim_neg->SetDirectory(0);
        hDz_sim_neg->SetDirectory(0);
        hEx_sim_neg->SetDirectory(0);
        hEy_sim_neg->SetDirectory(0);
        hEz_sim_neg->SetDirectory(0);

        if (fRepresentationType == "Splines_TH3") {
          int nBinsX = hDx_sim_pos_orig->GetNbinsX();
          int nBinsY = hDx_sim_pos_orig->GetNbinsY();
          int nBinsZ = hDx_sim_pos_orig->GetNbinsZ();
          for (int y = 1; y <= nBinsY; y++) {
            spline_dx_fwd_neg.push_back(std::vector<TSpline3*>());
            spline_dx_fwd_pos.push_back(std::vector<TSpline3*>());
            for (int z = 1; z <= nBinsZ; z++) {
              spline_dx_fwd_neg.back().push_back(MakeSpline(hDx_sim_neg, 1, y, z, 1, 1));
              spline_dx_fwd_pos.back().push_back(MakeSpline(hDx_sim_pos, 1, y, z, 1, 2));
            }
          }
          for (int x = 1; x <= nBinsX; x++) {
            spline_dy_fwd_neg.push_back(std::vector<TSpline3*>());
            spline_dy_fwd_pos.push_back(std::vector<TSpline3*>());
            for (int z = 1; z <= nBinsZ; z++) {
              spline_dy_fwd_neg.back().push_back(MakeSpline(hDy_sim_neg, 2, x, z, 1, 1));
              spline_dy_fwd_pos.back().push_back(MakeSpline(hDy_sim_pos, 2, x, z, 1, 2));
            }
          }
          for (int x = 1; x <= nBinsX; x++) {
            spline_dz_fwd_neg.push_back(std::vector<TSpline3*>());
            spline_dz_fwd_pos.push_back(std::vector<TSpline3*>());
            for (int y = 1; y <= nBinsY; y++) {
              spline_dz_fwd_neg.back().push_back(MakeSpline(hDz_sim_neg, 3, x, y, 1, 1));
              spline_dz_fwd_pos.back().push_back(MakeSpline(hDz_sim_pos, 3, x, y, 1, 2));
            }
          }

          nBinsX = hEx_sim_pos_orig->GetNbinsX();
          nBinsY = hEx_sim_pos_orig->GetNbinsY();
          nBinsZ = hEx_sim_pos_orig->GetNbinsZ();
          for (int y = 1; y <= nBinsY; y++) {
            spline_dEx_neg.push_back(std::vector<TSpline3*>());
            spline_dEx_pos.push_back(std::vector<TSpline3*>());
            for (int z = 1; z <= nBinsZ; z++) {
              spline_dEx_neg.back().push_back(MakeSpline(hEx_sim_neg, 1, y, z, 3, 1));
              spline_dEx_pos.back().push_back(MakeSpline(hEx_sim_pos, 1, y, z, 3, 2));
            }
          }
          for (int x = 1; x <= nBinsX; x++) {
            spline_dEy_neg.push_back(std::vector<TSpline3*>());
            spline_dEy_pos.push_back(std::vector<TSpline3*>());
            for (int z = 1; z <= nBinsZ; z++) {
              spline_dEy_neg.back().push_back(MakeSpline(hEy_sim_neg, 2, x, z, 3, 1));
              spline_dEy_pos.back().push_back(MakeSpline(hEy_sim_pos, 2, x, z, 3, 2));
            }
          }
          for (int x = 1; x <= nBinsX; x++) {
            spline_dEz_neg.push_back(std::vector<TSpline3*>());
            spline_dEz_pos.push_back(std::vector<TSpline3*>());
            for (int y = 1; y <= nBinsY; y++) {
              spline_dEz_neg.back().push_back(MakeSpline(hEz_sim_neg, 3, x, y, 3, 1));
              spline_dEz_pos.back().push_back(MakeSpline(hEz_sim_pos, 3, x, y, 3, 2));
            }
          }
          created_efield_splines = true;
        }

        SCEhistograms = {hDx_sim_pos, hDy_sim_pos, hDz_sim_pos, hEx_sim_pos, hEy_sim_pos,
                         hEz_sim_pos, hDx_sim_neg, hDy_sim_neg, hDz_sim_neg, hEx_sim_neg,
                         hEy_sim_neg, hEz_sim_neg};
      }
      infile->Close();
    }

    // ---- Calibration (backward) maps -----------------------------------
    if ((fEnableCalSpatialSCE == true) || (fEnableCalEfieldSCE == true)) {
      std::string fname;
      cet::search_path sp("FW_SEARCH_PATH");
      sp.find_file(fCalInputFilename, fname);

      std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
      if (!infile->IsOpen())
        throw cet::exception("CSUSpaceChargeMaps")
          << "Could not find the space charge effect file '" << fname << "'!\n";

      if ((fRepresentationType == "Voxelized_TH3") || (fRepresentationType == "Splines_TH3")) {

        TH3F* hDx_cal_pos_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_X_Pos");
        TH3F* hDy_cal_pos_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_Y_Pos");
        TH3F* hDz_cal_pos_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_Z_Pos");
        TH3F* hEx_cal_pos_orig = (TH3F*)infile->Get("Reco_ElecField_X_Pos");
        TH3F* hEy_cal_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Pos");
        TH3F* hEz_cal_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Pos");

        TH3F* hDx_cal_neg_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_X_Neg");
        TH3F* hDy_cal_neg_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_Y_Neg");
        TH3F* hDz_cal_neg_orig = (TH3F*)infile->Get("RecoBkwd_Displacement_Z_Neg");
        TH3F* hEx_cal_neg_orig = (TH3F*)infile->Get("Reco_ElecField_X_Neg");
        TH3F* hEy_cal_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Neg");
        TH3F* hEz_cal_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Neg");

        TH3F* hDx_cal_pos = (TH3F*)hDx_cal_pos_orig->Clone("hDx_pos");
        TH3F* hDy_cal_pos = (TH3F*)hDy_cal_pos_orig->Clone("hDy_pos");
        TH3F* hDz_cal_pos = (TH3F*)hDz_cal_pos_orig->Clone("hDz_pos");
        TH3F* hEx_cal_pos = (TH3F*)hEx_cal_pos_orig->Clone("hEx_pos");
        TH3F* hEy_cal_pos = (TH3F*)hEy_cal_pos_orig->Clone("hEy_pos");
        TH3F* hEz_cal_pos = (TH3F*)hEz_cal_pos_orig->Clone("hEz_pos");

        TH3F* hDx_cal_neg = (TH3F*)hDx_cal_neg_orig->Clone("hDx_neg");
        TH3F* hDy_cal_neg = (TH3F*)hDy_cal_neg_orig->Clone("hDy_neg");
        TH3F* hDz_cal_neg = (TH3F*)hDz_cal_neg_orig->Clone("hDz_neg");
        TH3F* hEx_cal_neg = (TH3F*)hEx_cal_neg_orig->Clone("hEx_neg");
        TH3F* hEy_cal_neg = (TH3F*)hEy_cal_neg_orig->Clone("hEy_neg");
        TH3F* hEz_cal_neg = (TH3F*)hEz_cal_neg_orig->Clone("hEz_neg");

        hDx_cal_pos->SetDirectory(0);
        hDy_cal_pos->SetDirectory(0);
        hDz_cal_pos->SetDirectory(0);
        hEx_cal_pos->SetDirectory(0);
        hEy_cal_pos->SetDirectory(0);
        hEz_cal_pos->SetDirectory(0);

        hDx_cal_neg->SetDirectory(0);
        hDy_cal_neg->SetDirectory(0);
        hDz_cal_neg->SetDirectory(0);
        hEx_cal_neg->SetDirectory(0);
        hEy_cal_neg->SetDirectory(0);
        hEz_cal_neg->SetDirectory(0);

        if (fRepresentationType == "Splines_TH3") {
          int nBinsX = hDx_cal_pos_orig->GetNbinsX();
          int nBinsY = hDx_cal_pos_orig->GetNbinsY();
          int nBinsZ = hDx_cal_pos_orig->GetNbinsZ();

          for (int y = 1; y <= nBinsY; y++) {
            spline_dx_bkwd_neg.push_back(std::vector<TSpline3*>());
            spline_dx_bkwd_pos.push_back(std::vector<TSpline3*>());
            for (int z = 1; z <= nBinsZ; z++) {
              spline_dx_bkwd_neg.back().push_back(MakeSpline(hDx_cal_neg, 1, y, z, 2, 1));
              spline_dx_bkwd_pos.back().push_back(MakeSpline(hDx_cal_pos, 1, y, z, 2, 2));
            }
          }
          for (int x = 1; x <= nBinsX; x++) {
            spline_dy_bkwd_neg.push_back(std::vector<TSpline3*>());
            spline_dy_bkwd_pos.push_back(std::vector<TSpline3*>());
            for (int z = 1; z <= nBinsZ; z++) {
              spline_dy_bkwd_neg.back().push_back(MakeSpline(hDy_cal_neg, 2, x, z, 2, 1));
              spline_dy_bkwd_pos.back().push_back(MakeSpline(hDy_cal_pos, 2, x, z, 2, 2));
            }
          }
          for (int x = 1; x <= nBinsX; x++) {
            spline_dz_bkwd_neg.push_back(std::vector<TSpline3*>());
            spline_dz_bkwd_pos.push_back(std::vector<TSpline3*>());
            for (int y = 1; y <= nBinsY; y++) {
              spline_dz_bkwd_neg.back().push_back(MakeSpline(hDz_cal_neg, 3, x, y, 2, 1));
              spline_dz_bkwd_pos.back().push_back(MakeSpline(hDz_cal_pos, 3, x, y, 2, 2));
            }
          }
          if (created_efield_splines == false) {
            nBinsX = hEx_cal_neg->GetNbinsX();
            nBinsY = hEx_cal_neg->GetNbinsY();
            nBinsZ = hEx_cal_neg->GetNbinsZ();
            for (int y = 1; y <= nBinsY; y++) {
              spline_dEx_neg.push_back(std::vector<TSpline3*>());
              spline_dEx_pos.push_back(std::vector<TSpline3*>());
              for (int z = 1; z <= nBinsZ; z++) {
                spline_dEx_neg.back().push_back(MakeSpline(hEx_cal_neg, 1, y, z, 3, 1));
                spline_dEx_pos.back().push_back(MakeSpline(hEx_cal_pos, 1, y, z, 3, 2));
              }
            }
            for (int x = 1; x <= nBinsX; x++) {
              spline_dEy_neg.push_back(std::vector<TSpline3*>());
              spline_dEy_pos.push_back(std::vector<TSpline3*>());
              for (int z = 1; z <= nBinsZ; z++) {
                spline_dEy_neg.back().push_back(MakeSpline(hEy_cal_neg, 2, x, z, 3, 1));
                spline_dEy_pos.back().push_back(MakeSpline(hEy_cal_pos, 2, x, z, 3, 2));
              }
            }
            for (int x = 1; x <= nBinsX; x++) {
              spline_dEz_neg.push_back(std::vector<TSpline3*>());
              spline_dEz_pos.push_back(std::vector<TSpline3*>());
              for (int y = 1; y <= nBinsY; y++) {
                spline_dEz_neg.back().push_back(MakeSpline(hEz_cal_neg, 3, x, y, 3, 1));
                spline_dEz_pos.back().push_back(MakeSpline(hEz_cal_pos, 3, x, y, 3, 2));
              }
            }
            created_efield_splines = true;
          }
        }

        CalSCEhistograms = {hDx_cal_pos, hDy_cal_pos, hDz_cal_pos, hEx_cal_pos, hEy_cal_pos,
                            hEz_cal_pos, hDx_cal_neg, hDy_cal_neg, hDz_cal_neg, hEx_cal_neg,
                            hEy_cal_neg, hEz_cal_neg};
      }
      infile->Close();
    }
  }

  //--------------------------------------------------------------------
  geo::Vector_t CSUSpaceChargeMaps::GetPosOffsets(geo::Point_t const& tmp_point) const
  {
    std::vector<double> thePosOffsets;
    geo::Point_t point = tmp_point;
    if (IsTooFarFromBoundaries(point)) { return {0.0, 0.0, 0.0}; }
    if (!IsInsideBoundaries(point) && !IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);

    if (point.X() > 0.) {
      thePosOffsets = GetOffsetsVoxel(point, SCEhistograms.at(0), SCEhistograms.at(1), SCEhistograms.at(2), 1, 2);
      thePosOffsets[0] = -1.0 * thePosOffsets[0];
    }
    else {
      thePosOffsets = GetOffsetsVoxel(point, SCEhistograms.at(6), SCEhistograms.at(7), SCEhistograms.at(8), 1, 1);
      thePosOffsets[0] = -1.0 * thePosOffsets[0];
    }

    // NOTE: electron-diverter offsets are intentionally NOT added here; they
    // are a separate chain link (spacecharge::ElectronDiverterDistortion).
    return {thePosOffsets[0], thePosOffsets[1], thePosOffsets[2]};
  }

  //--------------------------------------------------------------------
  geo::Vector_t CSUSpaceChargeMaps::GetCalPosOffsets(geo::Point_t const& tmp_point, int driftSide) const
  {
    std::vector<double> thePosOffsets;
    geo::Point_t point = tmp_point;

    if (IsTooFarFromBoundaries(point)) { return {0.0, 0.0, 0.0}; }
    if (!IsInsideBoundaries(point) && !IsTooFarFromBoundaries(point)) { point = PretendAtBoundary(point); }

    if ((driftSide == 1) && point.X() > -20.) {
      if (point.X() < 0.) point = {0.00001, point.Y(), point.Z()};
      thePosOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(0), CalSCEhistograms.at(1), CalSCEhistograms.at(2), 2, 2);
      thePosOffsets[0] = -1.0 * thePosOffsets[0];
    }
    else if ((driftSide == -1) && point.X() < 20.) {
      if (point.X() > 0.) point = {-0.00001, point.Y(), point.Z()};
      thePosOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(6), CalSCEhistograms.at(7), CalSCEhistograms.at(8), 2, 1);
      thePosOffsets[0] = -1.0 * thePosOffsets[0];
    }
    else
      thePosOffsets = {0., 0., 0.};

    return {thePosOffsets[0], thePosOffsets[1], thePosOffsets[2]};
  }

  //--------------------------------------------------------------------
  geo::Vector_t CSUSpaceChargeMaps::GetEfieldOffsets(geo::Point_t const& tmp_point) const
  {
    std::vector<double> theEfieldOffsets;
    geo::Point_t point = tmp_point;
    if (IsTooFarFromBoundaries(point)) { return {0.0, 0.0, 0.0}; }
    if (!IsInsideBoundaries(point) && !IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);

    if (point.X() > 0.)
      theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms.at(3), SCEhistograms.at(4), SCEhistograms.at(5), 3, 2);
    else
      theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms.at(9), SCEhistograms.at(10), SCEhistograms.at(11), 3, 1);
    theEfieldOffsets[0] = -1.0 * theEfieldOffsets[0];
    theEfieldOffsets[1] = -1.0 * theEfieldOffsets[1];
    theEfieldOffsets[2] = -1.0 * theEfieldOffsets[2];

    return {-theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2]};
  }

  //--------------------------------------------------------------------
  geo::Vector_t CSUSpaceChargeMaps::GetCalEfieldOffsets(geo::Point_t const& tmp_point, int driftSide) const
  {
    std::vector<double> theEfieldOffsets;
    geo::Point_t point = tmp_point;
    if (IsTooFarFromBoundaries(point)) { return {0.0, 0.0, 0.0}; }
    if (!IsInsideBoundaries(point) && !IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);

    if ((driftSide == 1) && point.X() > -20.) {
      if (point.X() < 0.) point = {0.00001, point.Y(), point.Z()};
      theEfieldOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(3), CalSCEhistograms.at(4), CalSCEhistograms.at(5), 3, 2);
    }
    else if ((driftSide == -1) && point.X() < 20.) {
      if (point.X() > 0.) point = {-0.00001, point.Y(), point.Z()};
      theEfieldOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(9), CalSCEhistograms.at(10), CalSCEhistograms.at(11), 3, 1);
    }
    else
      theEfieldOffsets = {0., 0., 0.};
    theEfieldOffsets[0] = -1.0 * theEfieldOffsets[0];
    theEfieldOffsets[1] = -1.0 * theEfieldOffsets[1];
    theEfieldOffsets[2] = -1.0 * theEfieldOffsets[2];

    return {-theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2]};
  }

  //--------------------------------------------------------------------
  std::vector<double> CSUSpaceChargeMaps::GetOffsetsVoxel(
    geo::Point_t const& point, TH3F* hX, TH3F* hY, TH3F* hZ, int maptype, int driftvol) const
  {
    if (fRepresentationType == "Voxelized_TH3") {
      return {hX->Interpolate(point.X(), point.Y(), point.Z()),
              hY->Interpolate(point.X(), point.Y(), point.Z()),
              hZ->Interpolate(point.X(), point.Y(), point.Z())};
    }
    else { // Splines_TH3
      return {InterpolateSplines(hX, point.X(), point.Y(), point.Z(), 1, maptype, driftvol),
              InterpolateSplines(hY, point.X(), point.Y(), point.Z(), 2, maptype, driftvol),
              InterpolateSplines(hZ, point.X(), point.Y(), point.Z(), 3, maptype, driftvol)};
    }
  }

  //--------------------------------------------------------------------
  bool CSUSpaceChargeMaps::IsInsideBoundaries(geo::Point_t const& point) const
  {
    return !((TMath::Abs(point.X()) <= 0.0) || (TMath::Abs(point.X()) >= fXMax) ||
             (point.Y() <= fYLow) || (point.Y() >= fYHigh) ||
             (point.Z() <= fZLow) || (point.Z() >= fZHigh));
  }

  //--------------------------------------------------------------------
  bool CSUSpaceChargeMaps::IsTooFarFromBoundaries(geo::Point_t const& point) const
  {
    return ((TMath::Abs(point.X()) < -20.0) || (TMath::Abs(point.X()) >= fXMax) ||
            (point.Y() < fYFarLow) || (point.Y() > fYFarHigh) ||
            (point.Z() < fZFarLow) || (point.Z() > fZFarHigh));
  }

  //--------------------------------------------------------------------
  geo::Point_t CSUSpaceChargeMaps::PretendAtBoundary(geo::Point_t const& point) const
  {
    double x = point.X(), y = point.Y(), z = point.Z();

    if (TMath::Abs(point.X()) == 0.0)
      x = -0.00001;
    else if (TMath::Abs(point.X()) < 0.00001)
      x = TMath::Sign(point.X(), 1) * 0.00001;
    else if (TMath::Abs(point.X()) >= fXMax)
      x = TMath::Sign(point.X(), 1) * (fXMax - 0.00001);

    if (point.Y() <= fYLow)
      y = fYLow + 0.00001;
    else if (point.Y() >= fYHigh)
      y = fYHigh - 0.00001;

    if (point.Z() <= fZLow)
      z = fZLow + 0.00001;
    else if (point.Z() >= fZHigh)
      z = fZHigh - 0.00001;

    return {x, y, z};
  }

  //--------------------------------------------------------------------
  TSpline3* CSUSpaceChargeMaps::MakeSpline(
    TH3F* spline_hist, int dim1, int dim2_bin, int dim3_bin, int maptype, int driftvol) const
  {
    TSpline3* spline = 0;

    std::vector<double> a, b;
    if (dim1 == 1) {
      for (int x = 1; x <= spline_hist->GetNbinsX(); ++x) {
        a.push_back(spline_hist->GetXaxis()->GetBinCenter(x));
        b.push_back(spline_hist->GetBinContent(x, dim2_bin, dim3_bin));
      }
    }
    else if (dim1 == 2) {
      for (int y = 1; y <= spline_hist->GetNbinsY(); ++y) {
        a.push_back(spline_hist->GetYaxis()->GetBinCenter(y));
        b.push_back(spline_hist->GetBinContent(dim2_bin, y, dim3_bin));
      }
    }
    else if (dim1 == 3) {
      for (int z = 1; z <= spline_hist->GetNbinsZ(); z++) {
        a.push_back(spline_hist->GetZaxis()->GetBinCenter(z));
        b.push_back(spline_hist->GetBinContent(dim2_bin, dim3_bin, z));
      }
    }
    else {
      cet::exception("CSUSpaceChargeMaps::MakeSpline") << "Unkown dimension " << dim1 << "\n";
    }

    spline = new TSpline3(
      Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin, dim3_bin, maptype, driftvol),
      &a[0], &b[0], a.size(), "b2e2", 0, 0);
    spline->SetName(
      Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin, dim3_bin, maptype, driftvol));

    return spline;
  }

  //--------------------------------------------------------------------
  double CSUSpaceChargeMaps::InterpolateSplines(
    TH3F* interp_hist, double xVal, double yVal, double zVal, int dim, int maptype, int driftvol) const
  {
    int bin_x = interp_hist->GetXaxis()->FindBin(xVal);
    int bin_y = interp_hist->GetYaxis()->FindBin(yVal);
    int bin_z = interp_hist->GetZaxis()->FindBin(zVal);

    int bincenter_x = interp_hist->GetXaxis()->GetBinCenter(bin_x);
    int bincenter_y = interp_hist->GetYaxis()->GetBinCenter(bin_y);
    int bincenter_z = interp_hist->GetZaxis()->GetBinCenter(bin_z);

    int max_x = interp_hist->GetNbinsX();
    int max_y = interp_hist->GetNbinsY();
    int max_z = interp_hist->GetNbinsZ();

    int low_x;
    int high_x;
    if (bin_x <= 1) {
      low_x = 1;
      high_x = 2;
    }
    else if (bin_x >= max_x) {
      low_x = max_x - 1;
      high_x = max_x;
    }
    else if (xVal > bincenter_x) {
      low_x = bin_x;
      high_x = bin_x + 1;
    }
    else {
      low_x = bin_x - 1;
      high_x = bin_x;
    }

    int low_y;
    int high_y;
    if (bin_y <= 1) {
      low_y = 1;
      high_y = 2;
    }
    else if (bin_y >= max_y) {
      low_y = max_y - 1;
      high_y = max_y;
    }
    else if (yVal > bincenter_y) {
      low_y = bin_y;
      high_y = bin_y + 1;
    }
    else {
      low_y = bin_y - 1;
      high_y = bin_y;
    }

    int low_z;
    int high_z;
    if (bin_z <= 1) {
      low_z = 1;
      high_z = 2;
    }
    else if (bin_z >= max_z) {
      low_z = max_z - 1;
      high_z = max_z;
    }
    else if (zVal > bincenter_z) {
      low_z = bin_z;
      high_z = bin_z + 1;
    }
    else {
      low_z = bin_z - 1;
      high_z = bin_z;
    }

    double interp_val = 0.0;

    if (dim == 1) {
      double a_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
      double a_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

      double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
      double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

      double f_11 = 0.0;
      double f_12 = 0.0;
      double f_21 = 0.0;
      double f_22 = 0.0;
      if (driftvol == 1) {
        if (maptype == 1) {
          f_11 = spline_dx_fwd_neg[low_y - 1][low_z - 1]->Eval(xVal);
          f_12 = spline_dx_fwd_neg[low_y - 1][high_z - 1]->Eval(xVal);
          f_21 = spline_dx_fwd_neg[high_y - 1][low_z - 1]->Eval(xVal);
          f_22 = spline_dx_fwd_neg[high_y - 1][high_z - 1]->Eval(xVal);
        }
        else if (maptype == 2) {
          f_11 = spline_dx_bkwd_neg[low_y - 1][low_z - 1]->Eval(xVal);
          f_12 = spline_dx_bkwd_neg[low_y - 1][high_z - 1]->Eval(xVal);
          f_21 = spline_dx_bkwd_neg[high_y - 1][low_z - 1]->Eval(xVal);
          f_22 = spline_dx_bkwd_neg[high_y - 1][high_z - 1]->Eval(xVal);
        }
        else if (maptype == 3) {
          f_11 = spline_dEx_neg[low_y - 1][low_z - 1]->Eval(xVal);
          f_12 = spline_dEx_neg[low_y - 1][high_z - 1]->Eval(xVal);
          f_21 = spline_dEx_neg[high_y - 1][low_z - 1]->Eval(xVal);
          f_22 = spline_dEx_neg[high_y - 1][high_z - 1]->Eval(xVal);
        }
      }
      else if (driftvol == 2) {
        if (maptype == 1) {
          f_11 = spline_dx_fwd_pos[low_y - 1][low_z - 1]->Eval(xVal);
          f_12 = spline_dx_fwd_pos[low_y - 1][high_z - 1]->Eval(xVal);
          f_21 = spline_dx_fwd_pos[high_y - 1][low_z - 1]->Eval(xVal);
          f_22 = spline_dx_fwd_pos[high_y - 1][high_z - 1]->Eval(xVal);
        }
        else if (maptype == 2) {
          f_11 = spline_dx_bkwd_pos[low_y - 1][low_z - 1]->Eval(xVal);
          f_12 = spline_dx_bkwd_pos[low_y - 1][high_z - 1]->Eval(xVal);
          f_21 = spline_dx_bkwd_pos[high_y - 1][low_z - 1]->Eval(xVal);
          f_22 = spline_dx_bkwd_pos[high_y - 1][high_z - 1]->Eval(xVal);
        }
        else if (maptype == 3) {
          f_11 = spline_dEx_pos[low_y - 1][low_z - 1]->Eval(xVal);
          f_12 = spline_dEx_pos[low_y - 1][high_z - 1]->Eval(xVal);
          f_21 = spline_dEx_pos[high_y - 1][low_z - 1]->Eval(xVal);
          f_22 = spline_dEx_pos[high_y - 1][high_z - 1]->Eval(xVal);
        }
      }

      interp_val = (f_11 * (a_2 - yVal) * (b_2 - zVal) + f_21 * (yVal - a_1) * (b_2 - zVal) +
                    f_12 * (a_2 - yVal) * (zVal - b_1) + f_22 * (yVal - a_1) * (zVal - b_1)) /
                   ((a_2 - a_1) * (b_2 - b_1));
    }
    else if (dim == 2) {
      double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
      double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

      double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
      double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

      double f_11 = 0.0;
      double f_12 = 0.0;
      double f_21 = 0.0;
      double f_22 = 0.0;
      if (driftvol == 1) {
        if (maptype == 1) {
          f_11 = spline_dy_fwd_neg[low_x - 1][low_z - 1]->Eval(yVal);
          f_12 = spline_dy_fwd_neg[low_x - 1][high_z - 1]->Eval(yVal);
          f_21 = spline_dy_fwd_neg[high_x - 1][low_z - 1]->Eval(yVal);
          f_22 = spline_dy_fwd_neg[high_x - 1][high_z - 1]->Eval(yVal);
        }
        else if (maptype == 2) {
          f_11 = spline_dy_bkwd_neg[low_x - 1][low_z - 1]->Eval(yVal);
          f_12 = spline_dy_bkwd_neg[low_x - 1][high_z - 1]->Eval(yVal);
          f_21 = spline_dy_bkwd_neg[high_x - 1][low_z - 1]->Eval(yVal);
          f_22 = spline_dy_bkwd_neg[high_x - 1][high_z - 1]->Eval(yVal);
        }
        else if (maptype == 3) {
          f_11 = spline_dEy_neg[low_x - 1][low_z - 1]->Eval(yVal);
          f_12 = spline_dEy_neg[low_x - 1][high_z - 1]->Eval(yVal);
          f_21 = spline_dEy_neg[high_x - 1][low_z - 1]->Eval(yVal);
          f_22 = spline_dEy_neg[high_x - 1][high_z - 1]->Eval(yVal);
        }
      }
      else if (driftvol == 2) {
        if (maptype == 1) {
          f_11 = spline_dy_fwd_pos[low_x - 1][low_z - 1]->Eval(yVal);
          f_12 = spline_dy_fwd_pos[low_x - 1][high_z - 1]->Eval(yVal);
          f_21 = spline_dy_fwd_pos[high_x - 1][low_z - 1]->Eval(yVal);
          f_22 = spline_dy_fwd_pos[high_x - 1][high_z - 1]->Eval(yVal);
        }
        else if (maptype == 2) {
          f_11 = spline_dy_bkwd_pos[low_x - 1][low_z - 1]->Eval(yVal);
          f_12 = spline_dy_bkwd_pos[low_x - 1][high_z - 1]->Eval(yVal);
          f_21 = spline_dy_bkwd_pos[high_x - 1][low_z - 1]->Eval(yVal);
          f_22 = spline_dy_bkwd_pos[high_x - 1][high_z - 1]->Eval(yVal);
        }
        else if (maptype == 3) {
          f_11 = spline_dEy_pos[low_x - 1][low_z - 1]->Eval(yVal);
          f_12 = spline_dEy_pos[low_x - 1][high_z - 1]->Eval(yVal);
          f_21 = spline_dEy_pos[high_x - 1][low_z - 1]->Eval(yVal);
          f_22 = spline_dEy_pos[high_x - 1][high_z - 1]->Eval(yVal);
        }
      }

      interp_val = (f_11 * (a_2 - xVal) * (b_2 - zVal) + f_21 * (xVal - a_1) * (b_2 - zVal) +
                    f_12 * (a_2 - xVal) * (zVal - b_1) + f_22 * (xVal - a_1) * (zVal - b_1)) /
                   ((a_2 - a_1) * (b_2 - b_1));
    }
    else if (dim == 3) {
      double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
      double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

      double b_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
      double b_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

      double f_11 = 0.0;
      double f_12 = 0.0;
      double f_21 = 0.0;
      double f_22 = 0.0;
      if (driftvol == 1) {
        if (maptype == 1) {
          f_11 = spline_dz_fwd_neg[low_x - 1][low_y - 1]->Eval(zVal);
          f_12 = spline_dz_fwd_neg[low_x - 1][high_y - 1]->Eval(zVal);
          f_21 = spline_dz_fwd_neg[high_x - 1][low_y - 1]->Eval(zVal);
          f_22 = spline_dz_fwd_neg[high_x - 1][high_y - 1]->Eval(zVal);
        }
        else if (maptype == 2) {
          f_11 = spline_dz_bkwd_neg[low_x - 1][low_y - 1]->Eval(zVal);
          f_12 = spline_dz_bkwd_neg[low_x - 1][high_y - 1]->Eval(zVal);
          f_21 = spline_dz_bkwd_neg[high_x - 1][low_y - 1]->Eval(zVal);
          f_22 = spline_dz_bkwd_neg[high_x - 1][high_y - 1]->Eval(zVal);
        }
        else if (maptype == 3) {
          f_11 = spline_dEz_neg[low_x - 1][low_y - 1]->Eval(zVal);
          f_12 = spline_dEz_neg[low_x - 1][high_y - 1]->Eval(zVal);
          f_21 = spline_dEz_neg[high_x - 1][low_y - 1]->Eval(zVal);
          f_22 = spline_dEz_neg[high_x - 1][high_y - 1]->Eval(zVal);
        }
      }
      else if (driftvol == 2) {
        if (maptype == 1) {
          f_11 = spline_dz_fwd_pos[low_x - 1][low_y - 1]->Eval(zVal);
          f_12 = spline_dz_fwd_pos[low_x - 1][high_y - 1]->Eval(zVal);
          f_21 = spline_dz_fwd_pos[high_x - 1][low_y - 1]->Eval(zVal);
          f_22 = spline_dz_fwd_pos[high_x - 1][high_y - 1]->Eval(zVal);
        }
        else if (maptype == 2) {
          f_11 = spline_dz_bkwd_pos[low_x - 1][low_y - 1]->Eval(zVal);
          f_12 = spline_dz_bkwd_pos[low_x - 1][high_y - 1]->Eval(zVal);
          f_21 = spline_dz_bkwd_pos[high_x - 1][low_y - 1]->Eval(zVal);
          f_22 = spline_dz_bkwd_pos[high_x - 1][high_y - 1]->Eval(zVal);
        }
        else if (maptype == 3) {
          f_11 = spline_dEz_pos[low_x - 1][low_y - 1]->Eval(zVal);
          f_12 = spline_dEz_pos[low_x - 1][high_y - 1]->Eval(zVal);
          f_21 = spline_dEz_pos[high_x - 1][low_y - 1]->Eval(zVal);
          f_22 = spline_dEz_pos[high_x - 1][high_y - 1]->Eval(zVal);
        }
      }

      interp_val = (f_11 * (a_2 - xVal) * (b_2 - yVal) + f_21 * (xVal - a_1) * (b_2 - yVal) +
                    f_12 * (a_2 - xVal) * (yVal - b_1) + f_22 * (xVal - a_1) * (yVal - b_1)) /
                   ((a_2 - a_1) * (b_2 - b_1));
    }

    return interp_val;
  }

} // namespace spacecharge
