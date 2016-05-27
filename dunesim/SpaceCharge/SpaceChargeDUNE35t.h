////////////////////////////////////////////////////////////////////////
// \file SpaceChargeDUNE35t.h
//
// \brief header of class for storing/accessing space charge distortions for DUNE 35ton
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef SPACECHARGE_SPACECHARGEDUNE35T_H
#define SPACECHARGE_SPACECHARGEDUNE35T_H

// LArSoft libraries
#include "larevt/SpaceCharge/SpaceCharge.h"

// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"

// C/C++ standard libraries
#include <string>
#include <vector>


namespace spacecharge {

  class SpaceChargeDUNE35t : public SpaceCharge {
 
    public:

      explicit SpaceChargeDUNE35t(fhicl::ParameterSet const& pset);
      SpaceChargeDUNE35t(SpaceChargeDUNE35t const&) = delete;
      virtual ~SpaceChargeDUNE35t() = default;
      
      bool Configure(fhicl::ParameterSet const& pset);
      bool Update(uint64_t ts=0);
      
      bool EnableSimulationSCE() const override;
      bool EnableCorrectionsSCE() const override;
      std::vector<double> GetPosOffsets(double xVal, double yVal, double zVal) const override;
 
    private:
    protected:

      std::vector<double> GetPosOffsetsParametric(double xVal, double yVal, double zVal) const;
      double GetOnePosOffsetParametric(double xVal, double yVal, double zVal, std::string axis) const;
      double TransformX(double xVal) const;
      double TransformY(double yVal) const;
      double TransformZ(double zVal) const;
      bool IsInsideBoundaries(double xVal, double yVal, double zVal) const;

      bool fEnableSimulationSCE;
      bool fEnableCorrectionsSCE;
      
      std::string fRepresentationType;
      std::string fInputFilename;
      
      TGraph **g1_x = new TGraph*[7];
      TGraph **g2_x = new TGraph*[7];
      TGraph **g3_x = new TGraph*[7];
      TGraph **g4_x = new TGraph*[7];
      TGraph **g5_x = new TGraph*[7];
      
      TGraph **g1_y = new TGraph*[7];
      TGraph **g2_y = new TGraph*[7];
      TGraph **g3_y = new TGraph*[7];
      TGraph **g4_y = new TGraph*[7];
      TGraph **g5_y = new TGraph*[7];
      TGraph **g6_y = new TGraph*[7];
      
      TGraph **g1_z = new TGraph*[7];
      TGraph **g2_z = new TGraph*[7];
      TGraph **g3_z = new TGraph*[7];
      TGraph **g4_z = new TGraph*[7];
      
      TF1 *f1_x = new TF1("f1_x","pol6");
      TF1 *f2_x = new TF1("f2_x","pol6");
      TF1 *f3_x = new TF1("f3_x","pol6");
      TF1 *f4_x = new TF1("f4_x","pol6");
      TF1 *f5_x = new TF1("f5_x","pol6");
      TF1 *fFinal_x = new TF1("fFinal_x","pol4");
      
      TF1 *f1_y = new TF1("f1_y","pol5");
      TF1 *f2_y = new TF1("f2_y","pol5");
      TF1 *f3_y = new TF1("f3_y","pol5");
      TF1 *f4_y = new TF1("f4_y","pol5");
      TF1 *f5_y = new TF1("f5_y","pol5");
      TF1 *f6_y = new TF1("f6_y","pol5");
      TF1 *fFinal_y = new TF1("fFinal_y","pol5");
      
      TF1 *f1_z = new TF1("f1_z","pol4");
      TF1 *f2_z = new TF1("f2_z","pol4");
      TF1 *f3_z = new TF1("f3_z","pol4");
      TF1 *f4_z = new TF1("f4_z","pol4");
      TF1 *fFinal_z = new TF1("fFinal_z","pol3");
    
  }; // class SpaceChargeDUNE35t
} //namespace spacecharge
#endif // SPACECHARGE_SPACECHARGEDUNE35T_H
