////////////////////////////////////////////////////////////////////////
// \file SpaceChargeProtoDUNEdp.h
//
// \brief header of class for storing/accessing space charge distortions for ProtoDUNE
//
// based on \author mrmooney@bnl.gov
//  \author jdawson@in2p3.fr
////////////////////////////////////////////////////////////////////////
#ifndef SPACECHARGE_SPACECHARGEPROTODUNEDP_H
#define SPACECHARGE_SPACECHARGEPROTODUNEDP_H
// LArSoft libraries
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"
// ROOT includes
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TH3.h"
#include "TTree.h"
#include "TLeaf.h"
// C/C++ standard libraries
#include <string>
#include <vector>
namespace spacecharge {
  class SpaceChargeProtoDUNEdp : public SpaceCharge {
 
    public:
      explicit SpaceChargeProtoDUNEdp(fhicl::ParameterSet const& pset);
      SpaceChargeProtoDUNEdp(SpaceChargeProtoDUNEdp const&) = delete;
      virtual ~SpaceChargeProtoDUNEdp() = default;
      
      bool Configure(fhicl::ParameterSet const& pset, detinfo::DetectorPropertiesData const&);
      bool Update(uint64_t ts=0);
      
      bool EnableSimSpatialSCE() const override;
      bool EnableSimEfieldSCE() const override;
      bool EnableCalSpatialSCE() const override;
      bool EnableCalEfieldSCE() const override;
      
      bool EnableCorrSCE() const override {return (EnableCalSpatialSCE()||EnableCalEfieldSCE()) ;}
      
      geo::Vector_t GetPosOffsets(geo::Point_t const& point) const override;
      geo::Vector_t GetEfieldOffsets(geo::Point_t const& point) const override;
      geo::Vector_t GetCalPosOffsets(geo::Point_t const& point, int const& TPCid) const override;
      geo::Vector_t GetCalEfieldOffsets(geo::Point_t const& point, int const& TPCid) const override;
 
    private:
    protected:
     
      std::vector<double> GetOffsetsVoxel(geo::Point_t const& point, TH3F* hX, TH3F* hY, TH3F* hZ) const;
     
      std::vector<TH3F*> SCEhistograms = std::vector<TH3F*>(7); //Histograms are Dx, Dy, Dz, dEx/E0, dEy/E0, dEz/E0 (positive; repeat for negative)
      std::vector<TH3F*> CalSCEhistograms = std::vector<TH3F*>(6); 
      
      short int driftcoordinate;
      Double_t  Anodebin;

      //    geo::Point_t PretendAtBoundary(geo::Point_t const& point) const; //?
      
      bool fEnableSimSpatialSCE;
      bool fEnableSimEfieldSCE;
      bool fEnableCalSpatialSCE;
      bool fEnableCalEfieldSCE;
      bool fEnableCorrSCE;
      
      double fEfield;
      
      std::string fRepresentationType;
      std::string fInputFilename;
      std::string fCalInputFilename; 
      

    
  }; // class SpaceChargeProtoDUNEdp
} //namespace spacecharge
#endif // SPACECHARGE_SPACECHARGEPROTODUNEDP_H
