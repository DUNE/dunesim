////////////////////////////////////////////////////////////////////////
// \file SpaceChargeProtoDUNE-DP.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for ProtoDUNE-DP
//
// based on\author mrmooney@colostate.edu
// \author jdawson@in2p3.fr
////////////////////////////////////////////////////////////////////////
// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"
// LArSoft includes
#include "dune/SpaceCharge/SpaceChargeProtoDUNEdp.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
// Framework includes
#include "cetlib_except/exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
// ROOT includes
#include "TFile.h"
#include "TH3.h"
#include "TTree.h"
#include "TLeaf.h"
//-----------------------------------------------
spacecharge::SpaceChargeProtoDUNEdp::SpaceChargeProtoDUNEdp(
  fhicl::ParameterSet const& pset
)
{
  //Configure(pset);
}
//------------------------------------------------
bool spacecharge::SpaceChargeProtoDUNEdp::Configure(fhicl::ParameterSet const& pset, 
							detinfo::DetectorProperties const* detprop)
{  

  fEnableSimSpatialSCE = pset.get<bool>("EnableSimSpatialSCE");
  fEnableSimEfieldSCE = pset.get<bool>("EnableSimEfieldSCE");
  fEnableCalSpatialSCE = pset.get<bool>("EnableCalSpatialSCE");
  fEnableCalEfieldSCE = pset.get<bool>("EnableCalEfieldSCE");
  
  //auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fEfield = detprop->Efield();
  
  if((fEnableSimSpatialSCE == true) || (fEnableSimEfieldSCE == true))
  {
    //only one representation type implemented for now
    //    fRepresentationType = pset.get<std::string>("RepresentationType");
    fRepresentationType="Voxelized_TH3";

    //x,y, z (in file z is drift... -> convert to x?)
    
    fInputFilename = pset.get<std::string>("InputFilename");
    
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fInputFilename,fname);
    
    std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
    if(!infile->IsOpen()) throw cet::exception("SpaceChargeProtoDUNEdp") << "Could not find the space charge effect file '" << fname << "'!\n";
    
   if (fRepresentationType == "Voxelized_TH3") { 

     //x position is in mm!! not cm...
      //Load in files //hEndPoint ou HEndPointDX (deltax??)
      TH3F* hDx_sim_orig = (TH3F*)infile->Get("hEndPointDZ");
      TH3F* hDy_sim_orig = (TH3F*)infile->Get("hEndPointDY");
      TH3F* hDz_sim_orig = (TH3F*)infile->Get("hEndPointDX");
      TH3F* hEx_sim_orig = (TH3F*)infile->Get("hEz"); //swap z->x
      TH3F* hEy_sim_orig = (TH3F*)infile->Get("hEy");
      TH3F* hEz_sim_orig = (TH3F*)infile->Get("hEx");
      
        
      TH3F* hDx_sim = (TH3F*)hDx_sim_orig->Clone("hDx");
      TH3F* hDy_sim = (TH3F*)hDy_sim_orig->Clone("hDy");
      TH3F* hDz_sim = (TH3F*)hDz_sim_orig->Clone("hDz");
      TH3F* hEx_sim = (TH3F*)hEx_sim_orig->Clone("hEx");
      TH3F* hEy_sim = (TH3F*)hEy_sim_orig->Clone("hEy");
      TH3F* hEz_sim = (TH3F*)hEz_sim_orig->Clone("hEz");
      
      
        
      hDx_sim->SetDirectory(0);
      hDy_sim->SetDirectory(0);
      hDz_sim->SetDirectory(0);
      hEx_sim->SetDirectory(0);
      hEy_sim->SetDirectory(0);
      hEz_sim->SetDirectory(0);
      
      
        
      SCEhistograms = {hDx_sim, hDy_sim, hDz_sim, hEx_sim, hEy_sim, hEz_sim};
                  
   } //other representations not included yet..

  }  
  fEnableCalEfieldSCE = false;
  fEnableCalSpatialSCE = false;
  
  /* if((fEnableCalSpatialSCE == true) || (fEnableCalEfieldSCE == true))
  {
   throw cet::exception("SpaceChargeProtoDUNEdp") << "Calibration files not include yet...'" << fname << "'!\n";
   
   }*/
  
   
  return true;
}
//------------------------------------------------
bool spacecharge::SpaceChargeProtoDUNEdp::Update(uint64_t ts) 
{
  if (ts == 0) return false;
  return true;
}
//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to turn simulation of SCE on for
/// spatial distortions
bool spacecharge::SpaceChargeProtoDUNEdp::EnableSimSpatialSCE() const
{
  return fEnableSimSpatialSCE;
}
//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to turn simulation of SCE on for
/// E field distortions
bool spacecharge::SpaceChargeProtoDUNEdp::EnableSimEfieldSCE() const
{
  return fEnableSimEfieldSCE;
}
//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to apply SCE corrections
//bool spacecharge::SpaceChargeProtoDUNE::EnableCorrSCE() const
//{
//  return fEnableCorrSCE;
//}

/// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceChargeProtoDUNEdp::EnableCalSpatialSCE() const
{
  return fEnableCalSpatialSCE;
}

/// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceChargeProtoDUNEdp::EnableCalEfieldSCE() const
{
  return fEnableCalEfieldSCE;
}
//----------------------------------------------------------------------------
/// Primary working method of service that provides position offsets to be
/// used in ionization electron drift
geo::Vector_t spacecharge::SpaceChargeProtoDUNEdp::GetPosOffsets(geo::Point_t const& tmp_point) const
{

  std::vector<double> thePosOffsets;
  geo::Point_t point = tmp_point;
  if(IsTooFarFromBoundaries(point)) {
    thePosOffsets.resize(3,0.0);
    return { -thePosOffsets[0], -thePosOffsets[1], -thePosOffsets[2] };
  }
  if(!IsInsideBoundaries(point)&&!IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point); //do pretend
  
  if (fRepresentationType=="Voxelized_TH3"){
    
    	thePosOffsets = GetOffsetsVoxel(point, SCEhistograms.at(0), SCEhistograms.at(1), SCEhistograms.at(2));
	//    	thePosOffsets[0] = -1.0*thePosOffsets[0]; //why the flip x??

      
  }  else thePosOffsets.resize(3,0.0); 
 
  return { thePosOffsets[0], thePosOffsets[1], thePosOffsets[2] };
}

//----------------------------------------------------------------------------
/// Primary working method of service that provides position offsets to be
/// used in calibration of space charge
geo::Vector_t spacecharge::SpaceChargeProtoDUNEdp::GetCalPosOffsets(geo::Point_t const& tmp_point, int const& TPCid) const
{
	
  std::vector<double> thePosOffsets;
   geo::Point_t point = tmp_point;

  
  if(IsTooFarFromBoundaries(point)) {
    thePosOffsets.resize(3,0.0);
    return { -thePosOffsets[0], -thePosOffsets[1], -thePosOffsets[2] };
  }
  if(!IsInsideBoundaries(point)&&!IsTooFarFromBoundaries(point)){ 
  	point = PretendAtBoundary(point); 
  }
  
  if (fRepresentationType == "Voxelized_TH3"){

      thePosOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(0), CalSCEhistograms.at(1), CalSCEhistograms.at(2));
      //      thePosOffsets[0] = -1.0*thePosOffsets[0]; why flip?
     } else thePosOffsets.resize(3,0.0);
  
  return { thePosOffsets[0], thePosOffsets[1], thePosOffsets[2] };
}

//----------------------------------------------------------------------------
/// Provides position offsets using voxelized interpolation
std::vector<double> spacecharge::SpaceChargeProtoDUNEdp::GetOffsetsVoxel
  (geo::Point_t const& point, TH3F* hX, TH3F* hY, TH3F* hZ) const
{
  //  if (fRepresentationType == "Voxelized_TH3"){
  
    return {
      hX->Interpolate(point.X(),point.Y(),point.Z()),
      hY->Interpolate(point.X(),point.Y(),point.Z()),
      hZ->Interpolate(point.X(),point.Y(),point.Z())
    };
    
    //  } 
  

}


//----------------------------------------------------------------------------

/// Primary working method of service that provides E field offsets to be
/// used in charge/light yield calculation (e.g.)
geo::Vector_t spacecharge::SpaceChargeProtoDUNEdp::GetEfieldOffsets(geo::Point_t const& tmp_point) const
{

  std::vector<double> theEfieldOffsets;
  geo::Point_t point = tmp_point;
  if(IsTooFarFromBoundaries(point)) {
    theEfieldOffsets.resize(3,0.0);
    return { -theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2] };
  }
  if(!IsInsideBoundaries(point)&&!IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);
  //  SCEhistograms = {hDx_sim, hDy_sim, hDz_sim, hEx_sim, hEy_sim, hEz_sim};
     
  if (fRepresentationType=="Voxelized_TH3"){
     theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms.at(3), SCEhistograms.at(4), SCEhistograms.at(5));
     /*    theEfieldOffsets[0] = -1.0*theEfieldOffsets[0];
    theEfieldOffsets[1] = -1.0*theEfieldOffsets[1];
    theEfieldOffsets[2] = -1.0*theEfieldOffsets[2]; */
  }  else theEfieldOffsets.resize(3,0.0);
    
   return { -theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2] };
}
//----------------------------------------------------------------------------
/// Primary working method of service that provides E field offsets to be
/// used in charge/light yield calculation (e.g.) for calibration
geo::Vector_t spacecharge::SpaceChargeProtoDUNEdp::GetCalEfieldOffsets(geo::Point_t const& tmp_point, int const& TPCid) const
{ 
  std::vector<double> theEfieldOffsets;
  geo::Point_t point = tmp_point;
  if(IsTooFarFromBoundaries(point)) {
    theEfieldOffsets.resize(3,0.0);
    return { -theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2] };
  }
  if(!IsInsideBoundaries(point)&&!IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);
  
  if (fRepresentationType == "Voxelized_TH3"){
 
      theEfieldOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(3), CalSCEhistograms.at(4), CalSCEhistograms.at(5));
   }else
    theEfieldOffsets.resize(3,0.0);
  
  return { -theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2] };
}

//----------------------------------------------------------------------------
//???

/// Transform X to SCE X coordinate:  [0.0,3.6] --> [0.0,3.6]
/*double spacecharge::SpaceChargeProtoDUNEdp::TransformX(double xVal) const
{
  double xValNew;
  xValNew = (fabs(xVal)/100.0);
  //xValNew -= 1.8;
  return xValNew;
}
//----------------------------------------------------------------------------
/// Transform Y to SCE Y coordinate:  [0.0,6.08] --> [0.0,6.0]
double spacecharge::SpaceChargeProtoDUNEdp::TransformY(double yVal) const
{
  double yValNew;
  yValNew = (6.00/6.08)*((yVal+0.2)/100.0);
  //yValNew -= 3.0;
  return yValNew;
}
//----------------------------------------------------------------------------
/// Transform Z to SCE Z coordinate:  [0.0,6.97] --> [0.0,7.2]
double spacecharge::SpaceChargeProtoDUNEdp::TransformZ(double zVal) const
{
  double zValNew;
  zValNew = (7.20/6.97)*((zVal+0.8)/100.0);
  return zValNew;
}
*/
//----------------------------------------------------------------------------
/// Check to see if point is inside boundaries of map (allow to go slightly out of range)
bool spacecharge::SpaceChargeProtoDUNEdp::IsInsideBoundaries(geo::Point_t const& point) const
{
  //point is in cm..

 double driftmin = SCEhistograms.at(3)->GetZaxis()->GetBinLowEdge(0)/10.0;
  double zmin =  SCEhistograms.at(3)->GetXaxis()->GetBinLowEdge(0)/10.0;
  double ymin =  SCEhistograms.at(3)->GetYaxis()->GetBinLowEdge(0)/10.0; //mm->cm
 
  double driftmax = SCEhistograms.at(3)->GetZaxis()->GetBinUpEdge(SCEhistograms.at(3)->GetZaxis()->GetNbins())/10.0;
  double zmax =  SCEhistograms.at(3)->GetXaxis()->GetBinUpEdge(SCEhistograms.at(3)->GetXaxis()->GetNbins())/10.0;
  double ymax =  SCEhistograms.at(3)->GetYaxis()->GetBinUpEdge(SCEhistograms.at(3)->GetYaxis()->GetNbins())/10.0;

  
  //  if(fRepresentationType=="Voxelized_TH3"){
  	return (
		(point.X() <= driftmax || point.X() >= driftmin)
		|| (point.Y() <= ymax || point.Y() >= ymin)
		|| (point.Z() <= zmax || point.Z() >= zmin)
		 );
	//} 
} 
  
bool spacecharge::SpaceChargeProtoDUNEdp::IsTooFarFromBoundaries(geo::Point_t const& point) const
{
  double driftmin = SCEhistograms.at(3)->GetZaxis()->GetBinLowEdge(0)/10.0;
  double zmin =  SCEhistograms.at(3)->GetXaxis()->GetBinLowEdge(0)/10.0;
  double ymin =  SCEhistograms.at(3)->GetYaxis()->GetBinLowEdge(0)/10.0; //mm->cm
 
  double driftmax = SCEhistograms.at(3)->GetZaxis()->GetBinUpEdge(SCEhistograms.at(3)->GetZaxis()->GetNbins())/10.0;
  double zmax =  SCEhistograms.at(3)->GetXaxis()->GetBinUpEdge(SCEhistograms.at(3)->GetXaxis()->GetNbins())/10.0;
  double ymax =  SCEhistograms.at(3)->GetYaxis()->GetBinUpEdge(SCEhistograms.at(3)->GetYaxis()->GetNbins())/10.0;

  
  // if(fRepresentationType=="Voxelized_TH3"){
  	return (
		(point.X() >= driftmax || point.X() < driftmin)
		|| (point.Y() >= ymax || point.Y() < ymin)
		|| (point.Z() >= zmax || point.Z() < zmin)
		 );
	// }
}


geo::Point_t spacecharge::SpaceChargeProtoDUNEdp::PretendAtBoundary(geo::Point_t const& point) const
{
  
  double x = point.X(), y = point.Y(), z = point.Z();
  //  
  /*  if(fRepresentationType=="Voxelized_TH3"){ 
  
    if      (TMath::Abs(point.X()) ==    0.0    ) x =                           -0.00001;
    else if (TMath::Abs(point.X()) <	 0.00001) x =   TMath::Sign(point.X(),1)*0.00001; 
    else if (TMath::Abs(point.X()) >=    360.0  ) x = TMath::Sign(point.X(),1)*359.99999;
  
    if      (point.Y() <=   5.2) y =   5.20001;
    else if (point.Y() >= 604.0) y = 603.99999;
  
    if      (point.Z() <=   -0.5) z =   -0.49999;
    else if (point.Z() >= 695.3) z = 695.29999;
    
  } else { 
  
    if      (TMath::Abs(point.X()) ==    0.0) x =                           -0.00001;
    else if (TMath::Abs(point.X()) <	 0.0) x =   TMath::Sign(point.X(),1)*0.00001; 
    else if (TMath::Abs(point.X()) >=  360.0) x = TMath::Sign(point.X(),1)*359.99999;
  
    if      (point.Y() <=  -0.2) y =  -0.19999;
    else if (point.Y() >= 607.8) y = 607.79999;
  
    if      (point.Z() <=  -0.8) z =  -0.79999;
    else if (point.Z() >= 696.2) z = 696.19999;
    
    }*/
  return {x, y, z};
}
