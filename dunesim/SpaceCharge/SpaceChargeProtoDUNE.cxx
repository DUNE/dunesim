////////////////////////////////////////////////////////////////////////
// \file SpaceChargeProtoDUNE.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for ProtoDUNE
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////
// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"
// LArSoft includes
#include "dune/SpaceCharge/SpaceChargeProtoDUNE.h"
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
spacecharge::SpaceChargeProtoDUNE::SpaceChargeProtoDUNE(
  fhicl::ParameterSet const& pset
)
{
  Configure(pset);
}
//------------------------------------------------
bool spacecharge::SpaceChargeProtoDUNE::Configure(fhicl::ParameterSet const& pset)
{  
  fEnableSimSpatialSCE = pset.get<bool>("EnableSimSpatialSCE");
  fEnableSimEfieldSCE = pset.get<bool>("EnableSimEfieldSCE");
  fEnableCorrSCE = pset.get<bool>("EnableCorrSCE");
  
  auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fEfield = detprop->Efield();
  
  if((fEnableSimSpatialSCE == true) || (fEnableSimEfieldSCE == true))
  {
    fRepresentationType = pset.get<std::string>("RepresentationType");
    fInputFilename = pset.get<std::string>("InputFilename");
    
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fInputFilename,fname);
    
    std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
    if(!infile->IsOpen()) throw cet::exception("SpaceChargeProtoDUNE") << "Could not find the space charge effect file '" << fname << "'!\n";
    
    if (fRepresentationType == "Voxelized") {
    
      //Load files and read in trees
      if (fInputFilename.find("NegX")<fInputFilename.length()){
        
        TTree* treeD_negX = (TTree*)infile->Get("SpaCEtree_fwdDisp");
        TTree* treeE_negX = (TTree*)infile->Get("SpaCEtree");
        
        fInputFilename.replace(fInputFilename.find("NegX"),3,"Pos");
        
        std::string fname2;
        sp.find_file(fInputFilename,fname2);
        std::unique_ptr<TFile> infile2(new TFile(fname2.c_str(), "READ"));
        if(!infile2->IsOpen()) throw cet::exception("SpaceChargeProtoDUNE") << "Could not find the space charge effect file '" << fname2 << "'!\n";
        
        TTree* treeD_posX = (TTree*)infile2->Get("SpaCEtree_fwdDisp");
        TTree* treeE_posX = (TTree*)infile2->Get("SpaCEtree");
        
        SCEhistograms_negX = Build_TH3(treeD_negX,treeE_negX,"x_true","y_true","z_true","fwd");
        SCEhistograms_posX = Build_TH3(treeD_posX,treeE_posX,"x_true","y_true","z_true","fwd");
      
      }else if (fInputFilename.find("PosX")<fInputFilename.length()){
      
        TTree* treeD_posX = (TTree*)infile->Get("SpaCEtree_fwdDisp");
        TTree* treeE_posX = (TTree*)infile->Get("SpaCEtree");
        
        fInputFilename.replace(fInputFilename.find("PosX"),3,"Neg");
        
        std::string fname2;
        sp.find_file(fInputFilename,fname2);
        std::unique_ptr<TFile> infile2(new TFile(fname2.c_str(), "READ"));
        if(!infile2->IsOpen()) throw cet::exception("SpaceChargeProtoDUNE") << "Could not find the space charge effect file '" << fname2 << "'!\n";
        
        TTree* treeD_negX = (TTree*)infile2->Get("SpaCEtree_fwdDisp");
        TTree* treeE_negX = (TTree*)infile2->Get("SpaCEtree");
      
        SCEhistograms_negX = Build_TH3(treeD_negX,treeE_negX,"x_true","y_true","z_true","fwd");
        SCEhistograms_posX = Build_TH3(treeD_posX,treeE_posX,"x_true","y_true","z_true","fwd");
      
      }else{
      
        TTree* treeD_negX = (TTree*)infile->Get("SpaCEtree_fwdDisp");
        TTree* treeE_negX = (TTree*)infile->Get("SpaCEtree");
               
        TTree* treeD_posX = (TTree*)infile->Get("SpaCEtree_fwdDisp");
        TTree* treeE_posX = (TTree*)infile->Get("SpaCEtree");
        
        SCEhistograms_negX = Build_TH3(treeD_negX,treeE_negX,"x_true","y_true","z_true","fwd");
        SCEhistograms_posX = Build_TH3(treeD_posX,treeE_posX,"x_true","y_true","z_true","fwd");
      
      }
      
    } else if(fRepresentationType == "Parametric")
    {      
    
      for(int i = 0; i < 5; i++)
      {
        g1_x[i] = (TGraph*)infile->Get(Form("deltaX/g1_%d",i));
        g2_x[i] = (TGraph*)infile->Get(Form("deltaX/g2_%d",i));
        g3_x[i] = (TGraph*)infile->Get(Form("deltaX/g3_%d",i));   
        g4_x[i] = (TGraph*)infile->Get(Form("deltaX/g4_%d",i));
        g5_x[i] = (TGraph*)infile->Get(Form("deltaX/g5_%d",i));
        g1_y[i] = (TGraph*)infile->Get(Form("deltaY/g1_%d",i));
        g2_y[i] = (TGraph*)infile->Get(Form("deltaY/g2_%d",i));
        g3_y[i] = (TGraph*)infile->Get(Form("deltaY/g3_%d",i));   
        g4_y[i] = (TGraph*)infile->Get(Form("deltaY/g4_%d",i));
        g5_y[i] = (TGraph*)infile->Get(Form("deltaY/g5_%d",i));
        g6_y[i] = (TGraph*)infile->Get(Form("deltaY/g6_%d",i));
        g1_z[i] = (TGraph*)infile->Get(Form("deltaZ/g1_%d",i));
        g2_z[i] = (TGraph*)infile->Get(Form("deltaZ/g2_%d",i));
        g3_z[i] = (TGraph*)infile->Get(Form("deltaZ/g3_%d",i));   
        g4_z[i] = (TGraph*)infile->Get(Form("deltaZ/g4_%d",i));
        g1_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g1_%d",i));
        g2_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g2_%d",i));
        g3_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g3_%d",i));
        g4_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g4_%d",i));
        g5_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g5_%d",i));
        g1_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g1_%d",i));
        g2_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g2_%d",i));
        g3_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g3_%d",i));
        g4_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g4_%d",i));
        g5_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g5_%d",i));
        g6_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g6_%d",i));
        g1_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g1_%d",i));
        g2_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g2_%d",i));
        g3_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g3_%d",i));
        g4_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g4_%d",i));
      }
      g1_x[5] = (TGraph*)infile->Get("deltaX/g1_5");
      g2_x[5] = (TGraph*)infile->Get("deltaX/g2_5");
      g3_x[5] = (TGraph*)infile->Get("deltaX/g3_5");   
      g4_x[5] = (TGraph*)infile->Get("deltaX/g4_5");
      g5_x[5] = (TGraph*)infile->Get("deltaX/g5_5");
      g1_y[5] = (TGraph*)infile->Get("deltaY/g1_5");
      g2_y[5] = (TGraph*)infile->Get("deltaY/g2_5");
      g3_y[5] = (TGraph*)infile->Get("deltaY/g3_5");   
      g4_y[5] = (TGraph*)infile->Get("deltaY/g4_5");
      g5_y[5] = (TGraph*)infile->Get("deltaY/g5_5");
      g6_y[5] = (TGraph*)infile->Get("deltaY/g6_5");
      
      g1_x[6] = (TGraph*)infile->Get("deltaX/g1_6");
      g2_x[6] = (TGraph*)infile->Get("deltaX/g2_6");
      g3_x[6] = (TGraph*)infile->Get("deltaX/g3_6");
      g4_x[6] = (TGraph*)infile->Get("deltaX/g4_6");
      g5_x[6] = (TGraph*)infile->Get("deltaX/g5_6");
      g1_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g1_5");
      g2_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g2_5");
      g3_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g3_5");
      g4_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g4_5");
      g5_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g5_5");
      g1_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g1_5");
      g2_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g2_5");
      g3_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g3_5");
      g4_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g4_5");
      g5_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g5_5");
      g6_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g6_5");
      g1_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g1_6");
      g2_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g2_6");
      g3_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g3_6");
      g4_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g4_6");
      g5_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g5_6");
    }
    infile->Close();
  }
  if(fEnableCorrSCE == true)
  {
    // Grab other parameters from pset  
  }
  return true;
}
//------------------------------------------------
bool spacecharge::SpaceChargeProtoDUNE::Update(uint64_t ts) 
{
  if (ts == 0) return false;
  return true;
}
//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to turn simulation of SCE on for
/// spatial distortions
bool spacecharge::SpaceChargeProtoDUNE::EnableSimSpatialSCE() const
{
  return fEnableSimSpatialSCE;
}
//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to turn simulation of SCE on for
/// E field distortions
bool spacecharge::SpaceChargeProtoDUNE::EnableSimEfieldSCE() const
{
  return fEnableSimEfieldSCE;
}
//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceChargeProtoDUNE::EnableCorrSCE() const
{
  return fEnableCorrSCE;
}
//----------------------------------------------------------------------------
/// Primary working method of service that provides position offsets to be
/// used in ionization electron drift
geo::Vector_t spacecharge::SpaceChargeProtoDUNE::GetPosOffsets(geo::Point_t const& tmp_point) const
{
  std::vector<double> thePosOffsets;
  geo::Point_t point = tmp_point;
  if(IsTooFarFromBoundaries(point)) {
    thePosOffsets.resize(3,0.0);
    return { -thePosOffsets[0], -thePosOffsets[1], -thePosOffsets[2] };
  }
  if(!IsInsideBoundaries(point)&&!IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);
  
  if (fRepresentationType == "Voxelized"){
    if (point.X()>0){
      thePosOffsets = GetOffsetsVoxel(point, SCEhistograms_posX.at(0), SCEhistograms_posX.at(1), SCEhistograms_posX.at(2));
      //thePosOffsets[0] = -1.0*thePosOffsets[0];
    } else {
      thePosOffsets = GetOffsetsVoxel(point, SCEhistograms_negX.at(0), SCEhistograms_negX.at(1), SCEhistograms_negX.at(2));
    }
      
  }else if(fRepresentationType == "Parametric")
    thePosOffsets = GetPosOffsetsParametric(point.X(), point.Y(), point.Z());
  else thePosOffsets.resize(3,0.0); 

  return { thePosOffsets[0], thePosOffsets[1], thePosOffsets[2] };
}

//----------------------------------------------------------------------------
/// Provides position offsets using voxelized interpolation
std::vector<double> spacecharge::SpaceChargeProtoDUNE::GetOffsetsVoxel
  (geo::Point_t const& point, TH3F* hX, TH3F* hY, TH3F* hZ) const
{
  double xnew = TransformX(point.X());
  double ynew = TransformY(point.Y());
  double znew = TransformZ(point.Z());
  
  return {
    hX->Interpolate(xnew,ynew,znew),
    hY->Interpolate(xnew,ynew,znew),
    hZ->Interpolate(xnew,ynew,znew)
   };
}

/// Build 3d histograms for voxelized interpolation
std::vector<TH3F*> spacecharge::SpaceChargeProtoDUNE::Build_TH3
  (TTree* tree, TTree* eTree, std::string xvar, std::string yvar, std::string zvar, std::string posLeaf) const
{

  //Define the protoDUNE detector
  double Lx = 3.6, Ly = 6.0, Lz = 7.2;
  double numDivisions_x = 18.0;
  double cell_size = Lx/numDivisions_x;
  double numDivisions_y = TMath::Nint((Ly/Lx)*((Double_t)numDivisions_x));
  double numDivisions_z = TMath::Nint((Lz/Lx)*((Double_t)numDivisions_x));
  
  double E_numDivisions_x = 18.0;
  double E_cell_size = Lx/E_numDivisions_x;
  double E_numDivisions_y = TMath::Nint((Ly/Lx)*((Double_t)E_numDivisions_x));
  double E_numDivisions_z = TMath::Nint((Lz/Lx)*((Double_t)E_numDivisions_x));
  
  //initialized histograms for Dx, Dy, Dz, and electric field
  TH3F* hDx = new TH3F("hDx", "", numDivisions_x+1, -0.5*cell_size, Lx+0.5*cell_size, numDivisions_y+1 ,-0.5*cell_size, Ly+0.5*cell_size, numDivisions_z+1, -0.5*cell_size, Lz+0.5*cell_size);
  TH3F* hDy = new TH3F("hDy", "", numDivisions_x+1, -0.5*cell_size, Lx+0.5*cell_size, numDivisions_y+1, -0.5*cell_size, Ly+0.5*cell_size, numDivisions_z+1, -0.5*cell_size, Lz+0.5*cell_size);
  TH3F* hDz = new TH3F("hDz", "", numDivisions_x+1, -0.5*cell_size, Lx+0.5*cell_size, numDivisions_y+1, -0.5*cell_size, Ly+0.5*cell_size, numDivisions_z+1, -0.5*cell_size, Lz+0.5*cell_size);
  
  TH3F* hEx = new TH3F("hEx", "", E_numDivisions_x+1, -0.5*E_cell_size, Lx+0.5*E_cell_size, E_numDivisions_y+1, -0.5*E_cell_size, Ly+0.5*E_cell_size, E_numDivisions_z+1, -0.5*E_cell_size, Lz+0.5*E_cell_size);
  TH3F* hEy = new TH3F("hEy", "", E_numDivisions_x+1, -0.5*E_cell_size, Lx+0.5*E_cell_size, E_numDivisions_y+1, -0.5*E_cell_size, Ly+0.5*E_cell_size, E_numDivisions_z+1, -0.5*E_cell_size, Lz+0.5*E_cell_size);
  TH3F* hEz = new TH3F("hez", "", E_numDivisions_x+1, -0.5*E_cell_size, Lx+0.5*E_cell_size, E_numDivisions_y+1, -0.5*E_cell_size, Ly+0.5*E_cell_size, E_numDivisions_z+1, -0.5*E_cell_size, Lz+0.5*E_cell_size);
 
  //For each event, read the tree and fill each histogram
  for (int ii = 0; ii<tree->GetEntries(); ii++){

    //Read the trees
    tree->GetEntry(ii);
    Double_t x = tree->GetBranch(xvar.c_str())->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t y = tree->GetBranch(yvar.c_str())->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t z = tree->GetBranch(zvar.c_str())->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t dx = tree->GetBranch("Dx")->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t dy = tree->GetBranch("Dy")->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t dz = tree->GetBranch("Dz")->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
		
    eTree->GetEntry(ii);
    Double_t Ex = eTree->GetBranch("Ex")->GetLeaf("data")->GetValue() / (100000.0*fEfield);
    Double_t Ey = eTree->GetBranch("Ey")->GetLeaf("data")->GetValue() / (100000.0*fEfield);
    Double_t Ez = eTree->GetBranch("Ez")->GetLeaf("data")->GetValue() / (100000.0*fEfield);
   
    //Fill the histograms		
    hDx->Fill(x,y,z,100.0*dx);
    hDy->Fill(x,y,z,100.0*dy);
    hDz->Fill(x,y,z,100.0*dz);
    hEx->Fill(x,y,z,Ex);
    hEy->Fill(x,y,z,Ey);
    hEz->Fill(x,y,z,Ez);
  }
  
  return  {hDx, hDy, hDz, hEx, hEy, hEz};

}

//----------------------------------------------------------------------------
/// Provides position offsets using a parametric representation
std::vector<double> spacecharge::SpaceChargeProtoDUNE::GetPosOffsetsParametric(double xVal, double yVal, double zVal) const
{
  std::vector<double> thePosOffsetsParametric;
  double xValNew = TransformX(xVal);
  double yValNew = TransformY(yVal);
  double zValNew = TransformZ(zVal);
  thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew,yValNew,zValNew,"X"));
  thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew,yValNew,zValNew,"Y"));
  thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew,yValNew,zValNew,"Z"));
  return thePosOffsetsParametric;
}
//----------------------------------------------------------------------------
/// Provides one position offset using a parametric representation, for a given
/// axis
double spacecharge::SpaceChargeProtoDUNE::GetOnePosOffsetParametric(double xValNew, double yValNew, double zValNew, std::string axis) const
{      
  double parA[6][7];
  double parB[6];
  
  xValNew -= 1.8;
  yValNew -= 3.0;
  
  for(int j = 0; j < 6; j++)
  {
    for(int i = 0; i < 7; i++)
      parA[j][i] = 0.0;
  
    parB[j] = 0.0;
  }
  
  if(axis == "X")
  {
    for(int j = 0; j < 7; j++)
    {
      parA[0][j] = g1_x[j]->Eval(zValNew);
      parA[1][j] = g2_x[j]->Eval(zValNew);
      parA[2][j] = g3_x[j]->Eval(zValNew);
      parA[3][j] = g4_x[j]->Eval(zValNew);
      parA[4][j] = g5_x[j]->Eval(zValNew);
    }
  
    f1_x->SetParameters(parA[0]);
    f2_x->SetParameters(parA[1]);
    f3_x->SetParameters(parA[2]);
    f4_x->SetParameters(parA[3]);
    f5_x->SetParameters(parA[4]);
  }
  else if(axis == "Y")
  {
    for(int j = 0; j < 6; j++)
    {
      parA[0][j] = g1_y[j]->Eval(zValNew);
      parA[1][j] = g2_y[j]->Eval(zValNew);
      parA[2][j] = g3_y[j]->Eval(zValNew);
      parA[3][j] = g4_y[j]->Eval(zValNew);
      parA[4][j] = g5_y[j]->Eval(zValNew);
      parA[5][j] = g6_y[j]->Eval(zValNew);
    }
  
    f1_y->SetParameters(parA[0]);
    f2_y->SetParameters(parA[1]);
    f3_y->SetParameters(parA[2]);
    f4_y->SetParameters(parA[3]);
    f5_y->SetParameters(parA[4]);
    f6_y->SetParameters(parA[5]);
  }
  else if(axis == "Z")
  {
    for(int j = 0; j < 5; j++)
    {
      parA[0][j] = g1_z[j]->Eval(zValNew);
      parA[1][j] = g2_z[j]->Eval(zValNew);
      parA[2][j] = g3_z[j]->Eval(zValNew);
      parA[3][j] = g4_z[j]->Eval(zValNew);
    }
  
    f1_z->SetParameters(parA[0]);
    f2_z->SetParameters(parA[1]);
    f3_z->SetParameters(parA[2]);
    f4_z->SetParameters(parA[3]);
  }
  
  double aValNew;
  double bValNew;
  
  if(axis == "Y")
  {
    aValNew = xValNew;
    bValNew = yValNew;
  }
  else
  {
    aValNew = yValNew;
    bValNew = xValNew;
  }
  
  double offsetValNew = 0.0;
  if(axis == "X")
  {
    parB[0] = f1_x->Eval(aValNew);
    parB[1] = f2_x->Eval(aValNew);
    parB[2] = f3_x->Eval(aValNew);
    parB[3] = f4_x->Eval(aValNew);
    parB[4] = f5_x->Eval(aValNew);
  
    fFinal_x->SetParameters(parB);
    offsetValNew = 100.0*fFinal_x->Eval(bValNew);
  }
  else if(axis == "Y")
  {
    parB[0] = f1_y->Eval(aValNew);
    parB[1] = f2_y->Eval(aValNew);
    parB[2] = f3_y->Eval(aValNew);
    parB[3] = f4_y->Eval(aValNew);
    parB[4] = f5_y->Eval(aValNew);
    parB[5] = f6_y->Eval(aValNew);
  
    fFinal_y->SetParameters(parB);
    offsetValNew = 100.0*fFinal_y->Eval(bValNew);
  }
  else if(axis == "Z")
  {
    parB[0] = f1_z->Eval(aValNew);
    parB[1] = f2_z->Eval(aValNew);
    parB[2] = f3_z->Eval(aValNew);
    parB[3] = f4_z->Eval(aValNew);
  
    fFinal_z->SetParameters(parB);
    offsetValNew = 100.0*fFinal_z->Eval(bValNew);
  }
  
  return offsetValNew;
}
//----------------------------------------------------------------------------
/// Primary working method of service that provides E field offsets to be
/// used in charge/light yield calculation (e.g.)
geo::Vector_t spacecharge::SpaceChargeProtoDUNE::GetEfieldOffsets(geo::Point_t const& tmp_point) const
{
  std::vector<double> theEfieldOffsets;
  geo::Point_t point = tmp_point;
  if(IsTooFarFromBoundaries(point)) {
    theEfieldOffsets.resize(3,0.0);
    return { -theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2] };
  }
  if(!IsInsideBoundaries(point)&&!IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);
  
  if (fRepresentationType == "Voxelized"){
    if (point.X() > 0.0) theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms_posX.at(3), SCEhistograms_posX.at(4), SCEhistograms_posX.at(5));
    else {
      theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms_negX.at(3), SCEhistograms_negX.at(4), SCEhistograms_negX.at(5));
      theEfieldOffsets[0] = theEfieldOffsets[0];
    }
  }else if(fRepresentationType == "Parametric")
    theEfieldOffsets = GetEfieldOffsetsParametric(point.X(), point.Y(), point.Z());
  else
    theEfieldOffsets.resize(3,0.0);
  
  return { -theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2] };
}
//----------------------------------------------------------------------------
/// Provides E field offsets using a parametric representation
std::vector<double> spacecharge::SpaceChargeProtoDUNE::GetEfieldOffsetsParametric(double xVal, double yVal, double zVal) const
{
  std::vector<double> theEfieldOffsetsParametric;
  double xValNew = TransformX(xVal);
  double yValNew = TransformY(yVal);
  double zValNew = TransformZ(zVal);
  theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew,yValNew,zValNew,"X"));
  theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew,yValNew,zValNew,"Y"));
  theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew,yValNew,zValNew,"Z"));
  return theEfieldOffsetsParametric;
}
//----------------------------------------------------------------------------
/// Provides one E field offset using a parametric representation, for a given
/// axis
double spacecharge::SpaceChargeProtoDUNE::GetOneEfieldOffsetParametric(double xValNew, double yValNew, double zValNew, std::string axis) const
{      
  xValNew -= 1.8;
  yValNew -= 3.0;

  double parA[6][7];
  double parB[6];
  
  for(int j = 0; j < 6; j++)
  {
    for(int i = 0; i < 7; i++)
      parA[j][i] = 0.0;
  
    parB[j] = 0.0;
  }
  
  if(axis == "X")
  {
    for(int j = 0; j < 7; j++)
    {
      parA[0][j] = g1_Ex[j]->Eval(zValNew);
      parA[1][j] = g2_Ex[j]->Eval(zValNew);
      parA[2][j] = g3_Ex[j]->Eval(zValNew);
      parA[3][j] = g4_Ex[j]->Eval(zValNew);
      parA[4][j] = g5_Ex[j]->Eval(zValNew);
    }
  
    f1_Ex->SetParameters(parA[0]);
    f2_Ex->SetParameters(parA[1]);
    f3_Ex->SetParameters(parA[2]);
    f4_Ex->SetParameters(parA[3]);
    f5_Ex->SetParameters(parA[4]);
  }
  else if(axis == "Y")
  {
    for(int j = 0; j < 6; j++)
    {
      parA[0][j] = g1_Ey[j]->Eval(zValNew);
      parA[1][j] = g2_Ey[j]->Eval(zValNew);
      parA[2][j] = g3_Ey[j]->Eval(zValNew);
      parA[3][j] = g4_Ey[j]->Eval(zValNew);
      parA[4][j] = g5_Ey[j]->Eval(zValNew);
      parA[5][j] = g6_Ey[j]->Eval(zValNew);
    }
  
    f1_Ey->SetParameters(parA[0]);
    f2_Ey->SetParameters(parA[1]);
    f3_Ey->SetParameters(parA[2]);
    f4_Ey->SetParameters(parA[3]);
    f5_Ey->SetParameters(parA[4]);
    f6_Ey->SetParameters(parA[5]);
  }
  else if(axis == "Z")
  {
    for(int j = 0; j < 5; j++)
    {
      parA[0][j] = g1_Ez[j]->Eval(zValNew);
      parA[1][j] = g2_Ez[j]->Eval(zValNew);
      parA[2][j] = g3_Ez[j]->Eval(zValNew);
      parA[3][j] = g4_Ez[j]->Eval(zValNew);
    }
  
    f1_Ez->SetParameters(parA[0]);
    f2_Ez->SetParameters(parA[1]);
    f3_Ez->SetParameters(parA[2]);
    f4_Ez->SetParameters(parA[3]);
  }
  
  double aValNew;
  double bValNew;
  
  if(axis == "Y")
  {
    aValNew = xValNew;
    bValNew = yValNew;
  }
  else
  {
    aValNew = yValNew;
    bValNew = xValNew;
  }
  
  double offsetValNew = 0.0;
  if(axis == "X")
  {
    parB[0] = f1_Ex->Eval(aValNew);
    parB[1] = f2_Ex->Eval(aValNew);
    parB[2] = f3_Ex->Eval(aValNew);
    parB[3] = f4_Ex->Eval(aValNew);
    parB[4] = f5_Ex->Eval(aValNew);
  
    fFinal_Ex->SetParameters(parB);
    offsetValNew = fFinal_Ex->Eval(bValNew);
  }
  else if(axis == "Y")
  {
    parB[0] = f1_Ey->Eval(aValNew);
    parB[1] = f2_Ey->Eval(aValNew);
    parB[2] = f3_Ey->Eval(aValNew);
    parB[3] = f4_Ey->Eval(aValNew);
    parB[4] = f5_Ey->Eval(aValNew);
    parB[5] = f6_Ey->Eval(aValNew);
  
    fFinal_Ey->SetParameters(parB);
    offsetValNew = fFinal_Ey->Eval(bValNew);
  }
  else if(axis == "Z")
  {
    parB[0] = f1_Ez->Eval(aValNew);
    parB[1] = f2_Ez->Eval(aValNew);
    parB[2] = f3_Ez->Eval(aValNew);
    parB[3] = f4_Ez->Eval(aValNew);
  
    fFinal_Ez->SetParameters(parB);
    offsetValNew = fFinal_Ez->Eval(bValNew);
  }
  
  return offsetValNew;
}
//----------------------------------------------------------------------------
/// Transform X to SCE X coordinate:  [0.0,3.6] --> [0.0,3.6]
double spacecharge::SpaceChargeProtoDUNE::TransformX(double xVal) const
{
  double xValNew;
  xValNew = (fabs(xVal)/100.0);
  //xValNew -= 1.8;
  return xValNew;
}
//----------------------------------------------------------------------------
/// Transform Y to SCE Y coordinate:  [0.0,6.08] --> [0.0,6.0]
double spacecharge::SpaceChargeProtoDUNE::TransformY(double yVal) const
{
  double yValNew;
  yValNew = (6.00/6.08)*((yVal+0.2)/100.0);
  //yValNew -= 3.0;
  return yValNew;
}
//----------------------------------------------------------------------------
/// Transform Z to SCE Z coordinate:  [0.0,6.97] --> [0.0,7.2]
double spacecharge::SpaceChargeProtoDUNE::TransformZ(double zVal) const
{
  double zValNew;
  zValNew = (7.20/6.97)*((zVal+0.8)/100.0);
  return zValNew;
}
//----------------------------------------------------------------------------
/// Check to see if point is inside boundaries of map (allow to go slightly out of range)
bool spacecharge::SpaceChargeProtoDUNE::IsInsideBoundaries(geo::Point_t const& point) const
{
  return !(
       (point.X() <= -360.0) || (point.X() >= 360.0)
    || (point.Y() <=   -0.2) || (point.Y() >= 607.8)
    || (point.Z() <=   -0.8) || (point.Z() >= 696.2)
    );
} 
  
bool spacecharge::SpaceChargeProtoDUNE::IsTooFarFromBoundaries(geo::Point_t const& point) const
{
  return (
       (point.X() < -360.851) || (point.X() > 360.851)
    || (point.Y() <    -10.2) || (point.Y() >   617.8)
    || (point.Z() <    -10.8) || (point.Z() >   706.2)
    );
}

geo::Point_t spacecharge::SpaceChargeProtoDUNE::PretendAtBoundary(geo::Point_t const& point) const
{
  double x = point.X(), y = point.Y(), z = point.Z();
  
  if      (point.X() <= -360.0) x = -359.99999; 
  else if (point.X() >=  360.0) x =  359.99999;
  
  if      (point.Y() <=  -0.2) y =  -0.19999;
  else if (point.Y() >= 607.8) y = 607.79999;
  
  if      (point.Z() <=  -0.8) z =  -0.79999;
  else if (point.Z() >= 696.2) z = 696.19999;
  
  return {x, y, z};
}
