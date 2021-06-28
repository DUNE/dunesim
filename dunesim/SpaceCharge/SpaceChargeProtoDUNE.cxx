////////////////////////////////////////////////////////////////////////
// \file SpaceChargeProtoDUNE.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for ProtoDUNE
//
// \author mrmooney@colostate.edu
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
#include "TSpline.h"
//-----------------------------------------------
spacecharge::SpaceChargeProtoDUNE::SpaceChargeProtoDUNE(
  fhicl::ParameterSet const& pset
)
{
  //Configure(pset);
}
//------------------------------------------------

//bool spacecharge::SpaceChargeProtoDUNE::Configure(fhicl::ParameterSet const& pset, detinfo::DetectorProperties const* detprop)
bool spacecharge::SpaceChargeProtoDUNE::Configure(fhicl::ParameterSet const& pset, 
                                                  detinfo::DetectorPropertiesData const& detProp)
{  

  fEnableSimSpatialSCE = pset.get<bool>("EnableSimSpatialSCE");
  fEnableSimEfieldSCE = pset.get<bool>("EnableSimEfieldSCE");
  fEnableCalSpatialSCE = pset.get<bool>("EnableCalSpatialSCE");
  fEnableCalEfieldSCE = pset.get<bool>("EnableCalEfieldSCE");
  
  fEnableElectronDiverterDistortions = pset.get<std::vector<bool>>("EnableElectronDiverterDistortions");
  fEDZCenter                         = pset.get<std::vector<double>>("EDZCenter");
  fEDAXPosOffs                       = pset.get<std::vector<double>>("EDAXPosOffs");
  fEDBZPosOffs                       = pset.get<std::vector<double>>("EDBZPosOffs");
  fEDs                               = pset.get<std::vector<double>>("EDs");
  fEDChargeLossZLow                  = pset.get<std::vector<double>>("EDChargeLossZLow");
  fEDChargeLossZHigh                 = pset.get<std::vector<double>>("EDChargeLossZHigh");

  size_t ieds = fEnableElectronDiverterDistortions.size();
  if (fEDZCenter.size() != ieds ||
      fEDAXPosOffs.size() != ieds ||
      fEDBZPosOffs.size() != ieds ||
      fEDs.size() != ieds ||
      fEDChargeLossZLow.size() != ieds ||
      fEDChargeLossZHigh.size() != ieds)
    {
      throw cet::exception("SpaceChargeProtoDUNE") << "Inconsistent configuration sizes: " <<
	ieds << " " << 
	fEDAXPosOffs.size() << " " <<
        fEDBZPosOffs.size() << " " <<
        fEDs.size() << " " <<
        fEDChargeLossZLow.size() << " " << 
	fEDChargeLossZHigh.size();
    }

  //auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  //fEfield = detprop->Efield();

  bool created_efield_splines = false;

  fEfield = detProp.Efield();
  
  if((fEnableSimSpatialSCE == true) || (fEnableSimEfieldSCE == true))
  {
    fRepresentationType = pset.get<std::string>("RepresentationType");
    fInputFilename = pset.get<std::string>("InputFilename");
    
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fInputFilename,fname);
    
    std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
    if(!infile->IsOpen()) throw cet::exception("SpaceChargeProtoDUNE") << "Could not find the space charge effect file '" << fname << "'!\n";
    
    if ((fRepresentationType == "Voxelized_TH3") || (fRepresentationType == "Splines_TH3")) { 
        
      //Load in files
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
        for(int y = 1; y <= nBinsY; y++){
          spline_dx_fwd_neg.push_back(std::vector<TSpline3*>());
          spline_dx_fwd_pos.push_back(std::vector<TSpline3*>());
          for(int z = 1; z <= nBinsZ; z++){
            spline_dx_fwd_neg.back().push_back(MakeSpline(hDx_sim_neg,1,y,z,1,1));
            spline_dx_fwd_pos.back().push_back(MakeSpline(hDx_sim_pos,1,y,z,1,2));
          }
        }
        for(int x = 1; x <= nBinsX; x++){
          spline_dy_fwd_neg.push_back(std::vector<TSpline3*>());
          spline_dy_fwd_pos.push_back(std::vector<TSpline3*>());

          for(int z = 1; z <= nBinsZ; z++){
            spline_dy_fwd_neg.back().push_back(MakeSpline(hDy_sim_neg,2,x,z,1,1));
            spline_dy_fwd_pos.back().push_back(MakeSpline(hDy_sim_pos,2,x,z,1,2));
          }
        }
        for(int x = 1; x <= nBinsX; x++){
          spline_dz_fwd_neg.push_back(std::vector<TSpline3*>());
          spline_dz_fwd_pos.push_back(std::vector<TSpline3*>());
          for(int y = 1; y <= nBinsY; y++){
            spline_dz_fwd_neg.back().push_back(MakeSpline(hDz_sim_neg,3,x,y,1,1));
            spline_dz_fwd_pos.back().push_back(MakeSpline(hDz_sim_pos,3,x,y,1,2));
          }
        }

        nBinsX = hEx_sim_pos_orig->GetNbinsX();
        nBinsY = hEx_sim_pos_orig->GetNbinsY();
        nBinsZ = hEx_sim_pos_orig->GetNbinsZ();
        for(int y = 1; y <= nBinsY; y++){
          spline_dEx_neg.push_back(std::vector<TSpline3*>());
          spline_dEx_pos.push_back(std::vector<TSpline3*>());
          for(int z = 1; z <= nBinsZ; z++){
            spline_dEx_neg.back().push_back(MakeSpline(hEx_sim_neg,1,y,z,3,1));
            spline_dEx_pos.back().push_back(MakeSpline(hEx_sim_pos,1,y,z,3,2));
          }
        }
        for(int x = 1; x <= nBinsX; x++){
          spline_dEy_neg.push_back(std::vector<TSpline3*>());
          spline_dEy_pos.push_back(std::vector<TSpline3*>());

          for(int z = 1; z <= nBinsZ; z++){
            spline_dEy_neg.back().push_back(MakeSpline(hEy_sim_neg,2,x,z,3,1));
            spline_dEy_pos.back().push_back(MakeSpline(hEy_sim_pos,2,x,z,3,2));
          }
        }
        for(int x = 1; x <= nBinsX; x++){
          spline_dEz_neg.push_back(std::vector<TSpline3*>());
          spline_dEz_pos.push_back(std::vector<TSpline3*>());
          for(int y = 1; y <= nBinsY; y++){
            spline_dEz_neg.back().push_back(MakeSpline(hEz_sim_neg,3,x,y,3,1));
            spline_dEz_pos.back().push_back(MakeSpline(hEz_sim_pos,3,x,y,3,2));
          }
        }
        created_efield_splines = true;
      }
                
      SCEhistograms = {hDx_sim_pos, hDy_sim_pos, hDz_sim_pos, hEx_sim_pos, hEy_sim_pos, hEz_sim_pos, hDx_sim_neg, hDy_sim_neg, hDz_sim_neg, hEx_sim_neg, hEy_sim_neg, hEz_sim_neg};
                  
   } else if (fRepresentationType == "Voxelized") {
    
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
        
        std::vector<TH3F*> temp = Build_TH3(treeD_posX, treeE_posX, treeD_negX, treeE_negX, "x_true", "y_true", "z_true", "fwd");
        
        for (size_t ii = 0; ii<temp.size(); ii++){
        	SCEhistograms.at(ii) = (TH3F*)temp.at(ii)->Clone(TString::Format("%s",temp.at(ii)->GetName()));
        	SCEhistograms.at(ii)->SetDirectory(0);
        }   
        
        
        infile2->Close();
      
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
      
        std::vector<TH3F*> temp = Build_TH3(treeD_posX, treeE_posX, treeD_negX, treeE_negX, "x_true", "y_true", "z_true", "fwd");
        
        for (size_t ii = 0; ii<temp.size(); ii++){
        	SCEhistograms.at(ii) = (TH3F*)temp.at(ii)->Clone(TString::Format("%s",temp.at(ii)->GetName()));
        	SCEhistograms.at(ii)->SetDirectory(0);
        }   
        
        infile2->Close();
      
      }else{
      
        TTree* treeD_negX = (TTree*)infile->Get("SpaCEtree_fwdDisp");
        TTree* treeE_negX = (TTree*)infile->Get("SpaCEtree");
               
        TTree* treeD_posX = (TTree*)infile->Get("SpaCEtree_fwdDisp");
        TTree* treeE_posX = (TTree*)infile->Get("SpaCEtree");
        
        std::vector<TH3F*> temp = Build_TH3(treeD_posX, treeE_posX, treeD_negX, treeE_negX, "x_true", "y_true", "z_true", "fwd");
        
        for (size_t ii = 0; ii<temp.size(); ii++){
        	SCEhistograms.at(ii) = (TH3F*)temp.at(ii)->Clone(TString::Format("%s",temp.at(ii)->GetName()));
        	SCEhistograms.at(ii)->SetDirectory(0);
        }           
            
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
  
  if((fEnableCalSpatialSCE == true) || (fEnableCalEfieldSCE == true))
  {
  
    fRepresentationType = pset.get<std::string>("RepresentationType");
    fCalInputFilename = pset.get<std::string>("CalibrationInputFilename");
    
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fCalInputFilename,fname);
    
    std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
    if(!infile->IsOpen()) throw cet::exception("SpaceChargeProtoDUNE") << "Could not find the space charge effect file '" << fname << "'!\n";
  
    if ((fRepresentationType == "Voxelized_TH3") || (fRepresentationType == "Splines_TH3")) { 
   
      //Load in files
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

        //TFile spline_file("splines.root", "RECREATE");
        //gROOT->SetBatch(1);
        for(int y = 1; y <= nBinsY; y++){
          spline_dx_bkwd_neg.push_back(std::vector<TSpline3*>());
          spline_dx_bkwd_pos.push_back(std::vector<TSpline3*>());
          for(int z = 1; z <= nBinsZ; z++){
            spline_dx_bkwd_neg.back().push_back(MakeSpline(hDx_cal_neg,1,y,z,2,1));
            spline_dx_bkwd_pos.back().push_back(MakeSpline(hDx_cal_pos,1,y,z,2,2));
          }
        }
        for(int x = 1; x <= nBinsX; x++){
          spline_dy_bkwd_neg.push_back(std::vector<TSpline3*>());
          spline_dy_bkwd_pos.push_back(std::vector<TSpline3*>());
          for(int z = 1; z <= nBinsZ; z++){
            spline_dy_bkwd_neg.back().push_back(MakeSpline(hDy_cal_neg,2,x,z,2,1));
            spline_dy_bkwd_pos.back().push_back(MakeSpline(hDy_cal_pos,2,x,z,2,2));
          }
        }
        for(int x = 1; x <= nBinsX; x++){
          spline_dz_bkwd_neg.push_back(std::vector<TSpline3*>());
          spline_dz_bkwd_pos.push_back(std::vector<TSpline3*>());
          for(int y = 1; y <= nBinsY; y++){
            spline_dz_bkwd_neg.back().push_back(MakeSpline(hDz_cal_neg,3,x,y,2,1));
            spline_dz_bkwd_pos.back().push_back(MakeSpline(hDz_cal_pos,3,x,y,2,2));
          }
        }
        if(created_efield_splines == false){
          nBinsX = hEx_cal_neg->GetNbinsX();
          nBinsY = hEx_cal_neg->GetNbinsY();
          nBinsZ = hEx_cal_neg->GetNbinsZ();
          for(int y = 1; y <= nBinsY; y++){
            spline_dEx_neg.push_back(std::vector<TSpline3*>());
            spline_dEx_pos.push_back(std::vector<TSpline3*>());
            for(int z = 1; z <= nBinsZ; z++){
              spline_dEx_neg.back().push_back(MakeSpline(hEx_cal_neg,1,y,z,3,1));
              spline_dEx_pos.back().push_back(MakeSpline(hEx_cal_pos,1,y,z,3,2));
            }
          }
          for(int x = 1; x <= nBinsX; x++){
            spline_dEy_neg.push_back(std::vector<TSpline3*>());
            spline_dEy_pos.push_back(std::vector<TSpline3*>());
            for(int z = 1; z <= nBinsZ; z++){
              spline_dEy_neg.back().push_back(MakeSpline(hEy_cal_neg,2,x,z,3,1));
              spline_dEy_pos.back().push_back(MakeSpline(hEy_cal_pos,2,x,z,3,2));
            }
          }
          for(int x = 1; x <= nBinsX; x++){
            spline_dEz_neg.push_back(std::vector<TSpline3*>());
            spline_dEz_pos.push_back(std::vector<TSpline3*>());
            for(int y = 1; y <= nBinsY; y++){
              spline_dEz_neg.back().push_back(MakeSpline(hEz_cal_neg,3,x,y,3,1));
              spline_dEz_pos.back().push_back(MakeSpline(hEz_cal_pos,3,x,y,3,2));
            }
          }
          created_efield_splines = true;
        }
        //spline_file.Close();
      }
                
      CalSCEhistograms = {hDx_cal_pos, hDy_cal_pos, hDz_cal_pos, hEx_cal_pos, hEy_cal_pos, hEz_cal_pos, hDx_cal_neg, hDy_cal_neg, hDz_cal_neg, hEx_cal_neg, hEy_cal_neg, hEz_cal_neg};
      
    } else if (fRepresentationType == "Voxelized") {
    
      //Load files and read in trees
      if (fCalInputFilename.find("NegX")<fCalInputFilename.length()){
        
        TTree* treeD_negX = (TTree*)infile->Get("SpaCEtree_bkwdDisp");
        TTree* treeE_negX = (TTree*)infile->Get("SpaCEtree");
        
        fCalInputFilename.replace(fCalInputFilename.find("NegX"),3,"Pos");
        
        std::string fname2;
        sp.find_file(fCalInputFilename,fname2);
        std::unique_ptr<TFile> infile2(new TFile(fname2.c_str(), "READ"));
        if(!infile2->IsOpen()) throw cet::exception("SpaceChargeProtoDUNE") << "Could not find the space charge effect file '" << fname2 << "'!\n";
        
        TTree* treeD_posX = (TTree*)infile2->Get("SpaCEtree_bkwdDisp");
        TTree* treeE_posX = (TTree*)infile2->Get("SpaCEtree");
        
        std::vector<TH3F*> temp = Build_TH3(treeD_posX, treeE_posX, treeD_negX,treeE_negX, "x_reco", "y_reco", "z_reco", "bkwd");        
        for (size_t ii = 0; ii<temp.size(); ii++){
        	CalSCEhistograms.at(ii) = (TH3F*)temp.at(ii)->Clone(TString::Format("%s",temp.at(ii)->GetName()));
        	CalSCEhistograms.at(ii)->SetDirectory(0);
        }   
        
        infile2->Close();
      
      }else if (fCalInputFilename.find("PosX")<fCalInputFilename.length()){
      
        TTree* treeD_posX = (TTree*)infile->Get("SpaCEtree_bkwdDisp");
        TTree* treeE_posX = (TTree*)infile->Get("SpaCEtree");
        
        fCalInputFilename.replace(fCalInputFilename.find("PosX"),3,"Neg");
        
        std::string fname2;
        sp.find_file(fCalInputFilename,fname2);
        std::unique_ptr<TFile> infile2(new TFile(fname2.c_str(), "READ"));
        if(!infile2->IsOpen()) throw cet::exception("SpaceChargeProtoDUNE") << "Could not find the space charge effect file '" << fname2 << "'!\n";
        
        TTree* treeD_negX = (TTree*)infile2->Get("SpaCEtree_bkwdDisp");
        TTree* treeE_negX = (TTree*)infile2->Get("SpaCEtree");
      
        std::vector<TH3F*> temp = Build_TH3(treeD_posX, treeE_posX, treeD_negX,treeE_negX, "x_reco", "y_reco", "z_reco", "bkwd");        
        for (size_t ii = 0; ii<temp.size(); ii++){
        	CalSCEhistograms.at(ii) = (TH3F*)temp.at(ii)->Clone(TString::Format("%s",temp.at(ii)->GetName()));
        	CalSCEhistograms.at(ii)->SetDirectory(0);
        }   
        
        infile2->Close();
      
      }else{
      
        TTree* treeD_negX = (TTree*)infile->Get("SpaCEtree_bkwdDisp");
        TTree* treeE_negX = (TTree*)infile->Get("SpaCEtree");
               
        TTree* treeD_posX = (TTree*)infile->Get("SpaCEtree_bkwdDisp");
        TTree* treeE_posX = (TTree*)infile->Get("SpaCEtree");
        
        std::vector<TH3F*> temp = Build_TH3(treeD_posX, treeE_posX, treeD_negX,treeE_negX, "x_reco", "y_reco", "z_reco", "bkwd");        
        for (size_t ii = 0; ii<temp.size(); ii++){
        	CalSCEhistograms.at(ii) = (TH3F*)temp.at(ii)->Clone(TString::Format("%s",temp.at(ii)->GetName()));
        	CalSCEhistograms.at(ii)->SetDirectory(0);
        }   
      
      }
      
   } else { std::cout << "No space charge representation type chosen." << std::endl;} 
    
    infile->Close();
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
//bool spacecharge::SpaceChargeProtoDUNE::EnableCorrSCE() const
//{
//  return fEnableCorrSCE;
//}

/// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceChargeProtoDUNE::EnableCalSpatialSCE() const
{
  return fEnableCalSpatialSCE;
}

/// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceChargeProtoDUNE::EnableCalEfieldSCE() const
{
  return fEnableCalEfieldSCE;
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
  
  if ((fRepresentationType=="Voxelized_TH3") || (fRepresentationType == "Splines_TH3")){
    if (point.X() > 0.) {
      thePosOffsets = GetOffsetsVoxel(point, SCEhistograms.at(0), SCEhistograms.at(1), SCEhistograms.at(2), 1, 2);
      thePosOffsets[0] = -1.0*thePosOffsets[0];
    } else {
      thePosOffsets = GetOffsetsVoxel(point, SCEhistograms.at(6), SCEhistograms.at(7), SCEhistograms.at(8), 1, 1);
      thePosOffsets[0] = -1.0*thePosOffsets[0];
    }
      
  }else if (fRepresentationType == "Voxelized"){
    if (point.X() > 0.) thePosOffsets = GetOffsetsVoxel(point, SCEhistograms.at(0), SCEhistograms.at(1), SCEhistograms.at(2), 1, 2);
    else thePosOffsets = GetOffsetsVoxel(point, SCEhistograms.at(6), SCEhistograms.at(7), SCEhistograms.at(8), 1, 1);
      
  }else if(fRepresentationType == "Parametric") thePosOffsets = GetPosOffsetsParametric(point.X(), point.Y(), point.Z());
  else thePosOffsets.resize(3,0.0); 
 
  geo::Point_t pafteroffset(point.X()+thePosOffsets[0], point.Y()+thePosOffsets[1], point.Z()+thePosOffsets[2]);
  geo::Vector_t edoffset = ElectronDiverterPosOffsets(pafteroffset);
  thePosOffsets[0] += edoffset.X();
  thePosOffsets[1] += edoffset.Y();
  thePosOffsets[2] += edoffset.Z();

  return { thePosOffsets[0], thePosOffsets[1], thePosOffsets[2] };
}

//----------------------------------------------------------------------------
/// Primary working method of service that provides position offsets to be
/// used in calibration of space charge
geo::Vector_t spacecharge::SpaceChargeProtoDUNE::GetCalPosOffsets(geo::Point_t const& tmp_point, int const& TPCid) const
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
  
  if ((fRepresentationType == "Voxelized_TH3") || (fRepresentationType == "Splines_TH3")){
    if ((TPCid == 2 || TPCid == 6 || TPCid == 10)&&point.X()>-20.){
      if (point.X()<0.) point = {0.00001, point.Y(), point.Z()};
      thePosOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(0), CalSCEhistograms.at(1), CalSCEhistograms.at(2), 2, 2);
      thePosOffsets[0] = -1.0*thePosOffsets[0];
    } else if((TPCid == 1 || TPCid == 5 || TPCid == 9)&&point.X()<20.) {
    	if (point.X()>0.) point= {-0.00001, point.Y(), point.Z()};
      thePosOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(6), CalSCEhistograms.at(7), CalSCEhistograms.at(8), 2, 1);
      thePosOffsets[0] = -1.0*thePosOffsets[0];
    } else thePosOffsets = {0., 0., 0.};
    
  } else if (fRepresentationType=="Voxelized"){
    if ((TPCid == 2 || TPCid == 6 || TPCid == 10)&&point.X()>-20.){
      if (point.X()<0.) point = {0.00001, point.Y(), point.Z()};
      thePosOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(0), CalSCEhistograms.at(1), CalSCEhistograms.at(2), 2, 2);
      thePosOffsets[0] = -1.0*thePosOffsets[0];
    } else if((TPCid == 1 || TPCid == 5 || TPCid == 9)&&point.X()<20.) {
    	if (point.X()>0.) point= {-0.00001, point.Y(), point.Z()};
      thePosOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(6), CalSCEhistograms.at(7), CalSCEhistograms.at(8), 2, 1);
    } else thePosOffsets = {0., 0., 0.};
      
  } else thePosOffsets.resize(3,0.0);
  
  return { thePosOffsets[0], thePosOffsets[1], thePosOffsets[2] };
}

//----------------------------------------------------------------------------
/// Provides position offsets using voxelized interpolation
std::vector<double> spacecharge::SpaceChargeProtoDUNE::GetOffsetsVoxel(geo::Point_t const& point, TH3F* hX, TH3F* hY, TH3F* hZ, int maptype, int driftvol) const
{
  if (fRepresentationType == "Voxelized_TH3"){
  
    return {
      hX->Interpolate(point.X(),point.Y(),point.Z()),
      hY->Interpolate(point.X(),point.Y(),point.Z()),
      hZ->Interpolate(point.X(),point.Y(),point.Z())
    };
    
  } else if (fRepresentationType == "Splines_TH3"){

    return {
      InterpolateSplines(hX,point.X(),point.Y(),point.Z(),1,maptype,driftvol),
      InterpolateSplines(hY,point.X(),point.Y(),point.Z(),2,maptype,driftvol),
      InterpolateSplines(hZ,point.X(),point.Y(),point.Z(),3,maptype,driftvol)
    };

  } else {
    double xnew = TransformX(point.X());
    double ynew = TransformY(point.Y());
    double znew = TransformZ(point.Z());
  
    return {
      hX->Interpolate(xnew,ynew,znew),
      hY->Interpolate(xnew,ynew,znew),
      hZ->Interpolate(xnew,ynew,znew)
    };
  }
}

//----------------------------------------------------------------------------
/// Build 3d histograms for voxelized interpolation
std::vector<TH3F*> spacecharge::SpaceChargeProtoDUNE::Build_TH3(TTree* tree_pos, TTree* eTree_pos, TTree* tree_neg, TTree* eTree_neg, std::string xvar, std::string yvar, std::string zvar, std::string posLeaf) const
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
  
  //initialized histograms for Dx, Dy, Dz, and electric field (pos x)
  TH3F* hDx_pos = new TH3F("hDx_pos", "", numDivisions_x+1, -0.5*cell_size, Lx+0.5*cell_size, numDivisions_y+1 ,-0.5*cell_size, Ly+0.5*cell_size, numDivisions_z+1, -0.5*cell_size, Lz+0.5*cell_size);
  TH3F* hDy_pos = new TH3F("hDy_pos", "", numDivisions_x+1, -0.5*cell_size, Lx+0.5*cell_size, numDivisions_y+1, -0.5*cell_size, Ly+0.5*cell_size, numDivisions_z+1, -0.5*cell_size, Lz+0.5*cell_size);
  TH3F* hDz_pos = new TH3F("hDz_pos", "", numDivisions_x+1, -0.5*cell_size, Lx+0.5*cell_size, numDivisions_y+1, -0.5*cell_size, Ly+0.5*cell_size, numDivisions_z+1, -0.5*cell_size, Lz+0.5*cell_size);
  
  TH3F* hEx_pos = new TH3F("hEx_pos", "", E_numDivisions_x+1, -0.5*E_cell_size, Lx+0.5*E_cell_size, E_numDivisions_y+1, -0.5*E_cell_size, Ly+0.5*E_cell_size, E_numDivisions_z+1, -0.5*E_cell_size, Lz+0.5*E_cell_size);
  TH3F* hEy_pos = new TH3F("hEy_pos", "", E_numDivisions_x+1, -0.5*E_cell_size, Lx+0.5*E_cell_size, E_numDivisions_y+1, -0.5*E_cell_size, Ly+0.5*E_cell_size, E_numDivisions_z+1, -0.5*E_cell_size, Lz+0.5*E_cell_size);
  TH3F* hEz_pos = new TH3F("hEz_pos", "", E_numDivisions_x+1, -0.5*E_cell_size, Lx+0.5*E_cell_size, E_numDivisions_y+1, -0.5*E_cell_size, Ly+0.5*E_cell_size, E_numDivisions_z+1, -0.5*E_cell_size, Lz+0.5*E_cell_size);
  
  //initialized histograms for Dx, Dy, Dz, and electric field (neg x)
  TH3F* hDx_neg = new TH3F("hDx_neg", "", numDivisions_x+1, -0.5*cell_size, Lx+0.5*cell_size, numDivisions_y+1 ,-0.5*cell_size, Ly+0.5*cell_size, numDivisions_z+1, -0.5*cell_size, Lz+0.5*cell_size);
  TH3F* hDy_neg = new TH3F("hDy_neg", "", numDivisions_x+1, -0.5*cell_size, Lx+0.5*cell_size, numDivisions_y+1, -0.5*cell_size, Ly+0.5*cell_size, numDivisions_z+1, -0.5*cell_size, Lz+0.5*cell_size);
  TH3F* hDz_neg = new TH3F("hDz_neg", "", numDivisions_x+1, -0.5*cell_size, Lx+0.5*cell_size, numDivisions_y+1, -0.5*cell_size, Ly+0.5*cell_size, numDivisions_z+1, -0.5*cell_size, Lz+0.5*cell_size);
  
  TH3F* hEx_neg = new TH3F("hEx_neg", "", E_numDivisions_x+1, -0.5*E_cell_size, Lx+0.5*E_cell_size, E_numDivisions_y+1, -0.5*E_cell_size, Ly+0.5*E_cell_size, E_numDivisions_z+1, -0.5*E_cell_size, Lz+0.5*E_cell_size);
  TH3F* hEy_neg = new TH3F("hEy_neg", "", E_numDivisions_x+1, -0.5*E_cell_size, Lx+0.5*E_cell_size, E_numDivisions_y+1, -0.5*E_cell_size, Ly+0.5*E_cell_size, E_numDivisions_z+1, -0.5*E_cell_size, Lz+0.5*E_cell_size);
  TH3F* hEz_neg = new TH3F("hEz_neg", "", E_numDivisions_x+1, -0.5*E_cell_size, Lx+0.5*E_cell_size, E_numDivisions_y+1, -0.5*E_cell_size, Ly+0.5*E_cell_size, E_numDivisions_z+1, -0.5*E_cell_size, Lz+0.5*E_cell_size);
 
  //For each event, read the tree and fill each histogram (pos x)
  for (int ii = 0; ii<tree_pos->GetEntries(); ii++){

    //Read the trees
    tree_pos->GetEntry(ii);
    Double_t x = tree_pos->GetBranch(xvar.c_str())->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t y = tree_pos->GetBranch(yvar.c_str())->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t z = tree_pos->GetBranch(zvar.c_str())->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t dx = tree_pos->GetBranch("Dx")->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t dy = tree_pos->GetBranch("Dy")->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t dz = tree_pos->GetBranch("Dz")->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    
    hDx_pos->Fill(x,y,z,100.0*dx);
    hDy_pos->Fill(x,y,z,100.0*dy);
    hDz_pos->Fill(x,y,z,100.0*dz);
  }
  
  for(int ii = 0; ii<eTree_pos->GetEntries(); ii++){
		
    eTree_pos->GetEntry(ii);
    Double_t x = eTree_pos->GetBranch("xpoint")->GetLeaf("data")->GetValue();
    Double_t y = eTree_pos->GetBranch("ypoint")->GetLeaf("data")->GetValue();
    Double_t z = eTree_pos->GetBranch("zpoint")->GetLeaf("data")->GetValue();
    Double_t Ex = eTree_pos->GetBranch("Ex")->GetLeaf("data")->GetValue() / (100000.0*fEfield);
    Double_t Ey = eTree_pos->GetBranch("Ey")->GetLeaf("data")->GetValue() / (100000.0*fEfield);
    Double_t Ez = eTree_pos->GetBranch("Ez")->GetLeaf("data")->GetValue() / (100000.0*fEfield);
   
    //Fill the histograms		
    hEx_pos->Fill(x,y,z,Ex);
    hEy_pos->Fill(x,y,z,Ey);
    hEz_pos->Fill(x,y,z,Ez);
  }
  
  //For each event, read the tree and fill each histogram (neg x)
  for (int ii = 0; ii<tree_neg->GetEntries(); ii++){

    //Read the trees
    tree_neg->GetEntry(ii);
    Double_t x = tree_neg->GetBranch(xvar.c_str())->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t y = tree_neg->GetBranch(yvar.c_str())->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t z = tree_neg->GetBranch(zvar.c_str())->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t dx = tree_neg->GetBranch("Dx")->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t dy = tree_neg->GetBranch("Dy")->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    Double_t dz = tree_neg->GetBranch("Dz")->GetLeaf(Form("data_%sDisp",posLeaf.c_str()))->GetValue();
    
    hDx_neg->Fill(x,y,z,100.0*dx);
    hDy_neg->Fill(x,y,z,100.0*dy);
    hDz_neg->Fill(x,y,z,100.0*dz);
  }
  
  for(int ii = 0; ii<eTree_neg->GetEntries(); ii++){
		
    eTree_neg->GetEntry(ii);
    Double_t x = eTree_neg->GetBranch("xpoint")->GetLeaf("data")->GetValue();
    Double_t y = eTree_neg->GetBranch("ypoint")->GetLeaf("data")->GetValue();
    Double_t z = eTree_neg->GetBranch("zpoint")->GetLeaf("data")->GetValue();
    Double_t Ex = eTree_neg->GetBranch("Ex")->GetLeaf("data")->GetValue() / (100000.0*fEfield);
    Double_t Ey = eTree_neg->GetBranch("Ey")->GetLeaf("data")->GetValue() / (100000.0*fEfield);
    Double_t Ez = eTree_neg->GetBranch("Ez")->GetLeaf("data")->GetValue() / (100000.0*fEfield);
   
    //Fill the histograms		
    hEx_neg->Fill(x,y,z,Ex);
    hEy_neg->Fill(x,y,z,Ey);
    hEz_neg->Fill(x,y,z,Ez);
  }
  
  return {hDx_pos, hDy_pos, hDz_pos, hEx_pos, hEy_pos, hEz_pos, hDx_neg, hDy_neg, hDz_neg, hEx_neg, hEy_neg, hEz_neg};

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
  
  if ((fRepresentationType=="Voxelized_TH3") || (fRepresentationType == "Splines_TH3")){
    if (point.X() > 0.) theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms.at(3), SCEhistograms.at(4), SCEhistograms.at(5), 3, 2);
    else theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms.at(9), SCEhistograms.at(10), SCEhistograms.at(11), 3, 1);
    theEfieldOffsets[0] = -1.0*theEfieldOffsets[0];
    theEfieldOffsets[1] = -1.0*theEfieldOffsets[1];
    theEfieldOffsets[2] = -1.0*theEfieldOffsets[2];
  }else if (fRepresentationType == "Voxelized"){
    if (point.X() > 0.) theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms.at(3), SCEhistograms.at(4), SCEhistograms.at(5), 3, 2);
    else theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms.at(9), SCEhistograms.at(10), SCEhistograms.at(11), 3, 1);
  }else if(fRepresentationType == "Parametric") theEfieldOffsets = GetEfieldOffsetsParametric(point.X(), point.Y(), point.Z());
  else theEfieldOffsets.resize(3,0.0);
    
   return { -theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2] };
}
//----------------------------------------------------------------------------
/// Primary working method of service that provides E field offsets to be
/// used in charge/light yield calculation (e.g.) for calibration
geo::Vector_t spacecharge::SpaceChargeProtoDUNE::GetCalEfieldOffsets(geo::Point_t const& tmp_point, int const& TPCid) const
{ 
  std::vector<double> theEfieldOffsets;
  geo::Point_t point = tmp_point;
  if(IsTooFarFromBoundaries(point)) {
    theEfieldOffsets.resize(3,0.0);
    return { -theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2] };
  }
  if(!IsInsideBoundaries(point)&&!IsTooFarFromBoundaries(point)) point = PretendAtBoundary(point);
  
  if ((fRepresentationType == "Voxelized_TH3") || (fRepresentationType == "Splines_TH3")){
    if ((TPCid == 2 || TPCid == 6 || TPCid == 10)&&point.X()>-20.){
      if (point.X()<0.) point = {0.00001, point.Y(), point.Z()};
      theEfieldOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(3), CalSCEhistograms.at(4), CalSCEhistograms.at(5), 3, 2);
    }else if ((TPCid == 1 || TPCid == 5 || TPCid == 9)&&point.X()<20.){
      if (point.X()>0.) point = {-0.00001, point.Y(), point.Z()};
      theEfieldOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(9), CalSCEhistograms.at(10), CalSCEhistograms.at(11), 3, 1);
    } else theEfieldOffsets = {0., 0., 0.};
    theEfieldOffsets[0] = -1.0*theEfieldOffsets[0];
    theEfieldOffsets[1] = -1.0*theEfieldOffsets[1];
    theEfieldOffsets[2] = -1.0*theEfieldOffsets[2];
  }else if (fRepresentationType == "Voxelized"){
    if ((TPCid == 2 || TPCid == 6 || TPCid == 10)&&point.X()>-20.){
      if (point.X()<0.) point = {0.00001, point.Y(), point.Z()};
      theEfieldOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(3), CalSCEhistograms.at(4), CalSCEhistograms.at(5), 3, 2);
    }else if ((TPCid == 1 || TPCid == 5 || TPCid == 9)&&point.X()<20.){
      if (point.X()>0.) point = {-0.00001, point.Y(), point.Z()};
      theEfieldOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(9), CalSCEhistograms.at(10), CalSCEhistograms.at(11), 3, 1);
    } else theEfieldOffsets = {0., 0., 0.};
  }else
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
  if((fRepresentationType=="Voxelized_TH3") || (fRepresentationType == "Splines_TH3")){
  	return !(
         (TMath::Abs(point.X()) <= 0.0) || (TMath::Abs(point.X()) >= 360.0)
      || (point.Y()             <= 5.2) || (point.Y()             >= 604.0)
      || (point.Z()             <= -0.5) || (point.Z()             >= 695.3)
    );
  } else{
  	return !(
         (TMath::Abs(point.X()) <=  0.0) || (TMath::Abs(point.X()) >= 360.0)
      || (point.Y()             <= -0.2) || (point.Y()             >= 607.8)
      || (point.Z()             <= -0.8) || (point.Z()             >= 696.2)
    );
  }
} 
  
bool spacecharge::SpaceChargeProtoDUNE::IsTooFarFromBoundaries(geo::Point_t const& point) const
{
  if((fRepresentationType=="Voxelized_TH3") || (fRepresentationType == "Splines_TH3")){
    return (
         (TMath::Abs(point.X()) < -20.0) || (TMath::Abs(point.X())  >= 360.0)
      || (point.Y()             < -14.8) || (point.Y()              >  624.0)
      || (point.Z()             < -20.5) || (point.Z()              >  715.3)
    );
  } else {
    return (
         (TMath::Abs(point.X()) < -20.0) || (TMath::Abs(point.X())  >= 360.0)
      || (point.Y()             < -20.2) || (point.Y()              >  627.8)
      || (point.Z()             < -20.8) || (point.Z()              >  716.2)
    );
  }
}

geo::Point_t spacecharge::SpaceChargeProtoDUNE::PretendAtBoundary(geo::Point_t const& point) const
{
  double x = point.X(), y = point.Y(), z = point.Z();
  
  if((fRepresentationType=="Voxelized_TH3") || (fRepresentationType == "Splines_TH3")){ 
  
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
    
  }
  return {x, y, z};
}

//----------------------------------------------------------------------------
/// Create one spline for later use in SCE map interpolation
TSpline3* spacecharge::SpaceChargeProtoDUNE::MakeSpline(TH3F* spline_hist, int dim1, int dim2_bin, int dim3_bin, int maptype, int driftvol) const
{
  TSpline3 *spline = 0;
  
  std::vector<double> a, b;
  if (dim1 == 1) {
    for (int x = 1; x <= spline_hist->GetNbinsX(); ++x) {
      a.push_back(spline_hist->GetXaxis()->GetBinCenter(x));
      b.push_back(spline_hist->GetBinContent(x, dim2_bin, dim3_bin));
    }
  }
  else if (dim1 == 2) {
    for(int y = 1; y <= spline_hist->GetNbinsY(); ++y) {
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
    cet::exception("SpaceChargeProtoDUNE::MakeSpline") << "Unkown dimension " << dim1 << "\n";
  }

  if(maptype == 1)
  {
    spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin,
                          dim3_bin, maptype, driftvol), &a[0], &b[0], a.size(),
                          "b2e2", 0, 0);
    spline->SetName(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin, dim3_bin,
                    maptype, driftvol));
  }
  else if(maptype == 2)
  {
    spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin,
                          dim3_bin, maptype, driftvol), &a[0], &b[0], a.size(),
                          "b2e2", 0, 0);
    spline->SetName(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin, dim3_bin,
                    maptype, driftvol));
  }
  else if(maptype == 3)
  {
    spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin,
                          dim3_bin, maptype, driftvol), &a[0], &b[0], a.size(),
                          "b2e2", 0, 0);
    spline->SetName(Form("spline_%d_%d_%d_%d_%d", dim1, dim2_bin, dim3_bin,
                    maptype, driftvol));
  }

  /*if(dim1 == 1)
  {
    double a[19];
    double b[19];
    for(int x = 1; x <= 19; x++)
    {
      a[x-1] = spline_hist->GetXaxis()->GetBinCenter(x);
      b[x-1] = spline_hist->GetBinContent(x,dim2_bin,dim3_bin);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }
  else if(dim1 == 2)
  {
    double a[31];
    double b[31];
    for(int y = 1; y <= 31; y++)
    {
      a[y-1] = spline_hist->GetYaxis()->GetBinCenter(y);
      b[y-1] = spline_hist->GetBinContent(dim2_bin,y,dim3_bin);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }
  else if(dim1 == 3)
  {
    double a[37];
    double b[37];
    for(int z = 1; z <= 37; z++)
    {
      a[z-1] = spline_hist->GetZaxis()->GetBinCenter(z);
      b[z-1] = spline_hist->GetBinContent(dim2_bin,dim3_bin,z);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }*/

  return spline;
}

//----------------------------------------------------------------------------
/// Interpolate given SCE map using splines
double spacecharge::SpaceChargeProtoDUNE::InterpolateSplines(TH3F* interp_hist, double xVal, double yVal, double zVal, int dim, int maptype, int driftvol) const
{

  //std::cout << "Interpolating " << interp_hist->GetName() << std::endl;
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
  if(bin_x <= 1)
  {
    low_x = 1;
    high_x = 2;
  }
  else if(bin_x >= max_x)
  {
    low_x = max_x-1;
    high_x = max_x;
  }
  else if(xVal > bincenter_x)
  {
    low_x = bin_x;
    high_x = bin_x+1;
  }
  else
  {
    low_x = bin_x-1;
    high_x = bin_x;
  }

  int low_y;
  int high_y;
  if(bin_y <= 1)
  {
    low_y = 1;
    high_y = 2;
  }
  else if(bin_y >= max_y)
  {
    low_y = max_y-1;
    high_y = max_y;
  }
  else if(yVal > bincenter_y)
  {
    low_y = bin_y;
    high_y = bin_y+1;
  }
  else
  {
    low_y = bin_y-1;
    high_y = bin_y;
  }

  int low_z;
  int high_z;
  if(bin_z <= 1)
  {
    low_z = 1;
    high_z = 2;
  }
  else if(bin_z >= max_z)
  {
    low_z = max_z-1;
    high_z = max_z;
  }
  else if(zVal > bincenter_z)
  {
    low_z = bin_z;
    high_z = bin_z+1;
  }
  else
  {
    low_z = bin_z-1;
    high_z = bin_z;
  }

  double interp_val = 0.0;
  
  if(dim == 1)
  {
    double a_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
    double a_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

    double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
    double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dx_fwd_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_fwd_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_fwd_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_fwd_neg[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dx_bkwd_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_bkwd_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_bkwd_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_bkwd_neg[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEx_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dEx_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dEx_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dEx_neg[high_y-1][high_z-1]->Eval(xVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dx_fwd_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_fwd_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_fwd_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_fwd_pos[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dx_bkwd_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_bkwd_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_bkwd_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_bkwd_pos[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEx_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dEx_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dEx_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dEx_pos[high_y-1][high_z-1]->Eval(xVal);
      }
    }

    interp_val = (f_11*(a_2-yVal)*(b_2-zVal) + f_21*(yVal-a_1)*(b_2-zVal) + f_12*(a_2-yVal)*(zVal-b_1) + f_22*(yVal-a_1)*(zVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }
  else if(dim == 2)
  {
    double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
    double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

    double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
    double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dy_fwd_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_fwd_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_fwd_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_fwd_neg[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dy_bkwd_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_bkwd_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_bkwd_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_bkwd_neg[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEy_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dEy_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dEy_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dEy_neg[high_x-1][high_z-1]->Eval(yVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dy_fwd_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_fwd_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_fwd_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_fwd_pos[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dy_bkwd_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_bkwd_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_bkwd_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_bkwd_pos[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEy_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dEy_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dEy_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dEy_pos[high_x-1][high_z-1]->Eval(yVal);
      }
    }

    interp_val = (f_11*(a_2-xVal)*(b_2-zVal) + f_21*(xVal-a_1)*(b_2-zVal) + f_12*(a_2-xVal)*(zVal-b_1) + f_22*(xVal-a_1)*(zVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }
  else if(dim == 3)
  {
    double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
    double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

    double b_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
    double b_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dz_fwd_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_fwd_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_fwd_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_fwd_neg[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dz_bkwd_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_bkwd_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_bkwd_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_bkwd_neg[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEz_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dEz_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dEz_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dEz_neg[high_x-1][high_y-1]->Eval(zVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dz_fwd_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_fwd_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_fwd_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_fwd_pos[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dz_bkwd_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_bkwd_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_bkwd_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_bkwd_pos[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEz_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dEz_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dEz_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dEz_pos[high_x-1][high_y-1]->Eval(zVal);
      }
    }

    interp_val = (f_11*(a_2-xVal)*(b_2-yVal) + f_21*(xVal-a_1)*(b_2-yVal) + f_12*(a_2-xVal)*(yVal-b_1) + f_22*(xVal-a_1)*(yVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }

  return interp_val;
}


geo::Vector_t spacecharge::SpaceChargeProtoDUNE::ElectronDiverterPosOffsets(geo::Point_t const& point) const
{
  double z = point.Z();
  double offset[3] = {0,0,0};

  for (size_t i=0; i<fEnableElectronDiverterDistortions.size(); ++i)
    {
      if (!fEnableElectronDiverterDistortions.at(i)) continue;
      if (point.X()>0) continue;
      if (z>fEDChargeLossZLow.at(i) && z<fEDChargeLossZHigh.at(i))
        {
          offset[0] = 2E9;
          offset[1] = 2E9;
          offset[2] = 2E9;
        }
      else
	{
	  double zdiff = z - fEDZCenter.at(i);
	  double zexp = TMath::Exp( -TMath::Sq(zdiff/fEDs.at(i)) );
	  offset[2] += fEDBZPosOffs.at(i) * zdiff * zexp;

	  // the timing offsets need to be computed after the z shift
	  double zdiffc = zdiff + offset[2];
	  double zexpc = TMath::Exp( -TMath::Sq(zdiffc/fEDs.at(i)) );
	  offset[0] += fEDAXPosOffs.at(i) * zexpc;
	}
    }
  return {offset[0], offset[1], offset[2]};
}
