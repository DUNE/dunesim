////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEBeamTPCMatching
// Plugin Type: analyzer (art v2_07_03)
// File:        ProtoDUNEBeamTPCMatching_module.cc
//
// Generated at Mon Sep  4 06:55:33 2017 by Leigh Whitehead using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////
//// Modified by Pablo and Leigh H. Howard, 
//Smear important variables of beam monitors, mimic Cherenkov monitors response
//and store some histograms through art utilities
//// pablo.fer@cern.ch
//// July 2018
///////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
//Histogram utilities
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include <fstream>
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"


#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"
#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"


namespace sim {
  class ProtoDUNEBeamTPCMatching;
}


class sim::ProtoDUNEBeamTPCMatching : public art::EDAnalyzer {
public:

  explicit ProtoDUNEBeamTPCMatching(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ProtoDUNEBeamTPCMatching(ProtoDUNEBeamTPCMatching const &) = delete;
  ProtoDUNEBeamTPCMatching(ProtoDUNEBeamTPCMatching &&) = delete;
  ProtoDUNEBeamTPCMatching & operator = (ProtoDUNEBeamTPCMatching const &) = delete;
  ProtoDUNEBeamTPCMatching & operator = (ProtoDUNEBeamTPCMatching &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;


  // Required functions.
  void analyze(art::Event const & e) override;




private:
float fBeamTPC_theta;
float fBeamTPC_phi;
float beam_true_bprof4_pos[4], beam_true_tof,    beam_true_dir[3];
float beam_smeared_bprof4_pos[4], beam_smeared_tof, beam_smeared_dir[3];

float tpc_true_startpos[4], tpc_true_dir[3];
//float tpc_true_E, tpc_true_mom, dummy;
float tpc_reco_startpos[4], tpc_reco_dir[3];
TVector3 tpc_pos, tpc_dir;


std::string fpandoraTrack;

TTree *beamtpc;
/*TH1F *fT0FHist_zoom;
TH1F *fT0FsHist_zoom;
TH1F *fUxHist_zoom;
TH1F *fUxsHist_zoom;
TH1F *fUyHist_zoom;
TH1F *fUysHist_zoom;
TH1F *fUzHist_zoom;
TH1F *fUzsHist_zoom;
TH1F *fT0FHist;
TH1F *fT0FsHist;
TH1F *fUxHist;
TH1F *fUxsHist;
TH1F *fUyHist;
TH1F *fUysHist;
TH1F *fUzHist;
TH1F *fUzsHist;
*/};


sim::ProtoDUNEBeamTPCMatching::ProtoDUNEBeamTPCMatching(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{

}

void sim::ProtoDUNEBeamTPCMatching::beginJob()
{
// Declaring histograms
  art::ServiceHandle<art::TFileService> tfs;
/*  fT0FHist_zoom = tfs->make<TH1F>("TOF_TRUE_zoom","T0F (ns)",100,95,96);
  fT0FsHist_zoom = tfs->make<TH1F>("TOF_SMEARED_zoom","T0F (ns)",100,95,96);
  fUxHist_zoom = tfs->make<TH1F>("X_DIRECTION_TRUE_zoom","Ux",100,-0.05,0.05);
  fUyHist_zoom = tfs->make<TH1F>("Y_DIRECTION_TRUE_zoom","Uy",100,-0.05,0.05);
  fUzHist_zoom = tfs->make<TH1F>("Z_DIRECTION_TRUE_zoom","Uz",100,-1,-0.95);
  fUxsHist_zoom = tfs->make<TH1F>("X_DIRECTION_SMEARED_zoom","Ux",100,-0.05,0.05);
  fUysHist_zoom = tfs->make<TH1F>("Y_DIRECTION_SMEARED_zoom","Uy",100,-0.05,0.05);
  fUzsHist_zoom = tfs->make<TH1F>("Z_DIRECTION_SMEARED_zoom","Uz",100,-1,-0.95);
  fT0FHist = tfs->make<TH1F>("TOF_TRUE","T0F (ns)",100,90,100);
  fT0FsHist = tfs->make<TH1F>("TOF_SMEARED","T0F (ns)",100,90,100);
  fUxHist = tfs->make<TH1F>("X_DIRECTION_TRUE","Ux",100,-1,1);
  fUyHist = tfs->make<TH1F>("Y_DIRECTION_TRUE","Uy",100,-1,1);
  fUzHist = tfs->make<TH1F>("Z_DIRECTION_TRUE","Uz",100,-1,1);
  fUxsHist = tfs->make<TH1F>("X_DIRECTION_SMEARED","Ux",100,-1,1);
  fUysHist = tfs->make<TH1F>("Y_DIRECTION_SMEARED","Uy",100,-1,1);
  fUzsHist = tfs->make<TH1F>("Z_DIRECTION_SMEARED","Uz",100,-1,1);
*/
  beamtpc = tfs->make<TTree>("beam-tpc", "Beam + TPC matiching Tree");
  beamtpc->Branch("beam_true_tof",      &beam_true_tof,          "beam_true_tof/F");
  beamtpc->Branch("beam_true_dir",      &beam_true_dir,          "beam_true_dir[3]/F");
  beamtpc->Branch("beam_true_bprof4_pos",&beam_true_bprof4_pos,  "beam_true_bprof4_pos[4]/F");
  beamtpc->Branch("beam_smeared_tof",   &beam_smeared_tof,       "beam_smeared_tof/F");
  beamtpc->Branch("beam_smeared_dir",   &beam_smeared_dir,       "beam_smeared_dir[3]/F");
  beamtpc->Branch("beam_smeared_bprof4_pos",&beam_smeared_bprof4_pos,  "beam_true_bprof4_pos[4]/F");

  beamtpc->Branch("tpc_true_startpos",  &tpc_true_startpos,      "tpc_true_startpos[4]/F");
  beamtpc->Branch("tpc_true_dir",       &tpc_true_dir,           "tpc_true_dir[3]/F");
  beamtpc->Branch("tpc_reco_startpos",  &tpc_reco_startpos,      "tpc_reco_startpos[4]/F");
  beamtpc->Branch("tpc_reco_dir",       &tpc_reco_dir,           "tpc_reco_dir[3]/F");
}

//--------------------------------------------------------------------
void sim::ProtoDUNEBeamTPCMatching::reconfigure(fhicl::ParameterSet const & pset)
{
  // Rotation angles from last beam section to TPC cordinates
//  fBeamTPC_theta = pset.get<float>("BeamTPC_theta");
//  fBeamTPC_phi   = pset.get<float>("BeamTPC_phi");
  // The name of the module that produced the tracks
  fpandoraTrack  = pset.get<std::string>("TrackLabel");
    } // reconfigure


void sim::ProtoDUNEBeamTPCMatching::analyze(art::Event const & evt)
{
  fBeamTPC_theta=0.27361019;
  fBeamTPC_phi=3.9651559;
//========Start of Beam part==========
  // Get the reconstructed tracks
  auto beamsim = evt.getValidHandle<std::vector<sim::ProtoDUNEbeamsim> >("generator");
  const sim::ProtoDUNEbeamsim temp = (*beamsim)[0];
//  unsigned short nInst = temp.NInstruments();
// Computing true and smeared TOF (between TOF1 and TRIG2)
  sim::ProtoDUNEBeamInstrument tof1 = temp.GetInstrument("TOF1");
  sim::ProtoDUNEBeamInstrument trig2 = temp.GetInstrument("TRIG2");
  beam_true_tof = trig2.GetT() - tof1.GetT();
  beam_smeared_tof = trig2.GetSmearedVar1() -tof1.GetSmearedVar1();
//  std::cout << "TOF " << tof << std::endl;
//  std::cout << "TOFS " << tof_smeared << std::endl;

  sim::ProtoDUNEBeamInstrument bprofext = temp.GetInstrument("BPROFEXT");
  sim::ProtoDUNEBeamInstrument bprof4 = temp.GetInstrument("BPROF4");

  beam_true_bprof4_pos[0]    = bprof4.GetX() - 0.0;//?
  beam_true_bprof4_pos[1]    = bprof4.GetY() - 0.0;//?
  beam_true_bprof4_pos[2]    = bprof4.GetZ() - 717242.5;
  beam_true_bprof4_pos[3]    = bprof4.GetT();
  beam_smeared_bprof4_pos[0] = bprof4.GetSmearedVar1() - 0.0;//?
  beam_smeared_bprof4_pos[1] = bprof4.GetSmearedVar2() - 0.0;//?
  beam_smeared_bprof4_pos[2] = bprof4.GetZ() - 717242.5;
  beam_smeared_bprof4_pos[3] = bprof4.GetT();

  double dummy;
  dummy = pow(pow(bprofext.GetX() -bprof4.GetX(),2)+pow(bprofext.GetY() -bprof4.GetY(),2)+pow(bprofext.GetZ() -bprof4.GetZ(),2),0.5);
  beam_true_dir[0] = (bprofext.GetX() -bprof4.GetX())/dummy;
  beam_true_dir[1] = (bprofext.GetY() -bprof4.GetY())/dummy;
  beam_true_dir[2] = (bprofext.GetZ() -bprof4.GetZ())/dummy;
  beam_true_dir[0] = cos(fBeamTPC_phi)*beam_true_dir[0] + sin(fBeamTPC_phi)*cos(fBeamTPC_theta)*beam_true_dir[1] - sin(fBeamTPC_phi)*sin(fBeamTPC_theta)*beam_true_dir[2];
  beam_true_dir[1] = -sin(fBeamTPC_phi)*beam_true_dir[0] + cos(fBeamTPC_phi)*cos(fBeamTPC_theta)*beam_true_dir[1] + sin(fBeamTPC_theta)*cos(fBeamTPC_phi)*beam_true_dir[2];
  beam_true_dir[2] = -sin(fBeamTPC_theta)*beam_true_dir[1] + cos(fBeamTPC_theta)*beam_true_dir[2];
//  std::cout << dir[2] << " " << sin(fBeamTPC_phi)*sin(fBeamTPC_theta)*dir[0] << " + " << - cos(fBeamTPC_phi)*sin(fBeamTPC_theta)*dir[1] << "+ " << cos(fBeamTPC_theta)*dir[2] << std::endl;
  dummy = pow(pow(bprofext.GetSmearedVar1() -bprof4.GetSmearedVar1(),2)+pow(bprofext.GetSmearedVar2() -bprof4.GetSmearedVar2(),2)+pow(bprofext.GetZ() -bprof4.GetZ(),2),0.5);
  beam_smeared_dir[0] = (bprofext.GetSmearedVar1() -bprof4.GetSmearedVar1())/dummy;
  beam_smeared_dir[1] = (bprofext.GetSmearedVar2() -bprof4.GetSmearedVar2())/dummy;
  beam_smeared_dir[2] = (bprofext.GetZ() -bprof4.GetZ())/dummy;
  beam_smeared_dir[0] = cos(fBeamTPC_phi)*beam_true_dir[0] + sin(fBeamTPC_phi)*cos(fBeamTPC_theta)*beam_true_dir[1] - sin(fBeamTPC_phi)*sin(fBeamTPC_theta)*beam_true_dir[2];
  beam_smeared_dir[1] = -sin(fBeamTPC_phi)*beam_true_dir[0] + cos(fBeamTPC_phi)*cos(fBeamTPC_theta)*beam_true_dir[1] + sin(fBeamTPC_theta)*cos(fBeamTPC_phi)*beam_true_dir[2];
  beam_smeared_dir[2] = -sin(fBeamTPC_theta)*beam_true_dir[1] + cos(fBeamTPC_theta)*beam_true_dir[2];
//  std::cout <<  fBeamTPC_phi << " " << fBeamTPC_theta << std::endl;

//  std::cout << "cos" << dir[0]*dir_smeared[0]+dir[1]*dir_smeared[1]+dir[2]*dir_smeared[2] << std::endl;
//  std::cout << bprofext.GetX() << "," << bprofext.GetSmearedVar1() << std::endl;
//  std::cout << dir[0] << "," << dir_smeared[0] << std::endl;

//Filling hisotgrams
/*  fT0FHist_zoom->Fill(beam_true_tof);
  fT0FsHist_zoom->Fill(beam_smeared_tof);
  fUxHist_zoom->Fill(beam_true_dir[0]);
  fUxsHist_zoom->Fill(beam_smeared_dir[0]);
  fUyHist_zoom->Fill(beam_true_dir[1]);
  fUysHist_zoom->Fill(beam_smeared_dir[1]);
  fUzHist_zoom->Fill(beam_true_dir[2]);
  fUzsHist_zoom->Fill(beam_smeared_dir[2]);
  fT0FHist->Fill(beam_true_tof);
  fT0FsHist->Fill(beam_smeared_tof);
  fUxHist->Fill(beam_true_dir[0]);
  fUxsHist->Fill(beam_smeared_dir[0]);
  fUyHist->Fill(beam_true_dir[1]);
  fUysHist->Fill(beam_smeared_dir[1]);
  fUzHist->Fill(beam_true_dir[2]);
  fUzsHist->Fill(beam_smeared_dir[2]);*/
//========End of Beam part==========
//
//========Start of TPC part==========
//Get tracks
/*  art::Handle< std::vector<simb::MCTruth> > mctruthhandle;
  std::vector< art::Ptr<simb::MCTruth> > mclist;
  art::Ptr<simb::MCTruth> mctruth;
  mctruth = mclist[0];
*/  
  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fpandoraTrack);
  protoana::ProtoDUNETruthUtils truthUtil;
  protoana::ProtoDUNETrackUtils trackUtil;
  dummy = 999999.;
  TVector3 pos, dir_start, dir_end;
  for(unsigned int t = 0; t < recoTracks->size(); ++t){
    const recob::Track thisTrack = (*recoTracks)[t];
    tpc_pos		= thisTrack.Vertex();
    tpc_dir		= thisTrack.VertexDirection();
    std::vector<anab::T0>  trackT0 = trackUtil.GetRecoTrackT0(thisTrack,evt,fpandoraTrack);
    tpc_reco_startpos[3] = trackT0[0].fTime;
    tpc_reco_startpos[0] = tpc_pos.X();
    tpc_reco_startpos[1] = tpc_pos.Y();
    tpc_reco_startpos[2] = tpc_pos.Z();
//    tpc_true_startpos = thisTrack->StartMomentumVector();
    if (dummy >= tpc_reco_startpos[3]){
      dummy = tpc_reco_startpos[3];
      const simb::MCParticle *trueMatch = truthUtil.GetMCParticleFromRecoTrack(thisTrack,evt,fpandoraTrack);
      bool hasTruth = (trueMatch != 0x0);
      if (hasTruth == true) {
//   	 std::cout << "hey" << std::endl;
         tpc_true_startpos[0]        = trueMatch->Vx(); 
	 tpc_true_startpos[1]        = trueMatch->Vy(); 
	 tpc_true_startpos[2]        = trueMatch->Vz(); 
	 tpc_true_startpos[3]        = trueMatch->T(); 
	 tpc_true_dir[0]             = trueMatch->Px() / trueMatch->P(); 
	 tpc_true_dir[1]             = trueMatch->Py() / trueMatch->P(); 
	 tpc_true_dir[2]             = trueMatch->Pz() / trueMatch->P(); 
  }
  }
  }


  

/*// Find the associations between tracks and T0
   const art::FindManyP<anab::T0> findTrackT0(trackHandle,evt,fpandoraTrack);
 
// Also look for cosmic tags so we can make a T0 plot for cosmic tagged events only
   const art::FindManyP<anab::CosmicTag> findCosmicTag(trackHandle,evt,fpandoraTrack);
*/
beamtpc->Fill();
}  




void sim::ProtoDUNEBeamTPCMatching::endJob()
{

}

DEFINE_ART_MODULE(sim::ProtoDUNEBeamTPCMatching)

