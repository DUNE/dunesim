////////////////////////////////////////////////////////////////////////
// Class:       PDSPDataDrivenBeam
// Plugin Type: producer (art v3_05_01)
// File:        PDSPDataDrivenBeam_module.cc
//
// Generated at Thu Aug  6 08:22:36 2020 by Jacob Calcutt using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <memory>
#include <map>

#include "TFile.h"
#include "THnSparse.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TH1D.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

class PDSPDataDrivenBeam;

class PDSPDataDrivenBeam : public art::EDProducer {
public:
  explicit PDSPDataDrivenBeam(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPDataDrivenBeam(PDSPDataDrivenBeam const&) = delete;
  PDSPDataDrivenBeam(PDSPDataDrivenBeam&&) = delete;
  PDSPDataDrivenBeam&
      operator=(PDSPDataDrivenBeam const&) = delete;
  PDSPDataDrivenBeam& operator=(PDSPDataDrivenBeam&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginJob() override;

private:
  
  TFile* fInputFile = 0x0;
  std::string fFileName;
  std::vector<std::string> fParticleTypes;
  int fNGenerate;
  int fSeed;
  TRandom3 fRNG;
  bool fVerbose;

  std::map<std::string, THnSparseD *> fPDFs;
  std::map<std::string, int> fParticleNums;
  double fTotalParticles = 0.;
  std::map<std::string, double> fParticleFracs;
  int fCurrentEvent = 0;

  std::vector<int> fRandPDG;
  std::vector<double> fRandMomentum;
  std::vector<double> fRandHUpstream, fRandVUpstream;
  std::vector<double> fRandHDownstream, fRandVDownstream;

  TTree * fOutputTree;
  int fOutputPDG;
  double fOutputMomentum;
  double fOutputHUpstream, fOutputVUpstream;
  double fOutputHDownstream, fOutputVDownstream;

  std::map<std::string, int> fNameToPDG = {
    {"Protons", 2212},
    {"Pions", 211},
    {"Electrons", -11},
  };

  TVector3 fBMBasisX = TVector3(1.,0.,0.);
  TVector3 fBMBasisY = TVector3(0.,1.,0.);
  TVector3 fBMBasisZ = TVector3(0.,0.,1.);

  double fRotateMonitorXZ, fRotateMonitorYZ;
  double fNP04FrontZ, fBeamX, fBeamY, fBeamZ;
  double fUpstreamZ, fDownstreamZ, fGeneratorZ;


  void Convert(double input_point[5], double minima[5],
               double maxima[5], double output_point[5]);

  void GenerateTrueEvent(simb::MCTruth &mcTruth);
  void SetPositionAndMomentum(TLorentzVector & position,
                              TLorentzVector & momentum);
  TVector3 ConvertProfCoordinates(double x, double y, double z, double zOffset);
  void BeamMonitorBasisVectors();
  void RotateMonitorVector(TVector3 & vec);
};


PDSPDataDrivenBeam::PDSPDataDrivenBeam(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fFileName(p.get<std::string>("FileName")),
    fNGenerate(p.get<int>("NGenerate")),
    fSeed(p.get<int>("Seed")),
    fRNG(fSeed),
    fRotateMonitorXZ(p.get<double>("RotateMonitorXZ")),
    fRotateMonitorYZ(p.get<double>("RotateMonitorYZ")),
    fNP04FrontZ(p.get<double>("NP04FrontZ")),
    fBeamX(p.get<double>("BeamX")),
    fBeamY(p.get<double>("BeamY")),
    fBeamZ(p.get<double>("BeamZ")),
    fUpstreamZ(p.get<double>("UpstreamZ")),
    fDownstreamZ(p.get<double>("DownstreamZ")),
    fGeneratorZ(p.get<double>("GeneratorZ")) {

  fParticleTypes = p.get<std::vector<std::string>>("ParticleTypes");
  produces<std::vector<simb::MCTruth>>();

}

void PDSPDataDrivenBeam::produce(art::Event& e) {

  // Define the truth collection for this event.
  auto mcTruths = std::make_unique<std::vector<simb::MCTruth>>();

  //Create and produce the true event
  simb::MCTruth truth;
  GenerateTrueEvent(truth);

  //Set in the MCTruth vector and into the event
  mcTruths->push_back(truth);
  e.put(std::move(mcTruths));
  
  if (fVerbose) {
  std::cout << fCurrentEvent << " " <<
               fRandPDG[fCurrentEvent] << " " <<
               fRandMomentum[fCurrentEvent] << " " <<
               fRandHUpstream[fCurrentEvent] << " " <<
               fRandVUpstream[fCurrentEvent] << " " <<
               fRandHDownstream[fCurrentEvent] << " " <<
               fRandVDownstream[fCurrentEvent] << std::endl; 
  }

  //Some output to the hist file
  fOutputPDG = fRandPDG[fCurrentEvent];
  fOutputMomentum = fRandMomentum[fCurrentEvent];
  fOutputHUpstream = fRandHUpstream[fCurrentEvent];
  fOutputVUpstream = fRandVUpstream[fCurrentEvent];
  fOutputHDownstream = fRandHDownstream[fCurrentEvent];
  fOutputVDownstream = fRandVDownstream[fCurrentEvent];
  fOutputTree->Fill();

  //Increment position in sampled particle list
  ++fCurrentEvent;
}

void PDSPDataDrivenBeam::beginJob() {
  //Open input file
  fInputFile = new TFile(fFileName.c_str(), "READ");

  //Get all the PDFs  
  //
  //Move to:
  //horiz/vert fibers in last XBPF
  //Px, Py, Pz -- or P, theta, phi (I like this more I think)
  //
  //Where the directions are taken from the reconstructed
  //tracks from the last 2 XBPFs.
  //
  //The momentum is a bit tricky: we need to sample in 
  //reco space, but generate in true space. We'll have to 
  //unfold from the reco momentum in the spectrometer to 
  //the true energy (taken from nominal MC?) at the generation 
  //point (what is called the face in the beam sim files I think)

  //Get all of the PDFs that we'll use to sample from
  //One for each particle type we're interested in
  //Also get the number of that type of particle
  for (size_t i = 0; i < fParticleTypes.size(); ++i) {
    std::string part_type = fParticleTypes[i];
    fPDFs[part_type] = (THnSparseD*)fInputFile->Get(part_type.c_str());

    std::string n_part_type = "n" + part_type; 
    TVectorD * n_parts = (TVectorD*)fInputFile->Get(n_part_type.c_str());
    fParticleNums[part_type] = (*n_parts)[0];

    fTotalParticles += (*n_parts)[0];
  }

  //Compute fractions of particles to sample from 
  double running_total = 0.;
  for (auto it = fParticleNums.begin(); it != fParticleNums.end(); ++it) {
    running_total += it->second / fTotalParticles;
    fParticleFracs[it->first] = running_total;
  }


  //Get the kinematic ranges from 1D projections of the first PDF (arbitrary)
  //To use for setting the values later on
  THnSparseD * temp_pdf = fPDFs.begin()->second;
  double minima[5];
  double maxima[5];
  for (size_t i = 0; i < 5; ++i) {
    TH1D * temp_hist = (TH1D*)temp_pdf->Projection(i);
    minima[i] = temp_hist->GetXaxis()->GetXmin();
    maxima[i] = temp_hist->GetXaxis()->GetXmax();

    //Free the memory used for the projection
    delete temp_hist;
  }
  
  //Will need to throw some random numbers
  //pdg type
  //reco momentum
  //vert/horiz fibers for both monitors
  //     --- Establish valid ranges?
  //
  //Then throw flat(0,1) 
  //and check against PDF in that 5D bin
  //
  //
  //If rejected, go back to throwing in phase space
  //and check against PDF again. This will fill the 
  //phase space according to the PDF
  //
  //
  //
  //Do this all in the beginJob step, 
  //where we throw all of the numbers until we get 
  //N events, where N is a fcl parameter describing how
  //many events we would like to generate. Save the 
  //kinematics, position, and PDG for these N events and then use
  //in the produce function to create the events. 
  int nSampled = 0;
  double kin_samples[5]; //the point in phase space to check against pdf
  double pdf_check; //the number used for the checking
  double pdg_sample = 0.; //point to sample fractions of particle types
  while (nSampled < fNGenerate) {
    //Get the pdg from the random num   
    pdg_sample = fRNG.Rndm();
    std::string part_type = "";
    for (auto it = fParticleFracs.begin(); it != fParticleFracs.end(); ++it) {
      if (pdg_sample <= it->second) {
        if (fVerbose) {
          std::cout << "Sampled " << it->first << " " <<
                       pdg_sample << std::endl;
        }
        fRandPDG.push_back(fNameToPDG[it->first]);
        part_type = it->first;
        break;
      }
    }
    
    //Now, find a random point in phase space, throw a number,
    //and check against the PDF
    //
    //If successful, store the kinematics and store for later production
    //If unsuccessful, select a new random point in phase, throw & check again
    bool sample_again = true;
    while (sample_again) {
      fRNG.RndmArray(5, &kin_samples[0]);
      pdf_check = fRNG.Rndm();  

      //Need to convert the numbers sampled for the kinematics (0, 1) 
      //to within the sampling range
      double kin_point[5];
      Convert(kin_samples, minima, maxima, kin_point);


      //Find the bin in the THnSparseD. If the bin has a value of 0,
      //then it would not have been allocated to save on memory.
      //The false parameter prevents that bin from being allocated here,
      //to save memory
      long long bin = fPDFs[part_type]->GetBin(&kin_point[0], false);

      //The bin has no chance of being populated, move on
      if (bin == -1) continue;

      //Find how likely we are to populate this bin
      double pdf_value = fPDFs[part_type]->GetBinContent(bin);

      //If successful, save info and move on
      if (pdf_check <= pdf_value) {
        if (fVerbose) {
          std::cout << "bin: " << bin << " PDF val: " << pdf_value <<
                       " Check: " << pdf_check << std::endl;
          std::cout << kin_samples[0] << " " << kin_samples[1] << " " <<
                       kin_samples[2] << " " << kin_samples[3] << " " <<
                       kin_samples[4] << std::endl;
          std::cout << kin_point[0] << " " << kin_point[1] << " " <<
                       kin_point[2] << " " << kin_point[3] << " " <<
                       kin_point[4] << std::endl;
          std::cout << "Will keep" << std::endl;
        }

        fRandMomentum.push_back(kin_point[0]);
        fRandVUpstream.push_back(int(kin_point[1]));
        fRandHUpstream.push_back(int(kin_point[2]));
        fRandVDownstream.push_back(int(kin_point[3]));
        fRandHDownstream.push_back(int(kin_point[4]));

        ++nSampled;
        sample_again = false;
      }
    }
  }

  //For quick output info
  art::ServiceHandle<art::TFileService> tfs;
  fOutputTree = tfs->make<TTree>("tree", ""); 
  fOutputTree->Branch("PDG", &fOutputPDG);
  fOutputTree->Branch("Event", &fCurrentEvent);
  fOutputTree->Branch("Momentum", &fOutputMomentum);
  fOutputTree->Branch("HUpstream", &fOutputHUpstream);
  fOutputTree->Branch("VUpstream", &fOutputVUpstream);
  fOutputTree->Branch("HDownstream", &fOutputHDownstream);
  fOutputTree->Branch("VDownstream", &fOutputVDownstream);


  //Setting up for positioning/momentum projection
  BeamMonitorBasisVectors();

}

//To turn the kinematic sample (range 0 --> 1)
//into the variable we use in the pdf (momentum, fiber numbers)
void PDSPDataDrivenBeam::Convert(
    double input_point[5], double minima[5],
    double maxima[5], double output_point[5]) {

  for (int i = 0; i < 5; ++i) {
    double delta = maxima[i] - minima[i];
    output_point[i] = minima[i] + delta * input_point[i];
  }
  
}

void PDSPDataDrivenBeam::GenerateTrueEvent(simb::MCTruth &mcTruth) {
  TLorentzVector position(0, 0, 0, 0);
  TLorentzVector momentum(0., 0., 0., 0.);
  SetPositionAndMomentum(position, momentum);
  
  if (fVerbose) {
    std::cout << "Position: ";
    position.Print();
    std::cout << "Momentum: ";
    momentum.Print();
  }

  mcTruth.SetOrigin(simb::kSingleParticle);

  std::string process = "primary";
  int trackID = 0; //Change?

  simb::MCParticle newParticle(trackID, fRandPDG[fCurrentEvent], process);
  newParticle.AddTrajectoryPoint(position, momentum);
  
  // Add the MCParticle to the MCTruth for the event.
  mcTruth.Add(newParticle);
}

void PDSPDataDrivenBeam::SetPositionAndMomentum(TLorentzVector & position, TLorentzVector & momentum) {
  //First, get the vector between the two sets of tracking monitors
  //
  //Shift to center of fibers. Middle is 96
  //Fiber 0 is most-positive: 95.5
  //Fiber 191 is most-negative: -95.5
  double upstreamX = 96. - int(fRandHUpstream[fCurrentEvent]) - .5;
  double upstreamY = 96. - int(fRandVUpstream[fCurrentEvent]) - .5;
  TVector3 upstream_point = ConvertProfCoordinates(upstreamX, upstreamY, 0., 
                                                   fUpstreamZ);

  double downstreamX = 96. - int(fRandHDownstream[fCurrentEvent]) - .5;
  double downstreamY = 96. - int(fRandVDownstream[fCurrentEvent]) - .5;
  TVector3 downstream_point = ConvertProfCoordinates(downstreamX, downstreamY, 0., 
                                                     fDownstreamZ);

  TVector3 dR = (downstream_point - upstream_point).Unit();

  if (fVerbose) {
  std::cout << "Upstream: " << upstream_point.X() << " " << 
               upstream_point.Y() << " " << upstream_point.Z() << std::endl;
  std::cout << "Downstream: " << downstream_point.X() << " " <<
               downstream_point.Y() << " " << downstream_point.Z() << std::endl;
  std::cout << "dR: " << dR.X() << " " << dR.Y() << " " << dR.Z() << std::endl;
  }

  //Project to generator point
  double deltaZ = (fGeneratorZ - downstream_point.Z());
  double deltaX = (dR.X() / dR.Z()) * deltaZ; 
  double deltaY = (dR.Y() / dR.Z()) * deltaZ; 

  TVector3 generator_point = downstream_point +
                             TVector3(deltaX, deltaY, deltaZ);
  //Set the position 4-vector
  //Time = 0 for now?
  position = TLorentzVector(generator_point, 0.);

  //Prints out the projected position at the face of the TPC
  if (fVerbose) {
    deltaZ = (-1.*fGeneratorZ);
    deltaX = (dR.X() / dR.Z()) * deltaZ; 
    deltaY = (dR.Y() / dR.Z()) * deltaZ; 

    TVector3 last_point = generator_point +
                          TVector3(deltaX, deltaY, deltaZ);

    std::cout << last_point.X() << " " << last_point.Y() << " " <<
                 last_point.Z() << std::endl;
  }
  
  //Scales the trajectory between monitors by the momentum value
  //Later: need to unfold from the sampled, reconstructed momentum
  //into 'true' momentum
  TVector3 mom_vec = fRandMomentum[fCurrentEvent]*dR;
  
  //Get the PDG and set the mass & energy accordingly
  int pdg = fRandPDG[fCurrentEvent];
  const TDatabasePDG * dbPDG = TDatabasePDG::Instance();
  const TParticlePDG * def = dbPDG->GetParticle(pdg);
  double mass = def->Mass();
  double energy = sqrt(mass*mass + mom_vec.Mag2());
  
  //Set the momentum 4-vector
  momentum = TLorentzVector(mom_vec, energy);
}

TVector3 PDSPDataDrivenBeam::ConvertProfCoordinates(
    double x, double y, double z, double zOffset) {

  double off = fNP04FrontZ - zOffset;

  TVector3 old(x,y,z);

  double newX = x*fBMBasisX.X() + y*fBMBasisY.X() + off*fabs(fBMBasisZ.X());
  double newY = x*fBMBasisX.Y() + y*fBMBasisY.Y() + off*fabs(fBMBasisZ.Y());
  double newZ = x*fBMBasisX.Z() + y*fBMBasisY.Z() - off*fabs(fBMBasisZ.Z());

  newX += fBeamX*10.;
  newY += fBeamY*10.;
  newZ += fBeamZ*10.;

  TVector3 result(newX/10., newY/10., newZ/10.);
  return result;
}

void PDSPDataDrivenBeam::BeamMonitorBasisVectors(){
  RotateMonitorVector(fBMBasisX);
  RotateMonitorVector(fBMBasisY);
  RotateMonitorVector(fBMBasisZ);
}

void PDSPDataDrivenBeam::RotateMonitorVector(TVector3 &vec){
  vec.RotateX( fRotateMonitorYZ * TMath::Pi()/180. );
  vec.RotateY( fRotateMonitorXZ * TMath::Pi()/180. );
}

DEFINE_ART_MODULE(PDSPDataDrivenBeam)
