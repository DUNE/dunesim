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
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "dunecore/DuneObj/ProtoDUNEBeamSpill.h"

#include <memory>
#include <utility>
#include <map>

#include "TFile.h"
#include "THnSparse.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
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
  std::string fInputFileName;
  TFile* fResolutionFile = 0x0;
  std::string fResolutionFileName;
  std::vector<std::string> fParticleTypes;
  int fNGenerate;
  int fSeed;
  TRandom3 fRNG;
  bool fVerbose;
  double fRotateMonitorXZ, fRotateMonitorYZ;
  double fNP04FrontZ, fBeamX, fBeamY, fBeamZ;
  double fUpstreamZ, fDownstreamZ, fGeneratorZ;
  int fUnsmearType;
  std::vector<double> fMinima, fMaxima;
  int fNominalMomentum;
  bool fDoWeights;


  std::map<std::string, THnSparseD *> fPDFs;
  std::map<std::string, TH1D *> fResolutionHists;
  std::map<std::string, TH2D *> fResolutionHists2D;
  std::map<std::string, TH2D *> fResolutionHists2DPlus;
  std::map<std::string, TH2D *> fResolutionHists2DMinus;
  std::map<std::string, TF1 *> fResolutions;
  std::map<std::string, double> fParticleNums;
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
  double fUnfoldedMomentum;
  double fOutputHUpstream, fOutputVUpstream;
  double fOutputHDownstream, fOutputVDownstream;
  double fPlusSigmaWeight, fMinusSigmaWeight;

  std::map<std::string, int> fNameToPDG = {
    {"Protons", 2212},
    {"Pions", 211},
    {"Electrons", -11},
    {"Muons", -13},
    {"Kaons", 321}
  };

  std::map<int, std::string> fPDGToName = {
    {2212, "Protons"},
    {211, "Pions"},
    {-11, "Electrons"},
    {-13, "Muons"},
    {321, "Kaons"}
  };

  TVector3 fBMBasisX = TVector3(1.,0.,0.);
  TVector3 fBMBasisY = TVector3(0.,1.,0.);
  TVector3 fBMBasisZ = TVector3(0.,0.,1.);


  void Convert(double input_point[5], std::vector<double> minima,
               std::vector<double> maxima, double output_point[5]);
  //void Convert(double input_point[5], double minima[5],
  //             double maxima[5], double output_point[5]);

  void GenerateTrueEvent(simb::MCTruth &mcTruth);
  void GenerateBeamEvent(beam::ProtoDUNEBeamEvent & beamEvent);
  void MakeTracks(beam::ProtoDUNEBeamEvent & beamEvent);
  beam::FBM MakeFiberMonitor(double pos);
  void SetPositionAndMomentum(TLorentzVector & position,
                              TLorentzVector & momentum);
  TVector3 ConvertProfCoordinates(double x, double y, double z, double zOffset);
  void BeamMonitorBasisVectors();
  void RotateMonitorVector(TVector3 & vec);
  void FitResolutions();
  double UnsmearMomentum1D();
  double UnsmearMomentum2D();
  void GetSystWeights();

  void Setup1GeV(); //Shared by 2GeV
  void Setup3GeV();
  void Setup6GeV(); //Shared by 7GeV
  
  void Scale2DRes();
};


PDSPDataDrivenBeam::PDSPDataDrivenBeam(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fInputFileName(p.get<std::string>("InputFileName")),
    fResolutionFileName(p.get<std::string>("ResolutionFileName")),
    fNGenerate(p.get<int>("NGenerate")),
    fSeed(p.get<int>("Seed")),
    fRNG(fSeed),
    fVerbose(p.get<bool>("Verbose")),
    fRotateMonitorXZ(p.get<double>("RotateMonitorXZ")),
    fRotateMonitorYZ(p.get<double>("RotateMonitorYZ")),
    fNP04FrontZ(p.get<double>("NP04FrontZ")),
    fBeamX(p.get<double>("BeamX")),
    fBeamY(p.get<double>("BeamY")),
    fBeamZ(p.get<double>("BeamZ")),
    fUpstreamZ(p.get<double>("UpstreamZ")),
    fDownstreamZ(p.get<double>("DownstreamZ")),
    fGeneratorZ(p.get<double>("GeneratorZ")),
    fUnsmearType(p.get<int>("UnsmearType")),
    fNominalMomentum(p.get<int>("NominalMomentum")),
    fDoWeights(p.get<bool>("DoWeights")) {

  std::vector<std::pair<std::string, double>> temp_vec =
      p.get<std::vector<std::pair<std::string, double>>>("ParticleTypes");

  for (auto it = temp_vec.begin(); it != temp_vec.end(); ++it) {
    fParticleTypes.push_back(it->first);
    fParticleNums[it->first] = it->second;
    fTotalParticles += it->second;
    if (fVerbose)
      std::cout << it->first << " " << it->second << " " << fTotalParticles <<
                   std::endl;
  }

  fMinima = p.get<std::vector<double>>("Minima");
  fMaxima = p.get<std::vector<double>>("Maxima");

  produces<std::vector<simb::MCTruth>>();
  produces<std::vector<beam::ProtoDUNEBeamEvent>>();

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

  //Define beam event for the event
  std::unique_ptr<std::vector<beam::ProtoDUNEBeamEvent>>
      beamData(new std::vector<beam::ProtoDUNEBeamEvent>);
  beam::ProtoDUNEBeamEvent beamEvent;

  GenerateBeamEvent(beamEvent);
  beamData->push_back(beamEvent);
  e.put(std::move(beamData));

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
  fInputFile = new TFile(fInputFileName.c_str(), "READ");
  fResolutionFile = new TFile(fResolutionFileName.c_str(), "READ");

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
  /*
  for (size_t i = 0; i < fParticleTypes.size(); ++i) {
    std::string part_type = fParticleTypes[i];
    fPDFs[part_type] = (THnSparseD*)fInputFile->Get(part_type.c_str());

    std::string n_part_type = "n" + part_type;
    TVectorD * n_parts = (TVectorD*)fInputFile->Get(n_part_type.c_str());
    fParticleNums[part_type] = (*n_parts)[0];

    fTotalParticles += (*n_parts)[0];


    //Also get the resolutions
    std::string res_name = "h" + part_type + "Res";
    fResolutionHists[part_type] = (TH1D*)fResolutionFile->Get(res_name.c_str());

    res_name += "2D";
    fResolutionHists2D[part_type] = (TH2D*)fResolutionFile->Get(res_name.c_str());
    //Scale according to the reco bin population
    TH2D * this_hist = fResolutionHists2D[part_type];
    for (int i = 1; i <= this_hist->GetNbinsX(); ++i) {
      double integral = this_hist->TH1::Integral(i, i);
      std::cout << "Integral: " << integral << std::endl;
      double total = 0.;
      for (int j = 1; j <= this_hist->GetNbinsY(); ++j) {
        double content = this_hist->GetBinContent(i, j);
        std::cout << "Content: " << content << std::endl;
        this_hist->SetBinContent(i, j,
            this_hist->GetBinContent(i, j) / integral);
        total += this_hist->GetBinContent(i, j);
      }
      std::cout << "Bin " << i << " total " << total << std::endl;
    }
  }
  */

  switch (fNominalMomentum) {
    case 1:
      Setup1GeV();
      break;
    case 2:
      Setup1GeV();
      break;
    case 3:
      Setup3GeV();
      break;
    case 6:
      Setup6GeV();
      break;
    case 7:
      Setup6GeV();
      break;
    default:
      Setup1GeV();
      break;
  }

  Scale2DRes();

  //Compute fractions of particles to sample from
  double running_total = 0.;
  for (auto it = fParticleNums.begin(); it != fParticleNums.end(); ++it) {
    running_total += it->second / fTotalParticles;
    fParticleFracs[it->first] = running_total;
    if (fVerbose) {
      std::cout << it->second << " " << running_total << std::endl;
      std::cout << fParticleFracs[it->first] << std::endl;
    }
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
      //Convert(kin_samples, minima, maxima, kin_point);
      Convert(kin_samples, fMinima, fMaxima, kin_point);


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
        fRandVUpstream.push_back(kin_point[1]);
        fRandHUpstream.push_back(kin_point[2]);
        fRandVDownstream.push_back(kin_point[3]);
        fRandHDownstream.push_back(kin_point[4]);

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
  fOutputTree->Branch("UnfoldedMomentum", &fUnfoldedMomentum);
  fOutputTree->Branch("HUpstream", &fOutputHUpstream);
  fOutputTree->Branch("VUpstream", &fOutputVUpstream);
  fOutputTree->Branch("HDownstream", &fOutputHDownstream);
  fOutputTree->Branch("VDownstream", &fOutputVDownstream);
  fOutputTree->Branch("PlusSigmaWeight", &fPlusSigmaWeight);
  fOutputTree->Branch("MinusSigmaWeight", &fMinusSigmaWeight);


  //Setting up for positioning/momentum projection
  BeamMonitorBasisVectors();

  //Fit the resolution hists
  FitResolutions();
}

//To turn the kinematic sample (range 0 --> 1)
//into the variable we use in the pdf (momentum, fiber numbers)
void PDSPDataDrivenBeam::Convert(
    double input_point[5], std::vector<double> minima,
    std::vector<double> maxima, double output_point[5]) {
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
  /*
  double upstreamX = 96. - int(fRandHUpstream[fCurrentEvent]) - .5;
  double upstreamY = 96. - int(fRandVUpstream[fCurrentEvent]) - .5;

  double downstreamX = 96. - int(fRandHDownstream[fCurrentEvent]) - .5;
  double downstreamY = 96. - int(fRandVDownstream[fCurrentEvent]) - .5;
  */

  double upstreamX = 96. - fRandHUpstream[fCurrentEvent];
  double upstreamY = 96. - fRandVUpstream[fCurrentEvent];
  double downstreamX = 96. - fRandHDownstream[fCurrentEvent];
  double downstreamY = 96. - fRandVDownstream[fCurrentEvent];

  TVector3 upstream_point = ConvertProfCoordinates(upstreamX, upstreamY, 0.,
                                                   fUpstreamZ);
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
  //Attempting to go from reco to true
  int pdg = fRandPDG[fCurrentEvent];
  switch (fUnsmearType) {
    case 1:
      fUnfoldedMomentum = UnsmearMomentum1D();
      break;
    case 2:
      fUnfoldedMomentum = UnsmearMomentum2D();
      if (fDoWeights) GetSystWeights();
      break;
    default:
      //Just do 1D
      fUnfoldedMomentum = UnsmearMomentum1D();
      break;
  }

  //TVector3 mom_vec = fRandMomentum[fCurrentEvent]*dR;
  //TVector3 mom_vec = unfolded_momentum*dR;
  TVector3 mom_vec = fUnfoldedMomentum*dR;

  //Get the PDG and set the mass & energy accordingly
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

  //Transform the positions relative to the center of the beam monitor to that
  //relative to the end of the beam line
  double newX = x*fBMBasisX.X() + y*fBMBasisY.X() + off*fabs(fBMBasisZ.X());
  double newY = x*fBMBasisX.Y() + y*fBMBasisY.Y() + off*fabs(fBMBasisZ.Y());
  double newZ = x*fBMBasisX.Z() + y*fBMBasisY.Z() - off*fabs(fBMBasisZ.Z());

  //Shift the position relative to the beam entrance on the cryostat
  newX += fBeamX*10.;
  newY += fBeamY*10.;
  newZ += fBeamZ*10.;

  TVector3 result(newX/10., newY/10., newZ/10.);
  return result;
}

void PDSPDataDrivenBeam::BeamMonitorBasisVectors() {
  RotateMonitorVector(fBMBasisX);
  RotateMonitorVector(fBMBasisY);
  RotateMonitorVector(fBMBasisZ);
}

void PDSPDataDrivenBeam::RotateMonitorVector(TVector3 &vec) {
  vec.RotateX(fRotateMonitorYZ * TMath::Pi()/180.);
  vec.RotateY(fRotateMonitorXZ * TMath::Pi()/180.);
}

void PDSPDataDrivenBeam::FitResolutions() {
  for (auto it = fResolutionHists.begin(); it != fResolutionHists.end(); ++it) {
    TH1D * hist = it->second;
    hist->Fit("gaus", "Q");
    fResolutions[it->first] = (TF1 *)hist->GetFunction("gaus");
  }
}

void PDSPDataDrivenBeam::GenerateBeamEvent(
    beam::ProtoDUNEBeamEvent & beamEvent) {
  beamEvent.SetTOFs(std::vector<double>{0.});
  beamEvent.SetTOFChans(std::vector<int>{0});
  beamEvent.SetUpstreamTriggers(std::vector<size_t>{0});
  beamEvent.SetDownstreamTriggers(std::vector<size_t>{0});
  beamEvent.SetCalibrations(0., 0., 0., 0.);
  beamEvent.DecodeTOF();

  beamEvent.SetMagnetCurrent(0.);
  beamEvent.SetTimingTrigger(12);

  beam::CKov dummy;
  dummy.trigger = 0;
  dummy.pressure = 0.;
  dummy.timeStamp = 0.;
  beamEvent.SetCKov0(dummy);
  beamEvent.SetCKov1(dummy);

  beamEvent.SetActiveTrigger(0);
  beamEvent.SetT0(std::make_pair(0.,0.));

  //Dummy positions for these
  beamEvent.SetFBMTrigger("XBPF022697", MakeFiberMonitor(.5));
  beamEvent.SetFBMTrigger("XBPF022698", MakeFiberMonitor(.5));
  beamEvent.SetFBMTrigger("XBPF022701", MakeFiberMonitor(.5));
  beamEvent.SetFBMTrigger("XBPF022702", MakeFiberMonitor(.5));

  double upstream_x = fRandHUpstream[fCurrentEvent];
  double upstream_y = fRandVUpstream[fCurrentEvent];
  double downstream_x = fRandHDownstream[fCurrentEvent];
  double downstream_y = fRandVDownstream[fCurrentEvent];

  beamEvent.SetFBMTrigger("XBPF022707", MakeFiberMonitor(96. - upstream_x));//X
  beamEvent.SetFBMTrigger("XBPF022708", MakeFiberMonitor(96. - upstream_y));//Y

  beamEvent.SetFBMTrigger("XBPF022716", MakeFiberMonitor(96. - downstream_x));//X
  beamEvent.SetFBMTrigger("XBPF022717", MakeFiberMonitor(96. - downstream_y));//Y

  MakeTracks(beamEvent);

  beamEvent.AddRecoBeamMomentum(fRandMomentum[fCurrentEvent]);
}

beam::FBM PDSPDataDrivenBeam::MakeFiberMonitor(double pos) {
  beam::FBM theFBM;

  //I should probably just make this into
  //a constructor for the FBM...
  theFBM.ID = -1;
  theFBM.glitch_mask = {};
  std::uninitialized_fill(std::begin(theFBM.fiberData),
                          std::end(theFBM.fiberData), 0.);
  std::uninitialized_fill(std::begin(theFBM.timeData),
                          std::end(theFBM.timeData), 0.);
  theFBM.timeStamp = 0.;

  short f = 96 -  short(floor(pos)) - 1;
  theFBM.fibers[f] = 1;
  theFBM.active.push_back(f);
  theFBM.decoded = true;

  return theFBM;
}

void PDSPDataDrivenBeam::MakeTracks(beam::ProtoDUNEBeamEvent & beamEvent) {

  //We should only have one active fiber at a time
  short fx1 = beamEvent.GetFBM("XBPF022707").active[0];
  short fy1 = beamEvent.GetFBM("XBPF022708").active[0];

  //Convert fiber number to position
  //p = 96. - f - .5
  double x1 = 96. - fx1 - .5;
  double y1 = 96. - fy1 - .5;

  TVector3 pos1 = ConvertProfCoordinates(x1, y1, 0., fUpstreamZ);

  short fx2 = beamEvent.GetFBM("XBPF022716").active[0];
  short fy2 = beamEvent.GetFBM("XBPF022717").active[0];

  double x2 = 96. - fx2 - .5;
  double y2 = 96. - fy2 - .5;

  TVector3 pos2 = ConvertProfCoordinates(x2, y2, 0., fDownstreamZ);

  TVector3 dR = (pos2 - pos1).Unit();

  double deltaZ = (-1.*pos2.Z());
  double deltaX = (dR.X() / dR.Z()) * deltaZ;
  double deltaY = (dR.Y() / dR.Z()) * deltaZ;

  TVector3 tpc_point = pos2 + TVector3(deltaX, deltaY, deltaZ);
  if (fVerbose) {
    std::cout << "TPC point: " << tpc_point.X() << " " << tpc_point.Y() <<
                 " " << tpc_point.Z() << std::endl;
  }

  std::vector< TVector3 > thePoints = {pos1, pos2, tpc_point};
  std::vector< TVector3 > theMomenta = {
    (pos2 - pos1).Unit(),
    (pos2 - pos1).Unit(),
    (pos2 - pos1).Unit()
  };

  recob::Track track(
      recob::TrackTrajectory(recob::tracking::convertCollToPoint(thePoints),
                             recob::tracking::convertCollToVector(theMomenta),
                             recob::Track::Flags_t(thePoints.size()),
                             false),
      0, -1., 0, recob::tracking::SMatrixSym55(),
      recob::tracking::SMatrixSym55(), 1);

  beamEvent.AddBeamTrack(track);
}

double PDSPDataDrivenBeam::UnsmearMomentum1D() {
  int pdg = fRandPDG[fCurrentEvent];
  TF1 * res = fResolutions[fPDGToName[pdg]];
  double mean = res->GetParameter(1);
  double sigma = res->GetParameter(2);
  double t = fRNG.Gaus(mean, sigma); //random number from momentum resolution
  return (fRandMomentum[fCurrentEvent]/(t + 1.));
}

double PDSPDataDrivenBeam::UnsmearMomentum2D() {

  if (fVerbose) {
    std::cout << "Using 2D Unsmear method" << std::endl;
  }

  int pdg = fRandPDG[fCurrentEvent];
  TH2D * res = fResolutionHists2D[fPDGToName[pdg]];

  int xBin = res->GetXaxis()->FindBin(fRandMomentum[fCurrentEvent]);
  if (fVerbose) {
    std::cout << "Momentum & bin: " << fRandMomentum[fCurrentEvent] << " " <<
                 xBin << std::endl;
  }

  double true_min = res->GetYaxis()->GetXmin();
  double true_max = res->GetYaxis()->GetXmax();
  for (int i = 1; i <= res->GetNbinsY(); ++i) {
    if (res->GetBinContent(xBin, i) > 0.) {
      true_min = res->GetYaxis()->GetBinLowEdge(i);
      break;
    }
  }
  for (int i = res->GetNbinsY(); i >= 1; --i) {
    if (res->GetBinContent(xBin, i) > 0.) {
      true_max = res->GetYaxis()->GetBinUpEdge(i);
      break;
    }
  }

  if (fVerbose)
    std::cout << "True min and max: " << true_min << " " << true_max << std::endl;
  
  double unsmeared_momentum = 0.;
  while (true) {
    double t = fRNG.Uniform(true_min, true_max);
    double pdf_check = fRNG.Rndm(); //random number to check against PDF 
    
    int yBin = res->GetYaxis()->FindBin(t);
    double pdf_value = res->GetBinContent(xBin, yBin);
    if (fVerbose) {
      std::cout << "True mom & bin: " << t << " " << yBin << std::endl;
      std::cout << "Check & val: " << pdf_check << " " << pdf_value <<
                   std::endl;
    }
    if (pdf_check < pdf_value) {
      unsmeared_momentum = t;
      if (fVerbose) std::cout << "Setting momentum to " << t << std::endl;
      break;
    }
  }

  return unsmeared_momentum;
}

void PDSPDataDrivenBeam::Scale2DRes() {
  for (auto it = fResolutionHists2D.begin();
       it != fResolutionHists2D.end(); ++it) {
    TH2D * this_hist = it->second;
    for (int i = 1; i <= this_hist->GetNbinsX(); ++i) {
      double integral = this_hist->TH1::Integral(i, i);
      // double total = 0.; // unused
      for (int j = 1; j <= this_hist->GetNbinsY(); ++j) {
        this_hist->SetBinContent(i, j,
            this_hist->GetBinContent(i, j) / integral);
        // total += this_hist->GetBinContent(i, j); // unused
      }
    }
  }

  if (fDoWeights) {
    for (auto it = fResolutionHists2DPlus.begin();
         it != fResolutionHists2DPlus.end(); ++it) {
      TH2D * this_hist = it->second;
      for (int i = 1; i <= this_hist->GetNbinsX(); ++i) {
        double integral = this_hist->TH1::Integral(i, i);
        // double total = 0.; // unused
        for (int j = 1; j <= this_hist->GetNbinsY(); ++j) {
          this_hist->SetBinContent(i, j,
              this_hist->GetBinContent(i, j) / integral);
          // total += this_hist->GetBinContent(i, j); // unused
        }
      }

      this_hist->Divide(fResolutionHists2D[it->first]);
    }

    for (auto it = fResolutionHists2DMinus.begin();
         it != fResolutionHists2DMinus.end(); ++it) {
      TH2D * this_hist = it->second;
      for (int i = 1; i <= this_hist->GetNbinsX(); ++i) {
        double integral = this_hist->TH1::Integral(i, i);
        // double total = 0.; // unused
        for (int j = 1; j <= this_hist->GetNbinsY(); ++j) {
          this_hist->SetBinContent(i, j,
              this_hist->GetBinContent(i, j) / integral);
          // total += this_hist->GetBinContent(i, j); // unused
        }
      }

      this_hist->Divide(fResolutionHists2D[it->first]);
    }
  }
}

void PDSPDataDrivenBeam::Setup1GeV() {
  for (size_t i = 0; i < fParticleTypes.size(); ++i) {
    std::string part_type = fParticleTypes[i];
    if (part_type == "Muons") {
      fPDFs[part_type] = (THnSparseD*)fInputFile->Get("Pions"); 
    }
    else {
      fPDFs[part_type] = (THnSparseD*)fInputFile->Get(part_type.c_str());
    }

    //Also get the resolutions
    std::string res_name = "";
    if (part_type == "Muons") {
      res_name = "hPionsRes";
    }
    else {
      res_name = "h" + part_type + "Res";
    }

    fResolutionHists[part_type] = (TH1D*)fResolutionFile->Get(res_name.c_str());

    res_name += "2D";
    fResolutionHists2D[part_type] = (TH2D*)fResolutionFile->Get(res_name.c_str());

    if (fDoWeights) {
      std::string plus_name = res_name + "Plus";
      fResolutionHists2DPlus[part_type] = (TH2D*)fResolutionFile->Get(plus_name.c_str());

      std::string minus_name = res_name + "Minus";
      fResolutionHists2DMinus[part_type] = (TH2D*)fResolutionFile->Get(minus_name.c_str());
    }
  }
}

void PDSPDataDrivenBeam::Setup3GeV() {


  for (size_t i = 0; i < fParticleTypes.size(); ++i) {
    std::string part_type = fParticleTypes[i];
    if (part_type == "Muons") {
      fPDFs[part_type] = (THnSparseD*)fInputFile->Get("Pions"); 
    }
    else if (part_type == "Kaons") {
      fPDFs[part_type] = (THnSparseD*)fInputFile->Get("Protons"); 
    }
    else {
      fPDFs[part_type] = (THnSparseD*)fInputFile->Get(part_type.c_str());
    }

    //Also get the resolutions
    std::string res_name = "";
    if (part_type == "Muons") {
      res_name = "hPionsRes";
    }
    /*
    else if (part_type == "Kaons") {
      res_name = "hProtonsRes";
    }*/
    else {
      res_name = "h" + part_type + "Res";
    }
    
    fResolutionHists[part_type] = (TH1D*)fResolutionFile->Get(res_name.c_str());

    res_name += "2D";
    fResolutionHists2D[part_type] = (TH2D*)fResolutionFile->Get(res_name.c_str());

    if (fDoWeights) {
      std::string plus_name = res_name + "Plus";
      fResolutionHists2DPlus[part_type] = (TH2D*)fResolutionFile->Get(plus_name.c_str());

      std::string minus_name = res_name + "Minus";
      fResolutionHists2DMinus[part_type] = (TH2D*)fResolutionFile->Get(minus_name.c_str());
    }
  }
}

void PDSPDataDrivenBeam::Setup6GeV() {
  for (size_t i = 0; i < fParticleTypes.size(); ++i) {
    std::string part_type = fParticleTypes[i];
    if (part_type == "Muons" || part_type == "Electrons") {
      fPDFs[part_type] = (THnSparseD*)fInputFile->Get("Pions"); 
    }
    else {
      fPDFs[part_type] = (THnSparseD*)fInputFile->Get(part_type.c_str());
    }

    //Also get the resolutions
    std::string res_name = "";
    if (part_type == "Muons") {
      res_name = "hPionsRes";
    }
    else {
      res_name = "h" + part_type + "Res";
    }

    fResolutionHists[part_type] = (TH1D*)fResolutionFile->Get(res_name.c_str());

    res_name += "2D";
    fResolutionHists2D[part_type] = (TH2D*)fResolutionFile->Get(res_name.c_str());

    if (fDoWeights) {
      std::string plus_name = res_name + "Plus";
      fResolutionHists2DPlus[part_type] = (TH2D*)fResolutionFile->Get(plus_name.c_str());

      std::string minus_name = res_name + "Minus";
      fResolutionHists2DMinus[part_type] = (TH2D*)fResolutionFile->Get(minus_name.c_str());
    }
  }
}

void PDSPDataDrivenBeam::GetSystWeights() {
  int pdg = fRandPDG[fCurrentEvent];

  TH2D * plus_map = fResolutionHists2DPlus[fPDGToName[pdg]];
  TH2D * minus_map = fResolutionHists2DMinus[fPDGToName[pdg]];
  int xBin = plus_map->GetXaxis()->FindBin(fRandMomentum[fCurrentEvent]);
  int yBin = plus_map->GetYaxis()->FindBin(fUnfoldedMomentum);
  fPlusSigmaWeight = plus_map->GetBinContent(xBin, yBin); 
  fMinusSigmaWeight = minus_map->GetBinContent(xBin, yBin); 
}

DEFINE_ART_MODULE(PDSPDataDrivenBeam)
