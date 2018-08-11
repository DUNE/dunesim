////////////////////////////////////////////////////////////////////////
// Class:       NEUTImport
// Module Type: producer
// File:        NEUTImport_module.cc
//
// Get input list of selected reco track, generates primary MCParticles
// matching each of them
////////////////////////////////////////////////////////////////////////

// ROOT includes
#include "TRandom3.h"
#include "TDatabasePDG.h"
#include "TString.h"
#include "TSystem.h" //need BaseName and DirName
#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// nutools includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"

// lar includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "TDatabasePDG.h"

#include "CLHEP/Random/RandFlat.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

class NEUTImport;

namespace evgendp{

  //----------------------------------------------------------------------------

  class NDecayParticle{
    //holds NDecayParticle input parameters from file

  public:

      NDecayParticle();
      ~NDecayParticle();

      TLorentzVector getPosition();
      TLorentzVector getMomentum();

      int event;
      int pdg;
      double mass;
      double startX;
      double startY;
      double startZ;
      double momX;
      double momY;
      double momZ;
      double energy;

  }; //end class NDecayParticle

  NDecayParticle::NDecayParticle(){

  }
  NDecayParticle::~NDecayParticle(){}

  TLorentzVector NDecayParticle::getPosition(){
    TLorentzVector position(startX, startY, startZ, 0);
    return position;
  }

  TLorentzVector NDecayParticle::getMomentum(){
    TLorentzVector momentum( 0.001*momX, 0.001*momY, 0.001*momZ, 0.001*energy );
    return momentum;
  }

  //----------------------------------------------------------------------------

  class NEUTImport : public art::EDProducer {

  public:

    struct Config {
      //holds configuration settings from the fcl

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> P0{
        Name("P0"),
        Comment("central momentum (GeV/c) to generate")
      };

      fhicl::Atom<double> SigmaP{
        Name("SigmaP"),
        Comment("variation in momenta (GeV/c)")
      };

      fhicl::Atom<std::string> FName{
        Name("FileName"),
        Comment("NEUT output")
      };
    }; //end struct Config

    using Parameters = art::EDProducer::Table<Config>;

    explicit NEUTImport(Parameters const& config);

    NEUTImport(NEUTImport const &) = delete;
    NEUTImport(NEUTImport &&) = delete;
    NEUTImport & operator = (NEUTImport const &) = delete;
    NEUTImport & operator = (NEUTImport &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    void beginJob() override;

    double fP0;
    double fSigmaP;
    std::string fFName;

    std::map< int, std::vector< NDecayParticle*> > NDecayEventMap;

  };

}//end namespace

////////////////////////////////////////////////////////////////////////////////

evgendp::NEUTImport::NEUTImport(Parameters const& config)
 :
   fP0           (config().P0()),
   fSigmaP       (config().SigmaP()),
   fFName        (config().FName())
{
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this);
    //art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "gen", p, { "Seed", "SeedGenerator" });

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();
}

////////////////////////////////////////////////////////////////////////////////



void evgendp::NEUTImport::beginJob(){

  const int NMaxParticlesPerEvent = 1000;

  TFile *neutFile = new TFile(fFName.c_str(), "READ");
  TTree *neutTree = (TTree*)neutFile->Get("nRooTracker"); //should always be the same
  int NEvents= (int)neutTree->GetEntries();
  std::cout << "Reading in " << NEvents << " NEUT events. " << std::endl;


  int tStdHepN = 0;
  int tStdHepStatus[NMaxParticlesPerEvent] = {0};
  int tStdHepPdg[NMaxParticlesPerEvent] = {0};
  double tStdHepP4[4*NMaxParticlesPerEvent] = {0.};

  neutTree->SetBranchAddress("StdHepN",&tStdHepN);
  neutTree->SetBranchAddress("StdHepStatus",&tStdHepStatus);
  neutTree->SetBranchAddress("StdHepPdg",&tStdHepPdg);
  neutTree->SetBranchAddress("StdHepP4",&tStdHepP4);


  for(int i=0; i<NEvents; i++) //Event loop
  {
    neutTree->GetEntry(i);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Event " << i << std::endl;
    std::cout << "NParticles: " << tStdHepN << std::endl;

    for(int j=0; j<tStdHepN; j++)
    {
      if( tStdHepPdg[j]!=0 && tStdHepStatus[j] == 1 && std::abs(tStdHepP4[4*j]) + std::abs(tStdHepP4[4*j+1]) + std::abs(tStdHepP4[4*j+2]) > 0 )
      {
    	static TDatabasePDG  pdgt;
    	TParticlePDG* pdgp = pdgt.GetParticle(tStdHepPdg[j]);

    	NDecayParticle *ndecayparticle = new NDecayParticle();
	ndecayparticle->pdg = tStdHepPdg[j];
	ndecayparticle->mass = pdgp->Mass();
	ndecayparticle->startX = 0;
	ndecayparticle->startY = -150;
	ndecayparticle->startZ = 150;
	ndecayparticle->momX = tStdHepP4[4*j];
	ndecayparticle->momY = tStdHepP4[4*j+1];
	ndecayparticle->momZ = tStdHepP4[4*j+2];
	ndecayparticle->energy = tStdHepP4[4*j+3];
    	NDecayEventMap[i].push_back(ndecayparticle);
    	std::cout << std::endl;
	std::cout << "Status " << j << ":\t" << tStdHepStatus[j] << std::endl;
	std::cout << "PDG " << j << ":\t" << ndecayparticle->pdg << std::endl;
	std::cout << "Mas " << j << ":\t" << ndecayparticle->mass << std::endl;
	std::cout << "Start.X " << j << ":\t" << ndecayparticle->startX << std::endl;
	std::cout << "Start.Y " << j << ":\t" <<ndecayparticle->startY << std::endl;
	std::cout << "Start.Z " << j << ":\t" << ndecayparticle->startZ << std::endl;
	std::cout << "Mom.X " << j << ":\t" << ndecayparticle->momX << std::endl;
	std::cout << "Mom.Y " << j << ":\t" <<ndecayparticle->momY << std::endl;
	std::cout << "Mom.Z " << j << ":\t" << ndecayparticle->momZ << std::endl;
	std::cout << "Energy " << j << ":\t" << tStdHepP4[4*j+3] << std::endl;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void evgendp::NEUTImport::produce(art::Event & e){


  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine();
  CLHEP::RandFlat flat(engine);

  std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

  simb::MCTruth truth;
  truth.SetOrigin(simb::kUnknown);

  //read the particle map and create the MCParticle
  int ndecayparticlecounter=0;
  for(NDecayParticle* ndecayparticle : NDecayEventMap[e.id().event()] ){

    simb::MCParticle particle( ndecayparticlecounter, ndecayparticle->pdg ,"primary",-200, ndecayparticle->mass , 1 );
    particle.AddTrajectoryPoint( ndecayparticle->getPosition() , ndecayparticle->getMomentum() );

    truth.Add(particle);

    ndecayparticlecounter++;
  }

  truthcol->push_back(truth);
  e.put(std::move(truthcol));

}

DEFINE_ART_MODULE(evgendp::NEUTImport)
