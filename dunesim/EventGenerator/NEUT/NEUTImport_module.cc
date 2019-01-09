////////////////////////////////////////////////////////////////////////
// Class:       NEUTImport
// Module Type: producer
// File:        NEUTImport_module.cc
//
// author: Christoph Alt
// email: christoph.alt@cern.ch
//
// Reads in events from NeutToRooTracker files and generates art events,
// selecting final state particles only.
// NeutToRooTracker: https://github.com/luketpickering/NeutToRooTracker (by Luke Pickering)
// 
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

      fhicl::Atom<double> LogLevel{
        Name("LogLevel"),
        Comment("LogLevel")
      };

      fhicl::Atom<double> StartEvent{
        Name("StartEvent"),
        Comment("StartEvent")
      };

      fhicl::Atom<double> NumberOfEvents{
        Name("NumberOfEvents"),
        Comment("NumberOfEvents")
      };

      fhicl::Atom<double> StartPositionX{
        Name("StartPositionX"),
        Comment("StartPositionX")
      };

      fhicl::Atom<double> StartPositionY{
        Name("StartPositionY"),
        Comment("StartPositionY")
      };

      fhicl::Atom<double> StartPositionZ{
        Name("StartPositionZ"),
        Comment("StartPositionZ")
      };

      fhicl::Atom<std::string> FileName{
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

    int fLogLevel;
    int fStartEvent;
    int fNumberOfEvents;
    double fStartPositionX;
    double fStartPositionY;
    double fStartPositionZ;
    std::string fFileName;

    std::map< int, std::vector< NDecayParticle*> > NDecayEventMap;

  };

}//end namespace

////////////////////////////////////////////////////////////////////////////////

evgendp::NEUTImport::NEUTImport(Parameters const& config)
 : EDProducer{config},
   fLogLevel		(config().LogLevel()),
   fStartEvent		(config().StartEvent()),
   fNumberOfEvents	(config().NumberOfEvents()),
   fStartPositionX	(config().StartPositionX()),
   fStartPositionY	(config().StartPositionY()),
   fStartPositionZ	(config().StartPositionZ()),
   fFileName		(config().FileName())
{
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this);
    //art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "gen", p, { "Seed", "SeedGenerator" });

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();
}

////////////////////////////////////////////////////////////////////////////////



void evgendp::NEUTImport::beginJob(){

  const int NMaxParticlesPerEvent = 1000;

  TFile *neutFile = new TFile(fFileName.c_str(), "READ");
  TTree *neutTree = (TTree*)neutFile->Get("nRooTracker"); //should always be the same
  int NEvents = (int)neutTree->GetEntries();
  std::cout << "Reading in " << NEvents << " NEUT events. " << std::endl;


  int tStdHepN = 0;
  int tStdHepStatus[NMaxParticlesPerEvent] = {0};
  int tStdHepPdg[NMaxParticlesPerEvent] = {0};
  double tStdHepP4[4*NMaxParticlesPerEvent] = {0.};

  neutTree->SetBranchAddress("StdHepN",&tStdHepN);
  neutTree->SetBranchAddress("StdHepStatus",&tStdHepStatus);
  neutTree->SetBranchAddress("StdHepPdg",&tStdHepPdg);
  neutTree->SetBranchAddress("StdHepP4",&tStdHepP4);


  for(int i=fStartEvent; i<fStartEvent+fNumberOfEvents; i++) //Event loop
  {
    neutTree->GetEntry(i);

    if(fLogLevel == 1)
    {
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "Event:\t" <<  i << std::endl;
      std::cout << "NParticles:\t" << tStdHepN << std::endl;
    }

    for(int j=0; j<tStdHepN; j++)
    {
      //Only take particles with real PDG code (!=0), status 1 and momentum > 0.
      //Also ignore first particle (j=0): this is the decaying nucleon in case of nucleon decay or incoming neutrino in case of atmospheric neutrino background.
      if( tStdHepPdg[j]!=0 && tStdHepStatus[j] == 1 && j>0 && std::abs(tStdHepP4[4*j]) + std::abs(tStdHepP4[4*j+1]) + std::abs(tStdHepP4[4*j+2]) > 0 )
      {
    	static TDatabasePDG  pdgt;
    	TParticlePDG* pdgp = pdgt.GetParticle(tStdHepPdg[j]);

    	NDecayParticle *ndecayparticle = new NDecayParticle();
	ndecayparticle->pdg = tStdHepPdg[j];
	ndecayparticle->mass = pdgp->Mass();
	ndecayparticle->startX = fStartPositionX;
	ndecayparticle->startY = fStartPositionY;
	ndecayparticle->startZ = fStartPositionZ;
	ndecayparticle->momX = tStdHepP4[4*j];
	ndecayparticle->momY = tStdHepP4[4*j+1];
	ndecayparticle->momZ = tStdHepP4[4*j+2];
	ndecayparticle->energy = tStdHepP4[4*j+3];
    	NDecayEventMap[i-fStartEvent].push_back(ndecayparticle);

	if(fLogLevel == 1)
	{
	  double fAbsoluteParticleMomentum = sqrt( pow(tStdHepP4[4*j],2) + pow(tStdHepP4[4*j+1],2) + pow(tStdHepP4[4*j+2],2) );
	  std::cout << std::endl;
	  std::cout << "Event #" << i-fStartEvent << " in LArSoft, event #" << i << " in ROOT file." << std::endl;
	  std::cout << "Status particle " << j << ":\t" << tStdHepStatus[j] << std::endl;
	  std::cout << "PDG particle " << j << ":\t" << ndecayparticle->pdg << std::endl;
	  std::cout << "Mass particle " << j << "\t" << ndecayparticle->mass << std::endl;
	  std::cout << "StartPositionX particle " << j << ":\t" << ndecayparticle->startX << std::endl;
	  std::cout << "StartPositionY particle " << j << ":\t" <<ndecayparticle->startY << std::endl;
	  std::cout << "StartPositionZ particle " << j << ":\t" << ndecayparticle->startZ << std::endl;
	  std::cout << "StartMomentumX particle " << j << ":\t" << ndecayparticle->momX << std::endl;
	  std::cout << "StartMomentumY particle " << j << ":\t" <<ndecayparticle->momY << std::endl;
	  std::cout << "StartMomentumZ particle " << j << ":\t" << ndecayparticle->momZ << std::endl;
	  std::cout << "StartEnergy particle " << j << ":\t" << tStdHepP4[4*j+3] << std::endl;
	  std::cout << "Absolute momentum particle " << j << ":\t" << fAbsoluteParticleMomentum << std::endl;
	}
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void evgendp::NEUTImport::produce(art::Event & e){


  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine(art::ScheduleID::first(),
                                                  moduleDescription().moduleLabel());
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
