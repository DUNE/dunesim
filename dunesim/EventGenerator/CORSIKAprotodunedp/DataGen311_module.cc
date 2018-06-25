////////////////////////////////////////////////////////////////////////
// Class:       DataGen311
// Module Type: producer
// File:        DataGen311_module.cc
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

class DataGen311;

namespace evgendp{

  //----------------------------------------------------------------------------

  class Track{
    //holds track input seetings from file

  public:

      Track();
      ~Track();

      TLorentzVector getPosition();
      TLorentzVector getMomentum();

      int run;
      int subrun;
      int event;
      int trackID;
      int pdg;
      double m;
      double theta;
      double phi;
      double length;
      double startX;
      double startY;
      double startZ;
      double mom;

  }; //end class Track

  Track::Track(){

  }
  Track::~Track(){}

  TLorentzVector Track::getPosition(){
    TLorentzVector position(startX, startY, startZ, 0);
    return position;
  }

  TLorentzVector Track::getMomentum(){

    double momX = mom*cos(theta);
    double momY = mom*sin(theta)*sin(phi);
    double momZ = mom*sin(theta)*sin(phi);

    static TDatabasePDG  pdgt;
    TParticlePDG* pdgp = pdgt.GetParticle(pdg);
    if (pdgp) m = pdgp->Mass();

    double energy = sqrt( mom*mom + m*m );

    TLorentzVector momentum( momX, momY, momZ, energy );
    return momentum;
  }

  //----------------------------------------------------------------------------

  class DataGen311 : public art::EDProducer {

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
        Comment("Track list")
      };
    }; //end struct Config

    using Parameters = art::EDProducer::Table<Config>;

    explicit DataGen311(Parameters const& config);

    DataGen311(DataGen311 const &) = delete;
    DataGen311(DataGen311 &&) = delete;
    DataGen311 & operator = (DataGen311 const &) = delete;
    DataGen311 & operator = (DataGen311 &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    void beginJob() override;

    double fP0;
    double fSigmaP;
    std::string fFName;

    std::map< int, std::vector< Track*> > fEventTrackMap;

  };

}//end namespace

////////////////////////////////////////////////////////////////////////////////

evgendp::DataGen311::DataGen311(Parameters const& config)
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



void evgendp::DataGen311::beginJob(){

  //read the external file and store the information on the map
  std::ifstream file;
  file.open( fFName );

  if( !file.is_open() )
    throw cet::exception("DataGen311") << "Can't open file " << fFName <<"\n";

  int eventBefore=0;
  int eventCounter=0;
  std::string line;
  while(std::getline( file, line )){
    std::stringstream sstream(line);

    Track *track = new Track();
    sstream >> track->run >> track->subrun >> track->event >> track->trackID >> track->theta >> track->phi >> track->startX >> track->startY >> track->startZ;
    track->pdg = 13; //<--Hardcoded!

    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat   flat(engine);
    track->mom = fP0 + fSigmaP*(2.0*flat.fire()-1.0);

    if( track->event != eventBefore ){
      eventBefore = track->event;
      eventCounter++;
    }

    //push back in map
    fEventTrackMap[eventCounter].push_back(track);
  }

  file.close();

  //initialize randomn number seed for

}

////////////////////////////////////////////////////////////////////////////////

void evgendp::DataGen311::produce(art::Event & e){


  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine();
  CLHEP::RandFlat flat(engine);

  std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

  simb::MCTruth truth;
  truth.SetOrigin(simb::kCosmicRay);

  //read the particle map and create the particle
  int trackCtr=0;
  for(Track* track : fEventTrackMap[e.id().event()] ){

    simb::MCParticle particle( trackCtr, track->pdg ,"primary",-200, track->m , 1 );
    particle.AddTrajectoryPoint( track->getPosition() , track->getMomentum() );

    truth.Add(particle);

    trackCtr++;
  }

  truthcol->push_back(truth);
  e.put(std::move(truthcol));

}

DEFINE_ART_MODULE(evgendp::DataGen311)
