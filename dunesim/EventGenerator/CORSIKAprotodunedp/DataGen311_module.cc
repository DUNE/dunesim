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
      //double theta;
      //double phi;
      double startDirectionX;
      double startDirectionY;
      double startDirectionZ;
      double length;
      double startX;
      double startY;
      double startZ;
      //double mom;
      double energy;

  }; //end class Track

  Track::Track(){

  }
  Track::~Track(){}

  TLorentzVector Track::getPosition(){
    TLorentzVector position(startX, startY, startZ, 0);
    return position;
  }

  TLorentzVector Track::getMomentum(){
    static TDatabasePDG  pdgt;
    TParticlePDG* pdgp = pdgt.GetParticle(pdg);
    if (pdgp) m = pdgp->Mass();

    double mom = sqrt( energy*energy - m*m );
    double momX = mom*startDirectionX;
    double momY = mom*startDirectionY;
    double momZ = mom*startDirectionZ;

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

      fhicl::Atom<int> EventsToProcess{
        Name("EventsToProcess"),
        Comment("Number of events to process")
      };

      fhicl::Atom<int> StartEvent{
        Name("StartEvent"),
        Comment("First event to process")
      };

      fhicl::Atom<int> PDG{
        Name("PDG"),
        Comment("Particle PDG code")
      };

      fhicl::Atom<bool> GetEnergyFromCORSIKA{
        Name("GetEnergyFromCORSIKA"),
        Comment("true: read particle energy from CORSIKA histograms as a fucntion of azimuth angle")
      };

      fhicl::Atom<double> FixedEnergy{
        Name("FixedEnergy"),
        Comment("Fixed particle energy")
      };

      fhicl::Atom<std::string> TrackFile{
        Name("TrackFile"),
        Comment("Track list")
      };

      fhicl::Atom<std::string> HistFile{
        Name("HistFile"),
        Comment("CORSIKA histogram file")
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

    int fEventsToProcess;
    int fStartEvent;
    int fPDG;
    bool fGetEnergyFromCORSIKA;
    int fFixedEnergy;
    std::string fTrackFile;
    std::string fHistFile;

    std::map< int, std::vector< Track*> > fEventTrackMap;

  };

}//end namespace

////////////////////////////////////////////////////////////////////////////////

evgendp::DataGen311::DataGen311(Parameters const& config)
 :
   fEventsToProcess 		(config().EventsToProcess()),
   fStartEvent      		(config().StartEvent()),
   fPDG             		(config().PDG()),
   fGetEnergyFromCORSIKA	(config().GetEnergyFromCORSIKA()),
   fFixedEnergy			(config().FixedEnergy()),
   fTrackFile       		(config().TrackFile()),
   fHistFile        		(config().HistFile())
{
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this);
    //art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "gen", p, { "Seed", "SeedGenerator" });

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();
}

////////////////////////////////////////////////////////////////////////////////



void evgendp::DataGen311::beginJob(){

  //read the files and store the information on the map
  std::ifstream trackfile;
  trackfile.open( fTrackFile );
  if( !trackfile.is_open() )
    throw cet::exception("DataGen311") << "Can't open file " << fTrackFile <<"\n";

  TFile *chistfile = new TFile(fHistFile.c_str());
  if( chistfile->IsOpen() == kFALSE )
    throw cet::exception("DataGen311") << "Can't open file " << fHistFile <<"\n";

  std::stringstream hetss;
  hetss << "hEnergyTheta" << fPDG;
  TH2D *hEnergyTheta = (TH2D*)chistfile->Get(hetss.str().c_str());

  int eventBefore=0;
  int eventCounter=0;
  int trackCounter=0, skippedTrackCounter=0;
  std::string line;
  std::getline(trackfile, line); // Skip header
  while(std::getline( trackfile, line )){
    std::stringstream sstream(line);

    int eventID = 9999999;
    sstream >> eventID;
    if (eventID < fStartEvent) {
      continue;
    } else if (eventID >= fStartEvent+fEventsToProcess) {
      break;
    }

    Track *track = new Track();
    sstream >> track->run >> track->subrun >> track->event >> track->trackID >> track->startDirectionX >> track->startDirectionY >> track->startDirectionZ >> track->length >> track->startX >> track->startY >> track->startZ;
    track->pdg = fPDG;
    trackCounter++;

    // Start direction components might not be perfectly normalized
    double sdirmag = sqrt(track->startDirectionX*track->startDirectionX + track->startDirectionY*track->startDirectionY + track->startDirectionZ*track->startDirectionZ);
    track->startDirectionX /= sdirmag;
    track->startDirectionY /= sdirmag;
    track->startDirectionZ /= sdirmag;


    if(fGetEnergyFromCORSIKA)
    {
      // Theta from +X to -X
      double theta = acos(track->startDirectionX);

      int biny = hEnergyTheta->GetYaxis()->FindBin(theta);
      TH1D *hEnergyAtTheta = hEnergyTheta->ProjectionX("", biny, biny);

      if( hEnergyAtTheta->GetEntries() < 1 ){
        skippedTrackCounter++;
        delete track;
        continue;
      }

      double energy = hEnergyAtTheta->GetRandom();
      track->energy = energy;
    }
    else //get energy from .fcl parameter
    {
      track->energy = fFixedEnergy;
    }

    if( track->event != eventBefore ){
      eventBefore = track->event;
      eventCounter++;
    }

    //push back in map
    fEventTrackMap[eventCounter].push_back(track);
  }

  trackfile.close();

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
