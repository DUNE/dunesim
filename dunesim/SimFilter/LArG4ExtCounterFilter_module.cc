#include <iostream>
#include <utility>
#include <set>


#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"

namespace filt{

  class LArG4ExtCounterFilter : public art::EDFilter {
    public:
      explicit LArG4ExtCounterFilter(fhicl::ParameterSet const & pset);
      virtual ~LArG4ExtCounterFilter() {};
      virtual bool filter(art::Event& e);
      void reconfigure(fhicl::ParameterSet const& pset);
      void beginJob() ;

    private:

      //Small class to hold the opposing set of counters
      struct CounterSetPair{
        std::vector<unsigned int> setA; //First set
        std::vector<unsigned int> setB; //Second set
        bool IsRequested; //Check if we care about a particular pairwise set of counters (e.g. does the user care about the EW counters etc)
      };
      std::vector<CounterSetPair> fCounterSetPairs;

      //Adjustable parameters
      bool fUseEWCounterPair; //Use the EW counter pair
      bool fUseNupSdownCounterPair; //Use North (up) South (down) counter pair
      bool fUseNdownSupCounterPair; //Use the North (down) South (up) counter pair
      std::vector<int> fInterestingPDGs; //A vector of particle PDGs which want to be filtered on
      double fParticleMinEnergy;  //The minimum energy of a particle to be filtered
      double fParticleMaxEnergy;  //The maximum energy of a particle to be filtered
  
      bool IsInterestingParticle(const art::Ptr<simb::MCParticle> particle);  //Define whether a particular particle is initially worth saving e.g. is it a muon, pion etc
      bool UsesCounterSetPair(const CounterSetPair &CSP, const std::set<unsigned int> &usedCounterIDs); //Check if a list of counter IDs matches with those in a particular set pair
      bool UsesCounterSet(const std::vector<unsigned int> &counterIDs, const std::set<unsigned int> &usedCounterIDs); //Check if any of a list of counters are counter in a particular counter set

  };

  LArG4ExtCounterFilter::LArG4ExtCounterFilter::LArG4ExtCounterFilter(fhicl::ParameterSet const & pset)
  {
    this->reconfigure(pset);
  }

  void LArG4ExtCounterFilter::reconfigure(fhicl::ParameterSet const& pset){
    fUseEWCounterPair = pset.get<bool>("UseEWCounterPair",1);
    std::cout<<"Use EW counter pair: " << fUseEWCounterPair<<std::endl;
    fUseNupSdownCounterPair = pset.get<bool>("UseNupSdownCounterPair",1);
    std::cout<<"Use N (up) S (down) counter pair: " << fUseNupSdownCounterPair << std::endl;
    fUseNdownSupCounterPair = pset.get<bool>("UseNdownSupCounterPair",1);
    std::cout<<"Use N (down) S (up) counter pair: " << fUseNdownSupCounterPair << std::endl;
    fInterestingPDGs = pset.get<std::vector<int> >("InterestingPDGs");
    std::cout<<"NInteresting PDGs: " << fInterestingPDGs.size() << std::endl;
    for (unsigned int i = 0; i < fInterestingPDGs.size(); i++){
      std::cout<<"-- PDG: " << fInterestingPDGs[i] << std::endl; 
    }
    fParticleMinEnergy = pset.get<double>("ParticleMinEnergy",0);
    std::cout<<"Particle min energy: " << fParticleMinEnergy << std::endl;
    fParticleMaxEnergy = pset.get<double>("ParticleMaxEnergy",999999999);
    std::cout<<"Particle max energy: " << fParticleMaxEnergy << std::endl;
  }

  bool LArG4ExtCounterFilter::filter(art::Event & e){

    art::ServiceHandle<geo::Geometry> geom;


    //Get the vector of particles
    art::Handle<std::vector<simb::MCParticle> > particles;
    e.getByLabel("largeant",particles);
    //Loop through the particles
    for (unsigned int part_i = 0; part_i < particles->size(); part_i++){
      //Get the particle
      const art::Ptr<simb::MCParticle> particle(particles,part_i);
      //Check if the particle is initially worth considering
      if (IsInterestingParticle(particle)){
        //OK so the particle matches the criteria we need.  Now we need to get the IDs of all external counters it passes through
        //To do this, get the trajectory points and find which counter (if any) said points are contained in
        //Use a set to store the counter IDs the particle passes through
        std::set<unsigned int> usedExtCounterIDs;
        //Now get the trajcectory points
        unsigned int npts = particle->NumberTrajectoryPoints();
        //Loop over those points
        for (unsigned int pt = 0; pt < npts; pt++){
          //Get the position at the point
          TLorentzVector pos4 = particle->Position(pt);
          //The function which checks whether a position is contained in a counter requires a double[3].  So convert the position to that format
          double pos[3];
          pos[0] = pos4.X();
          pos[1] = pos4.Y();
          pos[2] = pos4.Z();
          //If the position is not contained in a counter, the function throws an exception.  We don't want to end the process when this happens
          try{
            //std::cout<<"AuxDetID: " << geom->FindAuxDetAtPosition(pos) << "  pdg: " << particle->PdgCode() << std::endl;
            unsigned int counterID = geom->FindAuxDetAtPosition(pos);
            usedExtCounterIDs.insert(counterID);
          }
          catch(...){};
        }
        std::cout<<"NCounters (any) passed through: " << usedExtCounterIDs.size() << std::endl;
        //We now have a list (set) of the counter IDs which the particles passed through.  Now compare this list with the information stored in the CounterSetPairs
        for (unsigned int csp_i = 0; csp_i < fCounterSetPairs.size(); csp_i++){
          //Only bother comparing the information if a particlar counter set pair has been requested
          if (!fCounterSetPairs[csp_i].IsRequested) continue;
          if (UsesCounterSetPair(fCounterSetPairs[csp_i],usedExtCounterIDs)){
            std::cout<<"SELECTED EVENT"<<std::endl;
            return true;
          }
        }
      }
    }

    //Assume that the event is not worth saving
    return false;
  }

  void LArG4ExtCounterFilter::beginJob() {

    //Need to get the counter information.  By doing this at the start of the job, rather than per event, the code assumes the geomtry is not going to change between events
    art::ServiceHandle<geo::Geometry> geom;
    //Create the pairwise counter sets
    CounterSetPair EWCounterSetPair;
    CounterSetPair NupSdownCounterSetPair;
    CounterSetPair NdownSupCounterSetPair;
    //A stupid way of storing the IDs of the counters, this REALLY needs changing
    //The code loops through all of the counters in the geomtry, and if the number matches a particular counter number e.g. one of the east counters, store it in the correct, pairwise set
    for (unsigned int i = 0; i < geom->NAuxDets(); i++){
      //The WE counter pairs
      if (i >=6 && i <= 15) EWCounterSetPair.setA.push_back(i);
      else if (i >= 28 && i <=37) EWCounterSetPair.setB.push_back(i);
      //The N (up) S (down) counter pairs
      else if (i >= 22 && i <= 27) NupSdownCounterSetPair.setA.push_back(i);
      else if (i <= 5) NupSdownCounterSetPair.setB.push_back(i);
      //The N (down) S (up) counter pairs
      else if (i >= 16 && i <= 21) NdownSupCounterSetPair.setA.push_back(i);
      else if (i >= 38 && i <= 43) NdownSupCounterSetPair.setB.push_back(i);
    }

    //Enable/disable the counter set pairs
    EWCounterSetPair.IsRequested = fUseEWCounterPair;
    NupSdownCounterSetPair.IsRequested = fUseNupSdownCounterPair;
    NdownSupCounterSetPair.IsRequested = fUseNdownSupCounterPair;

    //Store them onto the vector of counter set pairs
    fCounterSetPairs.push_back(EWCounterSetPair);
    fCounterSetPairs.push_back(NupSdownCounterSetPair);
    fCounterSetPairs.push_back(NdownSupCounterSetPair);


    for (unsigned int i = 0; i < fCounterSetPairs.size(); i++){
      std::cout<<"Counter set pair: "<<i<<std::endl;
      std::cout<<"--setA size: " << fCounterSetPairs[i].setA.size() << std::endl;
      std::cout<<"--setB size: " << fCounterSetPairs[i].setB.size() << std::endl;

    }
  }

  bool LArG4ExtCounterFilter::IsInterestingParticle(const art::Ptr<simb::MCParticle> particle){
    //Loop over the list of requested PDGs.  See if that matches the particle under consideration
    for (unsigned int i = 0; i < fInterestingPDGs.size(); i++){
      //Check if the particle under consideration has a requested PDG
      if (particle->PdgCode() == fInterestingPDGs[i]){
        //Got a requested PDG,  now check that the energy matches the requested range
        TLorentzVector mom4 = particle->Momentum();
        if (mom4.T() > fParticleMinEnergy && mom4.T() < fParticleMaxEnergy){
          //std::cout<<"FOUND INTERESTING PARTICLE"<<std::endl;
          return true;
        }
      }
    }
    return false;
  }

  bool LArG4ExtCounterFilter::UsesCounterSetPair(const CounterSetPair &CSP, const std::set<unsigned int> &usedCounterIDs){

    bool usesSetA = UsesCounterSet(CSP.setA,usedCounterIDs);
    bool usesSetB = UsesCounterSet(CSP.setB,usedCounterIDs);
    if (usesSetA && usesSetB){
      std::cout<<"USES BOTH SETS OF COUNTERS"<<std::endl;
      return true;
    }

    return false;
  }

  bool LArG4ExtCounterFilter::UsesCounterSet(const std::vector<unsigned int> &counterIDs, const std::set<unsigned int> &userCounterIDs){
    //Iterate over the used counter IDs and compare them with every counter ID in the vector.  If a match is found, this means the particle of interest passed through one of the counters in the vector
    for (std::set<unsigned int>::iterator setIt = userCounterIDs.begin(); setIt != userCounterIDs.end(); setIt++){
      unsigned int usedCounterID = (*setIt);
      //Now loop through the IDs in the vector
      for (unsigned int i = 0; i < counterIDs.size(); i++){
        //Compare the elements.  If there is a match, return true
        if (usedCounterID == counterIDs[i]){
          std::cout<<"Counter ID match"<<std::endl;
          return true;
        }
      }
    }

    //None of the used counters matched the counters in the set
    return false;
  }

  DEFINE_ART_MODULE(LArG4ExtCounterFilter)
}
