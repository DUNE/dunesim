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

  class LArG4ParticleFilter : public art::EDFilter {
    public:
      explicit LArG4ParticleFilter(fhicl::ParameterSet const & pset);
      virtual ~LArG4ParticleFilter() {};
      virtual bool filter(art::Event& e);
      void reconfigure(fhicl::ParameterSet const& pset);
      void beginJob() ;

    private:

      bool IsInterestingParticle(const art::Ptr<simb::MCParticle> particle);
      double CalculateLength(const std::vector<TVector3> &position_segment);

      art::ServiceHandle<geo::Geometry> fGeom;


      std::vector<int> fInterestingPDGs;
      double fParticleMinEnergy;
      double fParticleMaxEnergy;
      bool fStopInTPC;
      double fParticleMinTPCLength;

  };

  LArG4ParticleFilter::LArG4ParticleFilter::LArG4ParticleFilter(fhicl::ParameterSet const & pset)
  {
    this->reconfigure(pset);
  }

  void LArG4ParticleFilter::reconfigure(fhicl::ParameterSet const& pset){
    fInterestingPDGs = pset.get<std::vector<int> >("InterestingPDGs");
    std::cout<<"NInteresting PDGs: " << fInterestingPDGs.size() << std::endl;
    for (unsigned int i = 0; i < fInterestingPDGs.size(); i++){
      std::cout<<"-- PDG: " << fInterestingPDGs[i] << std::endl; 
    }
    fParticleMinEnergy = pset.get<double>("ParticleMinEnergy",-1.);
    std::cout<<"Min particle energy: " << fParticleMinEnergy << std::endl;
    fParticleMaxEnergy = pset.get<double>("ParticleMaxEnergy",-1.);
    std::cout<<"Max particle energy: " << fParticleMaxEnergy << std::endl;
    fStopInTPC = pset.get<bool>("StopInTPC",false);
    std::cout<<"Stop in TPC: " << fStopInTPC << std::endl;
    fParticleMinTPCLength = pset.get<double>("ParticleMinTPCLength",-1.);
    std::cout<<"Min particle TPC length: " << fParticleMinTPCLength << std::endl;

    return;
  }

  bool LArG4ParticleFilter::filter(art::Event & e){

    //art::ServiceHandle<geo::Geometry> geom;


    //Get the vector of particles
    art::Handle<std::vector<simb::MCParticle> > particles;
    e.getByLabel("largeant",particles);
    //Loop through the particles
    for (unsigned int part_i = 0; part_i < particles->size(); part_i++){
      //Get the particle
      const art::Ptr<simb::MCParticle> particle(particles,part_i);
      //Check if the particle is initially worth considering
      if (IsInterestingParticle(particle)){
        std::cout<<"Found interesting particle"<<std::endl;
        return true;
      }
    }

    //Assume that the event is not worth saving
    return false;
  }

  void LArG4ParticleFilter::beginJob() {

    return;
  }

  bool LArG4ParticleFilter::IsInterestingParticle(const art::Ptr<simb::MCParticle> particle){
    //The bulk of the code goes here
    //Each check should probably be its own function, but for now everything can be crunched inline

    /*
    //Create a test particle for the purpose testing the filter checks
    simb::MCParticle *testParticle = new simb::MCParticle(10000,13,"test",-1,105,1);
    TLorentzVector test_position1(0.,30.,30.,30);
    TLorentzVector test_position2(0.,40.,40.,40);
    TLorentzVector test_position3(0.,50.,50.,50);
    TLorentzVector test_position4(0.,60.,60.,60);
    TLorentzVector test_position5(0.,70.,70.,70);
    TLorentzVector test_position6(0.,1000.,70.,70);
    TLorentzVector test_position7(0.,1002.,70.,70);
    TLorentzVector test_position8(0.,1003.,70.,70);
    TLorentzVector test_position9(0.,1005.,70.,70);

    //Particle should have several points in the TPC and then several out, followed by the same sequence again.  Finally the last point is within the TPC
    testParticle->AddTrajectoryPoint(test_position1,test_position1);
    testParticle->AddTrajectoryPoint(test_position2,test_position2);
    testParticle->AddTrajectoryPoint(test_position3,test_position3);
    testParticle->AddTrajectoryPoint(test_position4,test_position4);
    testParticle->AddTrajectoryPoint(test_position5,test_position5);
    testParticle->AddTrajectoryPoint(test_position6,test_position6);
    testParticle->AddTrajectoryPoint(test_position7,test_position7);
    testParticle->AddTrajectoryPoint(test_position8,test_position8);
    testParticle->AddTrajectoryPoint(test_position9,test_position9);
    testParticle->AddTrajectoryPoint(test_position1,test_position1);
    testParticle->AddTrajectoryPoint(test_position2,test_position2);
    testParticle->AddTrajectoryPoint(test_position3,test_position3);
    testParticle->AddTrajectoryPoint(test_position4,test_position4);
    testParticle->AddTrajectoryPoint(test_position5,test_position5);
    testParticle->AddTrajectoryPoint(test_position6,test_position6);
    testParticle->AddTrajectoryPoint(test_position7,test_position7);
    testParticle->AddTrajectoryPoint(test_position8,test_position8);
    testParticle->AddTrajectoryPoint(test_position9,test_position9);
    testParticle->AddTrajectoryPoint(test_position1,test_position1);
    */






    //Check the particle PDG
    bool OK = false;
    if (fInterestingPDGs.size() > 0){
      int pdg = particle->PdgCode();
      //Loop through the PDG vector and see if we have a match
      for (unsigned int i = 0; i < fInterestingPDGs.size(); i++){
        if (pdg == fInterestingPDGs[i]){
          OK = true;
          break;
        }
      }
      if (!OK) return false;
    }

    //Check the minimum particle energy
    if (fParticleMinEnergy > 0 && particle->Momentum(0).T() < fParticleMinEnergy) return false;

    //Check the max energy
    if (fParticleMaxEnergy > 0 && particle->Momentum(0).T() > fParticleMaxEnergy) return false;

    //Check if the particle stops in the TPC
    OK = false;
    if (fStopInTPC){
      //Get final position of particle
      TLorentzVector final_position_4vect = particle->Position(particle->NumberTrajectoryPoints()-1);
      double final_position[3];
      final_position[0] = final_position_4vect.X();
      final_position[1] = final_position_4vect.Y();
      final_position[2] = final_position_4vect.Z();

      //std::cout<<"X: " << final_position[0] << "  Y: " << final_position[1] << "  Z: " << final_position[2] << std::endl;
      geo::TPCID tpcid = fGeom->FindTPCAtPosition(final_position);
      //If the tpcid is NOT valid, then the particle did not stop in the TPC so reject this particle
      if (!tpcid.isValid) return false; 
    }

    //Check the length of a particle in the TPC
    OK = false;
    if (fParticleMinTPCLength > 0){
      //To do this, we need to collect the sequential particle positions which are contained in the TPC into segments.  The reason for doing this is that the particle may enter a TPC, leave it and enter it again and if this isn't taken into account, the length might be grossly overestimated
      //It is easiest to store the positions in a vector of vectors
      std::vector< std::vector<TVector3> > position_segments;
      //We are also going to need an empty vector to store in the above vector
      std::vector<TVector3> position_segment;
      //Loop through the trajectory points
      for (unsigned int i = 0; i < particle->NumberTrajectoryPoints(); i++){
        //Extract the current position of the particle
        double curr_pos[3];
        curr_pos[0] = particle->Position(i).X();
        curr_pos[1] = particle->Position(i).Y();
        curr_pos[2] = particle->Position(i).Z();
        geo::TPCID curr_tpcid = fGeom->FindTPCAtPosition(curr_pos);
        //std::cout<<curr_tpcid.isValid<<"  X: "<<curr_pos[0] << "  Y: " << curr_pos[1] << "  Z: " << curr_pos[2] << std::endl;
        //There are a couple of things to check here.  If the particle is currently in the TPC, then we need to store that particular position.  If it is NOT in the TPC, then its either exited the TPC or has not yet entered.  If it has just exited, then the position_segment should have some positions stored in it, it which case we now need to store this segment.  If it has not yet entered the TPC, then we don't need to do anything
        //If it is currently in the TPC
        if (curr_tpcid.isValid) position_segment.push_back(particle->Position(i).Vect());
        //It has just exited the TPC
        else if (position_segment.size() > 0){
          //Store the segment
          position_segments.push_back(position_segment);
          //Now reset the segment
          position_segment.clear();
        }
        //There is nothing to do because the particle has remained outside of the TPC
      }
      //We need to check once more if the position_segment vector has been filled
      if (position_segment.size() > 0){
        position_segments.push_back(position_segment);
        position_segment.clear();
      }
      //Now lets check the length of each segment
      //Firstly, if we didn't store a segment then the particle fails the check
      if (position_segments.size() == 0) return false;
      //Now loop through the segments and check if they are above threshold
      for (unsigned int i = 0; i < position_segments.size(); i++){
        double segment_length = CalculateLength(position_segments[i]);
        if (segment_length > fParticleMinTPCLength){
          //We found a track segment in the TPC which passes the length threshold so don't flag as bad
          OK = true;
          break;
        }
      }
      if (!OK) return false;
    }



    return true;
  }

  double LArG4ParticleFilter::CalculateLength(const std::vector<TVector3> &position_segment){
    double length = 0;
    //Check how many points we have in the segment.  If it is one or less, there is nothing to calculate so return 0
    if (position_segment.size() <= 1) return length;

    //Now we need to compare every adjacent pair of positions to work out the length of this segment
    for (unsigned int i = 1; i < position_segment.size(); i++){
      length += (position_segment[i] - position_segment[i-1]).Mag();
    }

    return length;
  }

  DEFINE_ART_MODULE(LArG4ParticleFilter)
}
