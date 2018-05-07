#include <iostream>
#include <utility>
#include <complex>
#include <algorithm>


#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"

namespace filt{

  class GenFilter : public art::EDFilter {
    public:
      explicit GenFilter(fhicl::ParameterSet const & pset);
      virtual ~GenFilter() {};
      virtual bool filter(art::Event& e);
      void reconfigure(fhicl::ParameterSet const& pset);
      void beginJob() ;

    private:

      struct CounterSetPair{
        bool isRequested;
        std::vector<geo::AuxDetGeo const*> setA;
        std::vector<geo::AuxDetGeo const*> setB;
        double normalVec[3];

      };
      std::vector<CounterSetPair> fCounterSetPairs;

      //Adjustable parameters
      bool fUseEWCounterPair; //Use the EW counter pair
      bool fUseNupSdownCounterPair; //Use North (up) South (down) counter pair
      bool fUseNdownSupCounterPair; //Use the North (down) South (up) counter pair
      std::vector<int> fInterestingPDGs; //A vector of particle PDGs which want to be filtered on
      double fParticleMinEnergy;  //The minimum energy of a particle to be filtered
      double fParticleMaxEnergy;  //The maximum energy of a particle to be filtered
      double fCounterSizeScaleFactor; //The scaling factor to increase/decrease thedimensions of the counter
 

      bool IsInterestingParticle(const simb::MCParticle &particle);
      bool ParticleHitsCounterSetPairs(const simb::MCParticle &particle, const CounterSetPair &CSP);
      bool ParticleHitsCounterSet(const simb::MCParticle &particle, const std::vector<geo::AuxDetGeo const*> &counters, const TVector3 &counter_norm);
  
  };


  GenFilter::GenFilter::GenFilter(fhicl::ParameterSet const & pset)
  {
    this->reconfigure(pset);
  }

  void GenFilter::reconfigure(fhicl::ParameterSet const& pset){
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

    fCounterSizeScaleFactor = pset.get<double>("CounterSizeScaleFactor",1.);
    std::cout<<"Counter size scale factor: " << fCounterSizeScaleFactor << std::endl;
  }


  bool GenFilter::filter(art::Event & e){

    std::vector< art::Handle< std::vector<simb::MCTruth> > > mclists;
    e.getManyByType(mclists);
    for (unsigned int i = 0; i < mclists.size() ; i++){
      for (unsigned int j = 0; j < mclists[i]->size(); j++){
        //Should have the truth record for the event now
        const art::Ptr<simb::MCTruth> mc_truth(mclists[i],j);
        for (int part = 0; part < mc_truth->NParticles(); part++){
          const simb::MCParticle particle = mc_truth->GetParticle(part);
          if (!IsInterestingParticle(particle)) continue;
          TVector3 particle_pos = particle.Position().Vect();
          TVector3 particle_dir = particle.Momentum().Vect().Unit();
          for (unsigned CSP_i = 0; CSP_i < fCounterSetPairs.size(); CSP_i++){
            if (!fCounterSetPairs[CSP_i].isRequested) continue;
            if (ParticleHitsCounterSetPairs(particle,fCounterSetPairs[i])){
              std::cout<<"HIT COUNTER SET"<<std::endl;
              return true;
            }
            /*
            TVector3 counter_norm(fCounterSetPairs[i].normalVec[0],fCounterSetPairs[i].normalVec[1],fCounterSetPairs[i].normalVec[2]);
            double counter_pos_array[3];
            fCounterSetPairs[CSP_i].setA.front()->GetCenter(counter_pos_array);
            TVector3 counter_pos(counter_pos_array[0], counter_pos_array[1], counter_pos_array[2]);
            double scale_factor = counter_norm.Dot(counter_pos - particle_pos)/(counter_norm.Dot(particle_dir));
            TVector3 particle_pos_in_plane = particle_pos + scale_factor * particle_dir;
            std::cout<<"Particle pos: " << "("<<particle_pos_in_plane.X()<<","<<particle_pos_in_plane.Y()<<","<<particle_pos_in_plane.Z()<<")"<<std::endl;
            std::cout<<"Counter pos: " << "("<<counter_pos.X()<<","<<counter_pos.Y()<<","<<counter_pos.Z()<<")"<<std::endl;
            */

          }
          /*
          geo::AuxDetGeo const*counter = fCounterSetPairs.front().setA.front();
          std::cout<<"Length: " << counter->Length() << std::endl;
          std::cout<<"HalfWidth1: " << counter->HalfWidth1() << std::endl;
          std::cout<<"HalfWidth2: " << counter->HalfWidth2() << std::endl;

          std::cout<<"HalfHeight: " << counter->HalfHeight() << std::endl;
          */
        }
      }
    }

    return false;
  }

  void GenFilter::beginJob() {
    art::ServiceHandle<geo::Geometry> geom;

    CounterSetPair EWCounterSetPair;
    CounterSetPair NupSdownCounterSetPair;
    CounterSetPair NdownSupCounterSetPair;

    for (unsigned int i = 0; i < geom->NAuxDets(); i++){
      //The WE counter pairs
      geo::AuxDetGeo const* counter = &(geom->AuxDet(i));
      if (i >=6 && i <= 15) EWCounterSetPair.setA.push_back(counter);
      else if (i >= 28 && i <=37) EWCounterSetPair.setB.push_back(counter);
      //The N (up) S (down) counter pairs
      else if (i >= 22 && i <= 27) NupSdownCounterSetPair.setA.push_back(counter);
      else if (i <= 5) NupSdownCounterSetPair.setB.push_back(counter);
      //The N (down) S (up) counter pairs
      else if (i >= 16 && i <= 21) NdownSupCounterSetPair.setA.push_back(counter);
      else if (i >= 38 && i <= 43) NdownSupCounterSetPair.setB.push_back(counter);
    }

    EWCounterSetPair.normalVec[0] = 0.;
    EWCounterSetPair.normalVec[1] = 0.;
    EWCounterSetPair.normalVec[2] = 1.;
    EWCounterSetPair.isRequested = fUseEWCounterPair;
    fCounterSetPairs.push_back(EWCounterSetPair);

    NupSdownCounterSetPair.normalVec[0] = 1.;
    NupSdownCounterSetPair.normalVec[1] = 0.;
    NupSdownCounterSetPair.normalVec[2] = 0.;
    NupSdownCounterSetPair.isRequested = fUseNupSdownCounterPair;
    fCounterSetPairs.push_back(NupSdownCounterSetPair);

    NdownSupCounterSetPair.normalVec[0] = 1.;
    NdownSupCounterSetPair.normalVec[1] = 0.;
    NdownSupCounterSetPair.normalVec[2] = 0.;
    NdownSupCounterSetPair.isRequested = fUseNdownSupCounterPair;
    fCounterSetPairs.push_back(NdownSupCounterSetPair);

    for (unsigned int i = 0; i < fCounterSetPairs.size(); i++){
      std::cout<<"Counter set pair: "<<i<<std::endl;
      std::cout<<"--setA size: " << fCounterSetPairs[i].setA.size() << std::endl;
      std::cout<<"--setB size: " << fCounterSetPairs[i].setB.size() << std::endl;

    }


  }

  bool GenFilter::IsInterestingParticle(const simb::MCParticle &particle){

    for (unsigned int i = 0; i < fInterestingPDGs.size(); i++){
      //Check if the particle under consideration has a requested PDG
      if (particle.PdgCode() == fInterestingPDGs[i]){
        //Got a requested PDG,  now check that the energy matches the requested range
        TLorentzVector mom4 = particle.Momentum();
        if (mom4.T() > fParticleMinEnergy && mom4.T() < fParticleMaxEnergy){
          //std::cout<<"FOUND INTERESTING PARTICLE"<<std::endl;
          return true;
        }
      }
    }

    return false;
  }

  bool GenFilter::ParticleHitsCounterSetPairs(const simb::MCParticle &particle, const CounterSetPair &CSP){
    //Need the normal to the counters
    TVector3 counter_norm(CSP.normalVec[0],CSP.normalVec[1],CSP.normalVec[2]);

    //Loop through one of the counter sets
    if (ParticleHitsCounterSet(particle, CSP.setA, counter_norm)){
      if (ParticleHitsCounterSet(particle, CSP.setB, counter_norm)){
        return true;
      }
    }
    return false;
  }

  bool GenFilter::ParticleHitsCounterSet(const simb::MCParticle &particle, const std::vector<geo::AuxDetGeo const*> &counters, const TVector3 &counter_norm){

    //Loop through the counters
    for (unsigned int i = 0; i < counters.size(); i++){
      //First step is to push the particle to the counter plane
      geo::AuxDetGeo const*counter = counters[i];

      double counter_pos_array[3];
      //fCounterSetPairs[CSP_i].setA.front()->GetCenter(counter_pos_array);
      counter->GetCenter(counter_pos_array);
      TVector3 counter_pos(counter_pos_array[0], counter_pos_array[1], counter_pos_array[2]);
      TVector3 particle_pos = particle.Position().Vect();
      TVector3 particle_dir = particle.Momentum().Vect().Unit();

      /*
          Length: 62.992
          HalfWidth1: 16.2814
          HalfWidth2: 13.5255
          HalfHeight: 0.475
          */

      double scale_factor = counter_norm.Dot(counter_pos - particle_pos)/(counter_norm.Dot(particle_dir));

      TVector3 particle_pos_in_plane = particle_pos + scale_factor * particle_dir;

      //We now have the particle position in the plane of the counter.  We now just need to calculate whether it sits inside the box
      //Create two TVector3s, each will hold the coordinates of opposing corners of the counter in question
      //We need to use the normals associated with the counter.  We already have two, one is a member of the C++ object and the other was passed to this function
      //Create the two TVector3s
      TVector3 pos_corner, neg_corner;
      //Lets start with the one passed to this function
      //The thin dimension of the counter is the one associated with normal passed to this function
      pos_corner += counter->HalfHeight()*counter_norm*fCounterSizeScaleFactor;
      neg_corner += -1.*counter->HalfHeight()*counter_norm*fCounterSizeScaleFactor; 

      //Now lets to the same for the normal stored in the C++ object (this "normal" actually points out the top of a counter)
      double counter_top_norm_array[3];
      counter->GetNormalVector(counter_top_norm_array);
      //Now package this up a TVector
      TVector3 counter_top_norm;
      counter_top_norm.SetX(counter_top_norm_array[0]);
      counter_top_norm.SetY(counter_top_norm_array[1]);
      counter_top_norm.SetZ(counter_top_norm_array[2]);
      //OK now we can add the dimension to the corner vectors.  The relevant dimension this time is Length/2
      pos_corner += counter->Length()*counter_top_norm*0.5*fCounterSizeScaleFactor;
      neg_corner += -1.*counter->Length()*counter_top_norm*0.5*fCounterSizeScaleFactor;
      //Now we need to the same for the vector pointing along the counter.  
      //The relevant dimension in this case in HalfWidth1 (going to assume the counter is a square and take the bigger width)
      //Because we have the other two vectors already, we can very easily get the final one by taking the cross product of them
      TVector3 counter_side_norm = counter_norm.Cross(counter_top_norm); 
      //now add the dimenions onto the corner vectors
      pos_corner += counter->HalfWidth1()*counter_side_norm*fCounterSizeScaleFactor;
      neg_corner += -1.*counter->HalfWidth1()*counter_side_norm*fCounterSizeScaleFactor;

      //Almost ready
      //We have calculated the corners assuming the centre of the counter is the origin.  Two choices, either translate the corners OR translate the propagated particle position
      //The latter is less lines of code so lets to that
      particle_pos_in_plane += -1.*counter_pos;

      //We are now ready to check if the particle sits in the counter
      //We don't know by default which way round the corners are oriented, but that doesn't matter as we can take abs, max and min
      //We are going to take abs of the particle position, the pos corner and the neg corner first
      /*
      pos_corner.SetXYZ(std::abs(pos_corner.X()),std::abs(pos_corner.Y()),std::abs(pos_corner.Z()));
      neg_corner.SetXYZ(std::abs(neg_corner.X()),std::abs(neg_corner.Y()),std::abs(neg_corner.Z()));
      particle_pos_in_plane.SetXYZ(std::abs(particle_pos_in_plane.X()),std::abs(particle_pos_in_plane.Y()),std::abs(particle_pos_in_plane.Z()));
      */

      /*
      //Dump the corners and particle pos
      std::cout<<"pos_corner: " << pos_corner.X()<<","<<pos_corner.Y()<<","<<pos_corner.Z()<<std::endl;
      std::cout<<"neg_corner: " << neg_corner.X()<<","<<neg_corner.Y()<<","<<neg_corner.Z()<<std::endl;
      std::cout<<"particle_pos_in_plane: " << particle_pos_in_plane.X()<<","<<particle_pos_in_plane.Y()<<","<<particle_pos_in_plane.Z()<<std::endl;
      */



      //Now we can check
      //This is going to be one huge if statement...
      if (particle_pos_in_plane.X() > std::min(pos_corner.X(), neg_corner.X()) && particle_pos_in_plane.X() < std::max(pos_corner.X(),neg_corner.X()) && particle_pos_in_plane.Y() > std::min(pos_corner.Y(), neg_corner.Y()) && particle_pos_in_plane.Y() < std::max(pos_corner.Y(),neg_corner.Y()) && particle_pos_in_plane.Z() > std::min(pos_corner.Z(), neg_corner.Z()) && particle_pos_in_plane.Z() < std::max(pos_corner.Z(),neg_corner.Z())){
        std::cout<<"Particle uses counter in set"<<std::endl;
        return true;
      }

    }

    return false;
  }


  DEFINE_ART_MODULE(GenFilter)

}
