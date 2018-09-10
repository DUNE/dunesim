////////////////////////////////////////////////////////////////////////
// Class:       CosmicFilter
// Plugin Type: filter (art v2_11_03)
// File:        CosmicFilter_module.cc
//
// Generated at Fri Aug 24 06:12:20 2018 by Philip Rodrigues using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcore/Geometry/Geometry.h"

#include <memory>

#include "TGeoManager.h"
#include "TLorentzVector.h"
#include "TVector3.h"


class CosmicFilter;


class CosmicFilter : public art::EDFilter {
public:
    explicit CosmicFilter(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CosmicFilter(CosmicFilter const &) = delete;
    CosmicFilter(CosmicFilter &&) = delete;
    CosmicFilter & operator = (CosmicFilter const &) = delete;
    CosmicFilter & operator = (CosmicFilter &&) = delete;

    // Required functions.
    bool filter(art::Event & e) override;

private:

    // Declare member data here.
    
    // Project the particle in a straight line: does it intersect the
    // volume with the given name?
    bool intersectsVolume(const simb::MCParticle& particle, const char* volumeName) const;

    art::ServiceHandle<geo::Geometry> geom;
};


CosmicFilter::CosmicFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
    // Call appropriate produces<>() functions here.
}

bool CosmicFilter::filter(art::Event & e)
{
    auto const& mctruths =
        *e.getValidHandle<std::vector<simb::MCTruth>>(art::InputTag{"generator"});
    
    const simb::MCTruth& truth=mctruths.at(0);
    const simb::MCParticle& cosmic_particle=truth.GetParticle(0);
    return intersectsVolume(cosmic_particle, "volDetEnclosure_0");
}

bool CosmicFilter::intersectsVolume(const simb::MCParticle& particle, const char* volumeName) const
{
    TGeoManager* geoManager=geom->ROOTGeoManager();
    geoManager->SetCurrentPoint(particle.Vx(), particle.Vy(), particle.Vz());
    TVector3 unitMom=particle.Momentum().Vect().Unit();
    geoManager->SetCurrentDirection(unitMom.X(), unitMom.Y(), unitMom.Z());
    geoManager->FindNextBoundaryAndStep();

    while(!geoManager->IsOutside()){
        geoManager->FindNextBoundaryAndStep();
        if(TString(geoManager->GetPath()).Contains(volumeName)) return true;
    }
    return false;
}

DEFINE_ART_MODULE(CosmicFilter)
