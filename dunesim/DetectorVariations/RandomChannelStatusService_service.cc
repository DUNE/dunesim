// Chris Backhouse - c.backhouse@ucl.ac.uk Dec 2017

#include "dune/DetectorVariations/RandomChannelStatusService.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/sim.h" // GetRandomNumberSeed()

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "CLHEP/Random/RandFlat.h"

namespace detvar
{
  // Implement Table 5 from dune docdb 4064
  void ChipAndChannelToWire(int chip, int chan,
                            geo::View_t& view,
                            int& wire)
  {
    assert(chip >= 1 && chip <= 8);
    assert(chan >= 0 && chan <= 15);

    int wireMin; // smallest value in block
    int chipMin, chanMin; // position of smallest value
    int height; // height of block
    int chipSign, chanSign; // which way the numbers count

    if(chip == 1 || chip == 2){
      wireMin = 1; chipMin = 2; chipSign = -1; chanSign = -1;
      /**/ if(chan <= 4){view = geo::kU; chanMin =  4; height = 5;}
      else if(chan <= 9){view = geo::kV; chanMin =  9; height = 5;}
      else/*          */{view = geo::kW; chanMin = 15; height = 6;}
    }

    if(chip == 3 || chip == 4){
      wireMin = 2; chipMin = 4; chipSign = -1; chanSign = +1;
      /**/ if(chan <= 5 ){view = geo::kW; chanMin =  0; height = 6;}
      else if(chan <= 10){view = geo::kV; chanMin =  6; height = 5;}
      else/*           */{view = geo::kU; chanMin = 11; height = 5;}
    }

    if(chip == 5 || chip == 6){
      chipMin = 5; chipSign = +1; chanSign = -1;
      /**/ if(chan <= 5){view = geo::kU; wireMin = 21; chanMin =  4; height = 5;}
      else if(chan <= 9){view = geo::kV; wireMin = 21; chanMin =  9; height = 5;}
      else/*          */{view = geo::kW; wireMin = 25; chanMin = 15; height = 6;}
    }

    if(chip == 7 || chip == 8){
      chipMin = 7; chipSign = +1; chanSign = +1;
      /**/ if(chan <= 6 ){view = geo::kW; wireMin = 26; chanMin = 0; height = 6;}
      else if(chan <= 11){view = geo::kV; wireMin = 22; chanMin = 6; height = 5;}
      else/*           */{view = geo::kU; wireMin = 22; chanMin = 11; height = 5;}
    }

    // Compute the wire number inside the block
    wire = wireMin+2*((chan-chanMin)*chanSign + (chip-chipMin)*height*chipSign);

    assert(wire >= 1 && wire <= 48);
  }


  RandomChannelStatusProvider::
  RandomChannelStatusProvider(const fhicl::ParameterSet& pset)
  {
    const double badfrac = pset.get<double>("BadChanFrac");

    art::ServiceHandle<geo::Geometry> geom;

    //    const unsigned int seed = pset.get<unsigned int>("Seed", sim::GetRandomNumberSeed());
    //    createEngine(seed);
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::RandFlat r(rng->getEngine());

    // Geometry doesn't have a way to iterate directly over channels. Iterate
    // over the wires and convert them. Use a set to remove duplicates
    std::set<raw::ChannelID_t> allchans;
    for(geo::WireID wire: geom->IterateWireIDs())
      allchans.insert(geom->PlaneWireToChannel(wire));

    // But a vector is much easier to pick from randomly
    const std::vector<raw::ChannelID_t> vchans(allchans.begin(),
                                               allchans.end());
    const int N = vchans.size();

    // Generate exactly the requested fraction of bad channels (rather than a
    // random sample with that probability). Should make results of studies
    // less noisy.
    while(fBadChans.size() < badfrac*N){
      // Insert a random element. There will be duplicates, but the set will
      // filter them out. Shouldn't be too inefficient for the low bad channel
      // fractions we'll use in practice.
      fBadChans.insert(vchans[r.shootInt(N)]);
    }

    // goodchans = allchans - badchans
    std::set_difference(allchans.begin(), allchans.end(),
                        fBadChans.begin(), fBadChans.end(),
                        std::inserter(fGoodChans, fGoodChans.begin()));
  }
}

DEFINE_ART_SERVICE_INTERFACE_IMPL(detvar::RandomChannelStatusService, lariov::ChannelStatusService)
