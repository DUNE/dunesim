// Chris Backhouse - c.backhouse@ucl.ac.uk Dec 2017

#include "dune/DetectorVariations/RandomChannelStatusService.h"

#include "larcore/Geometry/Geometry.h"

#include "TRandom3.h"

namespace detvar
{
  RandomChannelStatusProvider::
  RandomChannelStatusProvider(const fhicl::ParameterSet& pset)
  {
    const double badfrac = pset.get<double>("BadChanFrac");

    art::ServiceHandle<geo::Geometry> geom;

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
      fBadChans.insert(vchans[gRandom->Integer(N)]);
    }

    // goodchans = allchans - badchans
    std::set_difference(allchans.begin(), allchans.end(),
                        fBadChans.begin(), fBadChans.end(),
                        std::inserter(fGoodChans, fGoodChans.begin()));
  }
}

DEFINE_ART_SERVICE_INTERFACE_IMPL(detvar::RandomChannelStatusService, lariov::ChannelStatusService)
