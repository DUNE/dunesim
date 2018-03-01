// Chris Backhouse - c.backhouse@ucl.ac.uk Dec 2017

#include "dune/DetectorVariations/RandomChannelStatusService.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/sim.h" // GetRandomNumberSeed()

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "CLHEP/Random/RandFlat.h"

namespace detvar
{
  //......................................................................
  // Implement Table 5 from dune docdb 4064v2
  void ChipAndChannelToSpot(int chip, int chan, geo::View_t& view, int& wire)
  {
    // Ultimately this logic should be contained in some ChannelMap service

    assert(chip >= 1 && chip <= 8);
    assert(chan >= 0 && chan <= 15);

    int wireMin = -1; // smallest value in block
    int chipMin = -1, chanMin = -1; // position of smallest value
    int height = -1; // height of block
    int chipSign = 0, chanSign = 0; // which way the numbers count

    if(chip == 1 || chip == 2){
      wireMin = 1; chipMin = 2; chipSign = -1; chanSign = -1;
      /**/ if(chan <= 4){view = geo::kU; chanMin =  4; height = 5;}
      else if(chan <= 9){view = geo::kV; chanMin =  9; height = 5;}
      else/*          */{view = geo::kW; chanMin = 15; height = 6;}
    }

    if(chip == 3 || chip == 4){
      wireMin = 2; chipMin = 4; chipSign = -1; chanSign = +1;
      /**/ if(chan <=  5){view = geo::kW; chanMin =  0; height = 6;}
      else if(chan <= 10){view = geo::kV; chanMin =  6; height = 5;}
      else/*           */{view = geo::kU; chanMin = 11; height = 5;}
    }

    if(chip == 5 || chip == 6){
      chipMin = 5; chipSign = +1; chanSign = -1;
      /**/ if(chan <= 4){view = geo::kU; wireMin = 21; chanMin =  4; height = 5;}
      else if(chan <= 9){view = geo::kV; wireMin = 21; chanMin =  9; height = 5;}
      else/*          */{view = geo::kW; wireMin = 25; chanMin = 15; height = 6;}
    }

    if(chip == 7 || chip == 8){
      chipMin = 7; chipSign = +1; chanSign = +1;
      /**/ if(chan <=  5){view = geo::kW; wireMin = 26; chanMin =  0; height = 6;}
      else if(chan <= 10){view = geo::kV; wireMin = 22; chanMin =  6; height = 5;}
      else/*           */{view = geo::kU; wireMin = 22; chanMin = 11; height = 5;}
    }

    assert(wireMin >= 0 && chipMin >= 0 && chanMin >= 0 && height >= 0 && chipSign != 0 && chanSign != 0);

    // Compute the wire number inside the block
    wire = wireMin+2*((chan-chanMin)*chanSign + (chip-chipMin)*height*chipSign);

    assert(wire >= 1 && wire <= 48);
  }

  /// For whatever reason the iterators in Geometry can't be used to initialize
  /// a vector directly
  template<class T, class It_t> std::vector<T> VectorViaSet(It_t begin,
                                                            It_t end)
  {
    const std::set<T> s(begin, end);
    return std::vector<T>(s.begin(), s.end());
  }

  class RandomTPC
  {
  public:
    RandomTPC(const geo::GeometryCore* geom)
      : fCryos(VectorViaSet<geo::CryostatID>(geom->begin_cryostat_id(),
                                             geom->end_cryostat_id()))
    {
      for(geo::CryostatID c: fCryos){
        fTPCs[c] = VectorViaSet<geo::TPCID>(geom->begin_TPC_id(c),
                                            geom->end_TPC_id(c));

        fTPCsets[c] = VectorViaSet<readout::TPCsetID>(geom->begin_TPCset_id(c),
                                                      geom->end_TPCset_id(c));
      }
    }

    geo::TPCID GetTPC() const
    {
      const geo::CryostatID c = fCryos[gRandom->Integer(fCryos.size())];
      const std::vector<geo::TPCID>& ts = fTPCs.find(c)->second;
      return ts[gRandom->Integer(ts.size())];
    }

    readout::TPCsetID GetTPCset() const
    {
      const geo::CryostatID c = fCryos[gRandom->Integer(fCryos.size())];
      const std::vector<readout::TPCsetID>& ts = fTPCsets.find(c)->second;
      return ts[gRandom->Integer(ts.size())];
    }

  protected:
    std::vector<geo::CryostatID> fCryos;
    std::map<geo::CryostatID, std::vector<geo::TPCID>> fTPCs;
    std::map<geo::CryostatID, std::vector<readout::TPCsetID>> fTPCsets;
  };

  class SortChansByZ
  {
  public:
    bool operator()(const raw::ChannelID_t& a,
                    const raw::ChannelID_t& b) const
    {
      return GetZ(a) < GetZ(b);
    }
  protected:
    double GetZ(const raw::ChannelID_t& c) const
    {
      double z = 0;
      double y = -1e9;

      // Search through all the wires associated with this channel for the
      // top-most position (where it attaches to the APA frame) and quote that
      // Z position.
      const std::vector<geo::WireID> wires = geom->ChannelToWire(c);
      for(const geo::WireID& wire: wires){
        const geo::WireGeo& wg = geom->Wire(wire);
        if(wg.GetStart().Y() > y){
          y = wg.GetStart().Y();
          z = wg.GetStart().Z();
        }
        if(wg.GetEnd().Y() > y){
          y = wg.GetEnd().Y();
          z = wg.GetEnd().Z();
        }
      }

      return z;
    }

    art::ServiceHandle<geo::Geometry> geom;
  };

  //......................................................................
  // Three vectors, one for each view. Will be sorted by ID number
  std::vector<std::vector<raw::ChannelID_t>>
  ChannelsForTPC(const geo::GeometryCore* geom, geo::TPCID tpc)
  {
    // No good way of enumerating all the unique channels, or selecting them by
    // index. Have to figure them out from the wires.
    std::set<raw::ChannelID_t> chanset[3];
    for(geo::WireID wire: geom->IterateWireIDs(tpc)){
      const raw::ChannelID_t chan = geom->PlaneWireToChannel(wire);
      // But this also gives us wires that are actually attached to the other
      // face and just wrapped onto this face. So long as the order of the
      // vector returned from this function is meaningful, this should work
      // to keep just the ones we need.
      if(geom->ChannelToWire(chan)[0] == wire)
        chanset[geom->View(chan)].insert(chan);
    }
    // We'll want to index by channel number within the APAs. Sort spatially as
    // our best guess as to that mapping.
    std::vector<std::vector<raw::ChannelID_t>> chans;
    for(int i = 0; i < 3; ++i){
      chans.emplace_back(chanset[i].begin(), chanset[i].end());
      SortChansByZ scz;
      std::sort(chans.back().begin(), chans.back().end(), scz);
    }
    return chans;
  }

  //......................................................................
  RandomChannelStatusProvider::
  RandomChannelStatusProvider(const fhicl::ParameterSet& pset)
  {
    const double badfrac = pset.get<double>("BadChanFrac");

    EMode_t mode = kUnknown;

    const std::string modestr = pset.get<std::string>("Mode");
    if(modestr == "channels") mode = kRandomChans;
    if(modestr == "APAs")     mode = kRandomAPAs;
    if(modestr == "APAsides") mode = kRandomAPAsides;
    if(modestr == "boards")   mode = kRandomBoards;
    if(modestr == "chips")    mode = kRandomChips;
    if(mode == kUnknown){
      std::cout << "RandomChannelStatusService: unknown mode '"
                << modestr << "'" << std::endl;
      abort();
    }

    art::ServiceHandle<geo::Geometry> geom;

    // TODO don't seem to be able to access the RandomNumberGenerator service
    // from here, let alone anything to do with seeds.
    //    const unsigned int seed = pset.get<unsigned int>("Seed", sim::GetRandomNumberSeed());
    //    createEngine(seed);
    //    art::ServiceHandle<art::RandomNumberGenerator> rng;
    //    CLHEP::RandFlat r(rng->getEngine());

    const unsigned int target = badfrac*geom->Nchannels();

    switch(mode){
    case kRandomChans:    MarkChansBad(target);    break;

    case kRandomAPAs:     MarkAPAsBad(target);     break;

    case kRandomAPAsides: MarkAPASidesBad(target); break;

    case kRandomBoards:
    case kRandomChips:
      MarkBoardsOrChipsBad(mode, target);
      break;

    default:
      abort(); // impossible
    }

    // goodchans = allchans - badchans
    std::set<raw::ChannelID_t> allchans;
    for(geo::WireID wire: geom->IterateWireIDs())
      allchans.insert(geom->PlaneWireToChannel(wire));

    std::set_difference(allchans.begin(), allchans.end(),
                        fBadChans.begin(), fBadChans.end(),
                        std::inserter(fGoodChans, fGoodChans.begin()));
  }

  //......................................................................
  void RandomChannelStatusProvider::MarkChansBad(unsigned int target)
  {
    art::ServiceHandle<geo::Geometry> geom;

    // Geometry doesn't have a way to iterate directly over channels. Iterate
    // over the wires and convert them. Use a set to remove duplicates
    std::set<raw::ChannelID_t> allchans;
    for(geo::WireID wire: geom->IterateWireIDs())
      allchans.insert(geom->PlaneWireToChannel(wire));

    // But a vector is much easier to pick from randomly. This will be used for
    // the random chans mode.
    const std::vector<raw::ChannelID_t> vchans(allchans.begin(),
                                               allchans.end());

    // Generate exactly the requested fraction of bad channels (rather than a
    // random sample with that probability). Should make results of studies
    // less noisy.
    while(fBadChans.size() < target){
      // Insert a random element. There will be duplicates, but the set will
      // filter them out. Shouldn't be too inefficient for the low bad
      // channel fractions we'll use in practice.
      fBadChans.insert(vchans[gRandom->Integer(vchans.size())]);
    } // end while
  }

  //......................................................................
  void RandomChannelStatusProvider::MarkAPAsBad(unsigned int target)
  {
    art::ServiceHandle<geo::Geometry> geom;

    std::map<readout::TPCsetID, std::set<raw::ChannelID_t>> tpcset_to_chans;
    for(const readout::TPCsetID& ts: geom->IterateTPCsetIDs()){
      // There is no version of IterateWireIDs over a TPCset. Use another layer
      // of indirection.
      for(geo::TPCID t: geom->TPCsetToTPCs(ts)){
        for(const geo::WireID& wire: geom->IterateWireIDs(t)){
          const raw::ChannelID_t chan = geom->PlaneWireToChannel(wire);
          tpcset_to_chans[ts].insert(chan);
        }
      }
    }

    RandomTPC tpcs(geom.get());

    while(fBadChans.size() < target){
      // Mark all the channels in this TPCset (both sides of an APA) bad
      const readout::TPCsetID t = tpcs.GetTPCset();
      for(raw::ChannelID_t chan: tpcset_to_chans[t]){
        fBadChans.insert(chan);
      }

      std::cout << "RandomChannelStatusService: Generated "
                << fBadChans.size() << " bad channels of "
                << target << " required" << std::endl;
    }
  }

  //......................................................................
  void RandomChannelStatusProvider::MarkAPASidesBad(unsigned int target)
  {
    art::ServiceHandle<geo::Geometry> geom;

    std::map<geo::TPCID, std::vector<raw::ChannelID_t>> tpc_to_chans;
    for(const geo::TPCID& t: geom->IterateTPCIDs()){
      // Geometry doesn't provide a way to directly iterate the channels in the
      // TPC. Instead iterate the wires and convert to channels
      for(const geo::WireID& wire: geom->IterateWireIDs(t)){
        const raw::ChannelID_t chan = geom->PlaneWireToChannel(wire);
        // But this also gives us wires that are actually attached to the other
        // face and just wrapped onto this face. So long as the order of the
        // vector returned from this function is meaningful, this should work
        // to keep just the ones we need.
        if(geom->ChannelToWire(chan)[0] == wire)
          tpc_to_chans[t].push_back(chan);
      }
    }

    RandomTPC tpcs(geom.get());

    while(fBadChans.size() < target){
      const geo::TPCID t = tpcs.GetTPC();

      // Mark all the channels in this TPC (ie APA side) bad
      for(raw::ChannelID_t chan: tpc_to_chans[t]){
        fBadChans.insert(chan);
      }

      std::cout << "RandomChannelStatusService: Generated "
                << fBadChans.size() << " bad channels of "
                << target << " required" << std::endl;
    }
  }

  //......................................................................
  void RandomChannelStatusProvider::
  MarkBoardsOrChipsBad(EMode_t mode, unsigned int target)
  {
    art::ServiceHandle<geo::Geometry> geom;

    RandomTPC tpcs(geom.get());

    // Generate exactly the requested fraction of bad channels (rather than a
    // random sample with that probability). Should make results of studies
    // less noisy.
    while(fBadChans.size() < target){
      const geo::TPCID t = tpcs.GetTPC();

      const std::vector<std::vector<raw::ChannelID_t>> chans = ChannelsForTPC(geom.get(), t);

      // An APA has 20 boards total, 10 on each side. ie a side has 480 W wires
      // and 400+400 U+V wires.
      const int board = gRandom->Integer(10);

      if(mode == kRandomBoards) MarkBoardBad(board, chans);

      if(mode == kRandomChips){
        const int chip = 1+gRandom->Integer(8); // One of 8 chips
        MarkChipBad(board, chip, geom.get(), chans);
      }

      std::cout << "RandomChannelStatusService: Generated "
                << fBadChans.size() << " bad channels of "
                << target << " required" << std::endl;
    }
  }

  //......................................................................
  void RandomChannelStatusProvider::
  MarkBoardBad(int board,
               const std::vector<std::vector<raw::ChannelID_t>>& chans)
  {
    // Check we succesfully got a single side of an APA
    assert(chans[geo::kU].size() == 400);
    assert(chans[geo::kV].size() == 400);
    assert(chans[geo::kW].size() == 480);

    for(int i = 0; i < 40; ++i) fBadChans.insert(chans[geo::kU][board*40+i]);
    for(int i = 0; i < 40; ++i) fBadChans.insert(chans[geo::kV][board*40+i]);
    for(int i = 0; i < 48; ++i) fBadChans.insert(chans[geo::kW][board*48+i]);
  }

  //......................................................................
  void RandomChannelStatusProvider::
  MarkChipBad(int board, int chip,
              const geo::GeometryCore* geom,
              const std::vector<std::vector<raw::ChannelID_t>>& chans)
  {
    // Check we succesfully got a single side of an APA
    assert(chans[geo::kU].size() == 400);
    assert(chans[geo::kV].size() == 400);
    assert(chans[geo::kW].size() == 480);

    // Knock out all the channels in this chip
    for(int chan = 0; chan <= 15; ++chan){
      geo::View_t view = geo::kUnknown;
      // TODO the table calls this a "spot". We're using it as an index into
      // the sorted list of channel IDs for this view (by Z coordinate). That
      // could be wrong, but the Geometry doesn't seem to have any coresponding

      // concept.
      int spot;
      ChipAndChannelToSpot(chip, chan, view, spot);

      // How many wires to add on between each board
      const int stride = (view == geo::kW) ? 48 : 40;

      assert(spot+board*stride <= int(chans[view].size()));
      fBadChans.insert(chans[view][spot+board*stride-1]);
    } // end for chan
  }
}

DEFINE_ART_SERVICE_INTERFACE_IMPL(detvar::RandomChannelStatusService, lariov::ChannelStatusService)
