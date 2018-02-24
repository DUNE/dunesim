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

  //......................................................................
  geo::CryostatID RandomCryostat(const geo::GeometryCore* geom)
  {
    // For whatever reason the iterators here can't be used to initialize a
    // vector.
    const std::set<geo::CryostatID> css(geom->begin_cryostat_id(),
                                        geom->end_cryostat_id());
    const std::vector<geo::CryostatID> csv(css.begin(), css.end());
    return csv[gRandom->Integer(csv.size())];
  }

  //......................................................................
  geo::TPCID RandomTPC(const geo::GeometryCore* geom, geo::CryostatID c)
  {
    // For whatever reason the iterators here can't be used to initialize a
    // vector.
    const std::set<geo::TPCID> tss(geom->begin_TPC_id(c), geom->end_TPC_id(c));
    const std::vector<geo::TPCID> tsv(tss.begin(), tss.end());
    return tsv[gRandom->Integer(tsv.size())];
  }

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
      chanset[geom->View(chan)].insert(chan);
    }
    // Then we need to index by channel number within the APA. TODO: I
    // have no idea if the sorting of the channel IDs that the set did is
    // what we need. Maybe U and V get sorted in opposite order to each
    // other etc?
    std::vector<std::vector<raw::ChannelID_t>> chans;
    for(int i = 0; i < 3; ++i){
      chans.emplace_back(chanset[i].begin(), chanset[i].end());
    }
    return chans;
  }

  //......................................................................
  RandomChannelStatusProvider::
  RandomChannelStatusProvider(const fhicl::ParameterSet& pset)
  {
    const double badfrac = pset.get<double>("BadChanFrac");

    enum{
      kUnknown, kRandomChans, kRandomAPAs, kRandomBoards, kRandomChips
    } mode = kUnknown;

    const std::string modestr = pset.get<std::string>("Mode");
    if(modestr == "channels") mode = kRandomChans;
    if(modestr == "APAs")     mode = kRandomAPAs;
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

    // Geometry doesn't have a way to iterate directly over channels. Iterate
    // over the wires and convert them. Use a set to remove duplicates
    std::set<raw::ChannelID_t> allchans;
    for(geo::WireID wire: geom->IterateWireIDs())
      allchans.insert(geom->PlaneWireToChannel(wire));

    // But a vector is much easier to pick from randomly. This will be used for
    // the random chans mode.
    const std::vector<raw::ChannelID_t> vchans(allchans.begin(),
                                               allchans.end());
    const int N = vchans.size();

    // Generate exactly the requested fraction of bad channels (rather than a
    // random sample with that probability). Should make results of studies
    // less noisy.
    while(fBadChans.size() < badfrac*N){

      if(mode == kRandomChans){
        // Insert a random element. There will be duplicates, but the set will
        // filter them out. Shouldn't be too inefficient for the low bad
        // channel fractions we'll use in practice.
        fBadChans.insert(vchans[gRandom->Integer(N)]);
      }
      else{ // by APAs, boards, or chips

        const geo::CryostatID c = RandomCryostat(geom.get());
        const geo::TPCID t = RandomTPC(geom.get(), c);

        if(mode == kRandomAPAs){
          // Mark all the channels in this TPC bad
          for(geo::WireID wire: geom->IterateWireIDs(t)){
            fBadChans.insert(geom->PlaneWireToChannel(wire));
          }
        }
        else{ // boards or chips
          // We need to index by channel number within the APA. TODO: I have no
          // idea if the sorting by the channel IDs is what we need. Maybe U
          // and V get sorted in opposite order to each other etc?
          const std::vector<std::vector<raw::ChannelID_t>> chans = ChannelsForTPC(geom.get(), t);

          // Empirically we need to repeat the organization from the table 10
          // times to readout the whole APA (480 W wires and 400+400 U+V
          // wires).
          const int board = gRandom->Integer(10);

          if(mode == kRandomBoards) MarkBoardBad(board, chans);

          if(mode == kRandomChips){
            const int chip = 1+gRandom->Integer(8); // One of 8 chips
            MarkChipBad(board, chip, geom.get(), chans);
          }
        }

        // This procedure with all the enumerations and sets and so on is
        // slow. Give the user some feedback, and an incentive to fix it.
        std::cout << "RandomChannelStatusService: Generated " << fBadChans.size() << " bad channels of " << badfrac*N << " required" << std::endl;
      }
    } // end while

    // goodchans = allchans - badchans
    std::set_difference(allchans.begin(), allchans.end(),
                        fBadChans.begin(), fBadChans.end(),
                        std::inserter(fGoodChans, fGoodChans.begin()));
  }

  //......................................................................
  void RandomChannelStatusProvider::
  MarkBoardBad(int board,
               const std::vector<std::vector<raw::ChannelID_t>>& chans)
  {
    for(int i = 0; i < 48; ++i) fBadChans.insert(chans[geo::kW][board*48+i]);
    for(int i = 0; i < 40; ++i) fBadChans.insert(chans[geo::kU][board*40+i]);
    for(int i = 0; i < 40; ++i) fBadChans.insert(chans[geo::kV][board*40+i]);
  }

  //......................................................................
  void RandomChannelStatusProvider::
  MarkChipBad(int board, int chip,
              const geo::GeometryCore* geom,
              const std::vector<std::vector<raw::ChannelID_t>>& chans)
  {
    // Knock out all the channels in this chip
    for(int chan = 0; chan <= 15; ++chan){
      // TODO the table we're using just says "view W" with no statement about
      // which side of the APA it is, since it's a ProtoDUNE table. We ignore
      // that complication and I think in practice render one side dead
      // arbitrarily.
      geo::View_t view = geo::kUnknown;
      // TODO the table calls this a "spot". We're using it as an index
      // into the sorted list of channel IDs for this view. That could be
      // wrong, but the Geometry doesn't seem to have any coresponding
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
