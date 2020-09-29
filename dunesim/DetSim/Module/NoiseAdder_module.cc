////////////////////////////////////////////////////////////////////////
// Class:       NoiseAdder
// Plugin Type: producer (art v2_10_03)
// File:        NoiseAdder_module.cc
//
// Generated at Thu May 31 08:35:23 2018 by Philip Rodrigues using cetskelgen
// from cetlib version v3_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/JamesRandom.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "dune/DuneInterface/ChannelNoiseService.h"

#include "cetlib/search_path.h"

#include <memory>

class NoiseAdder;


class NoiseAdder : public art::EDProducer {
public:
    struct ElectronicsAddress
    {
        ElectronicsAddress(int _crate, int _slot, int _fiber, int _asic, int _asicChannel)
            : crate(_crate), slot(_slot), fiber(_fiber), asic(_asic), asicChannel(_asicChannel)
            {}
        int crate;
        int slot;
        int fiber;
        int asic;
        int asicChannel;
    };

    explicit NoiseAdder(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    NoiseAdder(NoiseAdder const &) = delete;
    NoiseAdder(NoiseAdder &&) = delete;
    NoiseAdder & operator = (NoiseAdder const &) = delete;
    NoiseAdder & operator = (NoiseAdder &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;


private:

    // Declare member data here.
    art::ServiceHandle<ChannelNoiseService> m_noiseService;
    std::map<int, ElectronicsAddress> m_channelMap;
};


NoiseAdder::NoiseAdder(fhicl::ParameterSet const & p) : EDProducer{p}

// Initialize member data here.
{
    produces< std::vector<raw::RawDigit> >();

    // The ProtoDUNE channel map lives in dune_raw_data:
    //
    // https://cdcvs.fnal.gov/redmine/projects/dune-raw-data/repository/revisions/develop/entry/dune-raw-data/Services/ChannelMap/protoDUNETPCChannelMap_v3.txt
    //
    // It's created by:
    //
    // https://cdcvs.fnal.gov/redmine/projects/dune-raw-data/repository/revisions/develop/entry/dune-raw-data/Services/ChannelMap/mapmakers/MakePdspChannelMap_v3.C
    cet::search_path sp("FW_SEARCH_PATH");
    std::string channelMapFile=sp.find_file("protoDUNETPCChannelMap_v3.txt");
    std::ifstream fin(channelMapFile.c_str());

    int crateNo, slotNo, fiberNo, FEMBChannel,
        StreamChannel, slotID, fiberID,
        chipNo, chipChannel, asicNo,
        asicChannel, planeType, offlineChannel;
    
    while(fin >> crateNo >> slotNo >> fiberNo >> FEMBChannel
          >> StreamChannel >> slotID >> fiberID
          >> chipNo >> chipChannel >> asicNo
          >> asicChannel >> planeType >> offlineChannel){
        m_channelMap.emplace(std::make_pair(offlineChannel,
                                            ElectronicsAddress(crateNo, slotNo, fiberNo, asicNo, asicChannel)));
    }
}

void NoiseAdder::produce(art::Event & e)
{
    // Implementation of required member function here.
    auto& digits_in =
        *e.getValidHandle<std::vector<raw::RawDigit>>("daq");
    std::unique_ptr<std::vector<raw::RawDigit>> digits_out(new std::vector<raw::RawDigit>);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);

    for(auto&& digit: digits_in){
        if(digit.Compression()!=0){
            // TODO: throw or just uncompress the stream and carry on
            std::cout << "Compression type " << digit.Compression() << std::endl;
        }
        std::vector<short> samples_out(digit.NADC(), 0);
        std::vector<float> samples_work(digit.NADC(), 0);
        for(size_t i=0; i<digit.NADC(); ++i){
            samples_out[i]=m_noiseService->addNoise(clockData, detProp, digit.Channel(), samples_work);
        }
        digits_out->push_back(raw::RawDigit(digit.Channel(),
                                           digit.Samples(),
                                           std::move(samples_out),
                                           raw::kNone));
    } // end loop over digits (=?channels)

    e.put(std::move(digits_out));
}

DEFINE_ART_MODULE(NoiseAdder)
