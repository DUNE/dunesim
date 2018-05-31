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

#include <memory>

class NoiseAdder;


class NoiseAdder : public art::EDProducer {
public:
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

};


NoiseAdder::NoiseAdder(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
    produces< std::vector<raw::RawDigit> >();
}

void NoiseAdder::produce(art::Event & e)
{
    CLHEP::HepJamesRandom rand;
    CLHEP::RandGauss gaus(rand, 5, 0);
    // Implementation of required member function here.
    auto& digits_in =
        *e.getValidHandle<std::vector<raw::RawDigit>>("daq");
    std::unique_ptr<std::vector<raw::RawDigit>> digits_out(new std::vector<raw::RawDigit>);

    for(auto&& digit: digits_in){
        if(digit.Compression()!=0){
            std::cout << "Compression type " << digit.Compression() << std::endl;
        }
        std::vector<short> samples_out(digit.NADC(), 0);
        for(size_t i=0; i<digit.NADC(); ++i){
            samples_out[i]=digit.ADC(i)+gaus.fire();
        }
        digits_out->push_back(raw::RawDigit(digit.Channel(),
                                           digit.Samples(),
                                           std::move(samples_out),
                                           raw::kNone));
    } // end loop over digits (=?channels)

    e.put(std::move(digits_out));
}

DEFINE_ART_MODULE(NoiseAdder)
