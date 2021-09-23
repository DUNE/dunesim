// NoisyAdcSimulator.h
//
// PLasorak
// 05/2020
//
// Implementation of AdcSimulator for an noisy ADC with linear response.
//
// The ADC is characterized by its sensitivity Vsen (aka step size or LSB)
// and # bits Nbit.
// The sensitivity is the voltage difference between adjacent ADC bins.
// The voltage range is Vmax = Vsen*(2^Nbit-1) and the count is
//                      0 for Vin < 0.5*Vsen
//               2^Nbit-1 for Vin >= Vmax - 0.5*Vsen
//    int(Vin/Vsen + 0.5) otherwise

#ifndef NoisyAdcSimulator_H
#define NoisyAdcSimulator_H

#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "dune/DuneInterface/Data/AdcSimulator.h"

#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandBinomial.h"

class NoisyAdcSimulator : public AdcSimulator {

public:

  // Ctor from params.
  explicit NoisyAdcSimulator(double vsen, unsigned int nbit, int noise);

  // Ctor for  art tool.
  explicit NoisyAdcSimulator(fhicl::ParameterSet const& ps);

  // Evaluate an ADC count.
  Count count(double vin, Channel chan =0, Tick tick =0) const override;

private:
  int m_noise;
  double m_vsen;
  double m_vmax;
  Count m_adcmax;
  CLHEP::HepRandomEngine* m_random_engine;


};

#endif
