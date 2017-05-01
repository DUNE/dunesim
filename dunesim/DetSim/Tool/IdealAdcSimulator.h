// IdealAdcSimulator.h
//
// David Adams
// April 2017
//
// Implementation of AdcSimulator for an ideal ADC with linear response.
//
// The ADC is characterized by its sensitivity Vsen (aka step size or LSB)
// and # bits Nbit.
// The sensitivity is the voltage difference between adjacent ADC bins.
// The voltage range is Vmax = Vsen*(2^Nbit-1) and the count is
//                      0 for Vin < 0.5*Vsen
//               2^Nbit-1 for Vin >= Vmax - 0.5*Vsen
//    int(Vin/Vsen + 0.5) otherwise

#ifndef IdealAdcSimulator_H
#define IdealAdcSimulator_H

#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "dune/DuneInterface/AdcSimulator.h"

class IdealAdcSimulator : public AdcSimulator {

public:

  // Ctor from params.
  explicit IdealAdcSimulator(double vsen, unsigned int nbit =12);

  // Ctor for  art tool.
  explicit IdealAdcSimulator(fhicl::ParameterSet const& ps);

  // Evaluate an ADC count.
  Count count(double vin, Channel chan =0, Tick tick =0) const override;

private:

  double m_vsen;
  double m_vmax;
  Count m_adcmax;

};

DEFINE_ART_CLASS_TOOL(IdealAdcSimulator)

#endif
