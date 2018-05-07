// AdcCodeHelper.h
//
// David Adams
// November 2015
//
// Utility class to help interprete ADC signals.
// Parameters:
//   SupressedSignalMax - Signals with |sig-ped| >= this are not suppressed.
//                        If 0, all are suppressed
//
// Legacy behavior, i.e. that for ADCStickyCodeCheck in
// dunetpc v04.29.02 is recovered with parameter value
//   SuppressedSignalMax = 64

#ifndef AdcCodeHelper_H
#define AdcCodeHelper_H

#include "dune/DuneInterface/AdcTypes.h"

class AdcCodeHelper {

public:

  // Ctor.
  AdcCodeHelper(AdcCount aSuppressedSignalMax =0);

  // Return if the sticky bits are set.
  bool hasStickyBits(AdcCount sig);

  // Return if |pedsig| < SuppressedSignalMax.
  // Returns true if SuppressedSignalMax <= 0.
  bool isSmall(AdcCount sig, AdcSignal ped =0.0);

  // Convert a float signal to an integer signal.
  AdcCount intSignal(AdcSignal fsig) const;

  // Return a pedestal-subtracted signal.
  AdcSignal subtract(AdcCount sig, AdcSignal ped) const;

  // Return a pedestal-subtracted signal as an integer.
  AdcCount intSubtract(AdcCount sig, AdcSignal ped) const;

private:

  // Properties.
  //AdcCount m_SupresssedValue; // unused
  AdcSignal m_SuppressedSignalMax;

};
  
#endif
