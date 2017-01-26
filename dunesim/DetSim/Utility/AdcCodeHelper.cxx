// AdcCodeHelper.cxx

#include "dune/DetSim/Utility/AdcCodeHelper.h"
#include <cmath>

//**********************************************************************

AdcCodeHelper::AdcCodeHelper(AdcCount aSuppressedSignalMax)
: m_SuppressedSignalMax(aSuppressedSignalMax) { }

//**********************************************************************

bool AdcCodeHelper::hasStickyBits(AdcCount sig) {
  const unsigned int onemask = 0x003f;
  unsigned int lsb= sig & onemask;
  bool allzero = lsb == 0;
  bool allone = lsb == onemask;
  return allzero || allone;
}

//**********************************************************************

bool AdcCodeHelper::isSmall(AdcCount sig, AdcSignal ped) {
  if ( m_SuppressedSignalMax <= 0.0 ) return false;
  AdcSignal pedsignal = sig - ped;
  return std::abs(pedsignal) < m_SuppressedSignalMax;
}

//**********************************************************************

AdcCount AdcCodeHelper::intSignal(AdcSignal fsig) const {
  if ( fsig > 0.1 ) return AdcCount(fsig + 0.5);
  if ( fsig < -0.1 ) return -AdcCount(-fsig + 0.5);
  return 0;
}

//**********************************************************************

AdcSignal AdcCodeHelper::subtract(AdcCount sig, AdcSignal ped) const {
  return sig - ped;
}

//**********************************************************************

AdcCount AdcCodeHelper::intSubtract(AdcCount sig, AdcSignal ped) const {
  return intSignal(subtract(sig, ped));
}

//**********************************************************************
