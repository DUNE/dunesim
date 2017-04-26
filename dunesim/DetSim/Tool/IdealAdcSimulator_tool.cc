// IdealAdcSimulator_tool.cc

#include "IdealAdcSimulator.h"

//**********************************************************************

IdealAdcSimulator::IdealAdcSimulator(const fhicl::ParameterSet&) { }

//**********************************************************************

AdcSimulator::Count
IdealAdcSimulator::count(double vin, Channel, Tick) const {
  return 1234;
}

//**********************************************************************
