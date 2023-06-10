#ifndef BEAMMISCUTILS_H
#define BEAMMISCUTILS_H
#include <string>
#include "fhiclcpp/ParameterSetRegistry.h"
#include "fhiclcpp/ParameterSet.h"

namespace beammisc {
struct CountConfig {
  CountConfig(fhicl::ParameterSet & generator)
    : fNP04frontTreeName(generator.get<std::string>("NP04frontTreeName")),
      fTOF1TreeName(generator.get<std::string>("TOF1TreeName")),
      fBPROF1TreeName(generator.get<std::string>("BPROF1TreeName")),
      fBPROF2TreeName(generator.get<std::string>("BPROF2TreeName")),
      fBPROF3TreeName(generator.get<std::string>("BPROF3TreeName")),
      fTRIG1TreeName(generator.get<std::string>("TRIG1TreeName")),
      fBPROFEXTTreeName(generator.get<std::string>("BPROFEXTTreeName", "")),
      fBPROF4TreeName(generator.get<std::string>("BPROF4TreeName")),
      fTRIG2TreeName(generator.get<std::string>("TRIG2TreeName")) {};

  std::string fNP04frontTreeName;
  std::string fTOF1TreeName;
  std::string fBPROF1TreeName;
  std::string fBPROF2TreeName;
  std::string fBPROF3TreeName;
  std::string fTRIG1TreeName;
  std::string fBPROFEXTTreeName;
  std::string fBPROF4TreeName;
  std::string fTRIG2TreeName;
};

fhicl::ParameterSet GetGenerator(fhicl::ParameterSet & pset);
fhicl::ParameterSet GetPars(std::string fcl_file);

}

#endif
