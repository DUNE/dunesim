#include "BeamMiscUtils.h"
fhicl::ParameterSet beammisc::GetGenerator(fhicl::ParameterSet & pset) {
  return pset.get<fhicl::ParameterSet>("physics.producers.generator");
}

fhicl::ParameterSet beammisc::GetPars(std::string fcl_file) {
  std::cout << "Fcl file: " << fcl_file << std::endl;
  fhicl::ParameterSet pset;

  // Configuration file lookup policy.
  char const* fhicl_env = getenv("FHICL_FILE_PATH");
  std::string search_path;

  if (fhicl_env == nullptr) {
    std::cerr << "Expected environment variable FHICL_FILE_PATH is missing or empty: using \".\"\n";
    search_path = ".";
  }
  else {
    search_path = std::string{fhicl_env};
  }

  cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};

  //fhicl::make_ParameterSet(fcl_file, lookupPolicy, pset);
  pset = fhicl::ParameterSet::make(fcl_file, lookupPolicy);
  return pset;

}
