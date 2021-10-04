# dunetpc/dune/DetSim/Tool/CMakeLists.txt

# David Adams
# April 2017

# Add plugin for each tool.

cet_build_plugin(IdealAdcSimulator   art::tool
  art_Utilities canvas
  cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  )

cet_build_plugin(NoisyAdcSimulator   art::tool
  nurandom_RandomUtils_NuRandomService_service
  CLHEP
  art_Utilities canvas
  cetlib cetlib_except
  messagefacility::MF_MessageLogger
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  ROOT_BASIC_LIB_LIST
  ART_FRAMEWORK_CORE
  )

add_subdirectory(test)
