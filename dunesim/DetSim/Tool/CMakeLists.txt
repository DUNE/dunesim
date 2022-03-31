# dunetpc/dune/DetSim/Tool/CMakeLists.txt

# David Adams
# April 2017

# Add plugin for each tool.

cet_build_plugin(IdealAdcSimulator   
  art::tool
  art::Utilities 
  canvas::canvas
  cetlib::cetlib 
  cetlib_except::cetlib_except
  ROOT_BASIC_LIB_LIST
  )

cet_build_plugin(NoisyAdcSimulator   
  art::tool
  nurandom::RandomUtils_NuRandomService_service
  CLHEP
  art::Utilities 
  canvas::canvas
  cetlib::cetlib cetlib_except
  messagefacility::MF_MessageLogger
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core
  art_root_io::TFileService_service
  art::Framework_Core
  ROOT_BASIC_LIB_LIST

  )

add_subdirectory(test)
