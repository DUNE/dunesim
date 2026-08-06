# dunetpc/dune/larsim_DetSim/Service/test/CMakeLists.txt

include(CetTest)

cet_enable_asserts()

cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR})
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})
cet_test_env_prepend(FW_SEARCH_PATH ${PROJECT_BINARY_DIR}/../dunecore/gdml)
cet_test_env_prepend(ROOT_LIBRARY_PATH ${PROJECT_BINARY_DIR}/lib)

MESSAGE( STATUS "CMAKE_LIBRARY_PATH: " ${CMAKE_LIBRARY_PATH} )

cet_test(test_LarsoftHuffmanCompressService SOURCE test_LarsoftHuffmanCompressService.cxx
  LIBRARIES
    dunesim::DetSim_Service_LarsoftHuffmanCompressService_service
    lardataobj::RawData
)

cet_test(test_FixedDetPedestalService SOURCE test_FixedDetPedestalService.cxx
  LIBRARIES
    dunecore::ArtSupport
    larcorealg::Geometry
    dunesim::DetSim_Service_FixedDetPedestalService_service
)

cet_test(test_ExponentialChannelNoiseService SOURCE test_ExponentialChannelNoiseService.cxx
  LIBRARIES
    dunesim::DetSim_Service_ExponentialChannelNoiseService_service
    dunecore::ArtSupport
    larcorealg::Geometry
)
