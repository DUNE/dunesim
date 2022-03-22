# dunetpc/dune/larsim_DetSim/Service/test/CMakeLists.txt

include(CetTest)

cet_enable_asserts()

MESSAGE( STATUS "CMAKE_LIBRARY_PATH: " ${CMAKE_LIBRARY_PATH} )

cet_test(test_LarsoftHuffmanCompressService SOURCES test_LarsoftHuffmanCompressService.cxx
  LIBRARIES
    dunesim::DetSim_Service_LarsoftHuffmanCompressService_service
    lardataobj::RawData
)

cet_test(test_FixedDetPedestalService SOURCES test_FixedDetPedestalService.cxx
  LIBRARIES
    dunecore::ArtSupport
    larcorealg::Geometry
    dunesim::DetSim_Service_FixedDetPedestalService_service
)

cet_test(test_ExponentialChannelNoiseService SOURCES test_ExponentialChannelNoiseService.cxx
  LIBRARIES
    dunesim::DetSim_Service_ExponentialChannelNoiseService_service
    dunecore::ArtSupport
    larcorealg::Geometry
)
