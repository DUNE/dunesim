# dunetpc/dune/larsim_DetSim/Service/test/CMakeLists.txt

include(CetTest)

cet_enable_asserts()

MESSAGE( STATUS "CMAKE_LIBRARY_PATH: " ${CMAKE_LIBRARY_PATH} )

cet_test(test_LarsoftHuffmanCompressService SOURCES test_LarsoftHuffmanCompressService.cxx
  LIBRARIES
    dune_DetSim_Service_LarsoftHuffmanCompressService_service
    lardataobj_RawData
)

cet_test(test_FixedDetPedestalService SOURCES test_FixedDetPedestalService.cxx
  LIBRARIES
    dune_ArtSupport
    larcorealg_Geometry
)

cet_test(test_ExponentialChannelNoiseService SOURCES test_ExponentialChannelNoiseService.cxx
  LIBRARIES
    dune_DetSim_Service_ExponentialChannelNoiseService_service
    dune_ArtSupport
    larcorealg_Geometry
)
