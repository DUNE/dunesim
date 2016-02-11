# dunetpc/dune/DetSim/Service/test/CMakeLists.txt

include(CetTest)

MESSAGE( STATUS "CMAKE_LIBRARY_PATH: " ${CMAKE_LIBRARY_PATH} )

cet_test(test_LarsoftHuffmanCompressService SOURCES test_LarsoftHuffmanCompressService.cxx
  LIBRARIES
    dune_DetSim_Service_LarsoftHuffmanCompressService_service
    RawData
)

cet_test(test_ExponentialChannelNoiseService SOURCES test_ExponentialChannelNoiseService.cxx
  LIBRARIES
    dune_DetSim_Service_ExponentialChannelNoiseService_service
    dune_ArtSupport
)
