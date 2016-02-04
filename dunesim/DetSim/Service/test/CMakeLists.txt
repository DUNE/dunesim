# dunetpc/dune/DetSim/Service/test/CMakeLists.txt

cet_test(test_LarsoftHuffmanCompressService SOURCES test_LarsoftHuffmanCompressService.cxx
  LIBRARIES
    dune_DetSim_Service_LarsoftHuffmanCompressService_service
    RawData
)
