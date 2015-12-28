# dunezs/DZSService/CMakeLists.txt
#
# David Adams
# November 2015
#
# CMake instructions for dunezs services.

# Build the library excluding services.
#art_make(NO_PLUGINS
#         LIBRARY_NAME DZSService
#         LIB_LIBRARIES
##          ZeroSuppression35t_service
#           DZSCore
#        )

#art_make( BASENAME_ONLY
#          SERVICE_LIBRARIES
#            DZSCore
#            Geometry_service
#            ${FHICLCPP}
#)

# Add plugin for each service.
#simple_plugin(SimChannelExtract35tService   "service"
#                DZSCore
#                Geometry Simulation Utilities
#                SignalShapingServiceDUNE35t_service
#                ${ART_UTILITIES}
#                ${ROOT_BASIC_LIB_LIST} BASENAME_ONLY
#             )

#simple_plugin(ZeroSuppress35tService        "service" BASENAME_ONLY
#                DZSCore
#)
#simple_plugin(ZeroSuppress35tLegacyService  "service" BASENAME_ONLY
#                DZSCore
#)
#simple_plugin(ZeroSuppressFixedService      "service" BASENAME_ONLY
#)
#simple_plugin(CompressReplaceService        "service" BASENAME_ONLY
#)
#simple_plugin(SimChannelExtractAllService   "service"
#                DZSCore
#                Simulation
#                Utilities
#                SignalShapingServiceDUNE35t_service
#                ${ART_UTILITIES}
#                ${ROOT_BASIC_LIB_LIST} BASENAME_ONLY
#             )
#simple_plugin(ExponentialChannelNoiseService   "service"
#                Geometry
#                SignalShapingServiceDUNE35t_service
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
#                ${SEEDSERVICE_SERVICE}
#                ${ART_FRAMEWORK_CORE}
#                ${ART_UTILITIES}
#                ${CLHEP}
#                ${ROOT_BASIC_LIB_LIST}
#             )
#simple_plugin(WhiteChannelNoiseService   "service"
#                Geometry
#                SignalShapingServiceDUNE35t_service
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
#                ${SEEDSERVICE_SERVICE}
#                ${ART_FRAMEWORK_CORE}
#                ${ART_UTILITIES}
#                ${CLHEP}
#                ${ROOT_BASIC_LIB_LIST}
#             )
#simple_plugin(ProvidedPedestalAdditionService "service" BASENAME_ONLY
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
#                ${SEEDSERVICE_SERVICE}
#                ${ART_FRAMEWORK_CORE}
#                ${ART_UTILITIES}
#                ${CLHEP}
#                ${ROOT_BASIC_LIB_LIST}
#             )
#simple_plugin(StuckBitsService   "service"
#                Geometry
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
#                ${SEEDSERVICE_SERVICE}
#                ${ART_FRAMEWORK_CORE}
#                ${ART_UTILITIES}
#                ${MF_MESSAGELOGGER}
#                ${MF_UTILITIES}
#                ${CLHEP}
#                ${ROOT_BASIC_LIB_LIST}
#                BASENAME_ONLY
#             )
