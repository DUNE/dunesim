
art_make(
          MODULE_LIBRARIES
                        larcoreobj_SummaryData
                        larcore_Geometry_Geometry_service
			                  larcorealg_Geometry
                        nurandom::RandomUtils_NuRandomService_service
                        larsim_Simulation lardataobj_Simulation
                        nusimdata::SimulationBase
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger

                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        CLHEP
                        ROOT_BASIC_LIB_LIST
                        SQLITE3
                        ROOT_EG
                        ifdh::ifdh
        )

install_headers()
install_fhicl()
install_source()
