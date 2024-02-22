art_make(
          MODULE_LIBRARIES
                        larcoreobj_SummaryData
                        larcore_Geometry_Geometry_service
			larcorealg_Geometry
                        nurandom::RandomUtils_NuRandomService_service
                        larsim_Simulation lardataobj_Simulation
                        nusimdata::SimulationBase
                        ROOT::Core
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
			CLHEP::CLHEP
			SQLite::SQLite3
                        ROOT::EG
                        ifdh::ifdh
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
