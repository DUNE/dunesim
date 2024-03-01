art_make(
          MODULE_LIBRARIES
                        larcoreobj::SummaryData
                        larcore::Geometry_Geometry_service
			larcorealg::Geometry
                        nurandom::RandomUtils_NuRandomService_service
                        nugen::EventGeneratorBase_GENIE
                        larsim::Simulation lardataobj::Simulation
                        nusimdata::SimulationBase
                        ROOT::Core
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        art_root_io::tfile_support
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except::cetlib_except
			CLHEP::CLHEP
			SQLite::SQLite3
                        ROOT::EG
                        ifdh::ifdh
        )

install_headers()
install_fhicl()
install_source()
