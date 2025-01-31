art_make(
          MODULE_LIBRARIES
                        larcoreobj::SummaryData
                        larcore::Geometry_Geometry_service
			larcorealg::Geometry
                        nurandom::RandomUtils_NuRandomService_service
                        larsim::Simulation lardataobj::Simulation
                        nusimdata::SimulationBase
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except::cetlib_except
                        CLHEP::CLHEP
			ROOT::Core
			ROOT::Hist
			SQLite::SQLite3
                        ROOT::EG
        )

install_headers()
install_fhicl()
install_source()
