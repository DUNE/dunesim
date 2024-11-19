art_make(
          MODULE_LIBRARIES
                        larcoreobj::SummaryData
                        larcore::Geometry_Geometry_service
                        larcorealg::Geometry
                        larsim::Simulation lardataobj::Simulation
                        nusimdata::SimulationBase
                        ROOT::Core
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except::cetlib_except
			CLHEP::CLHEP
                        ROOT::Geom
        )

install_headers()
install_fhicl()
install_source()
