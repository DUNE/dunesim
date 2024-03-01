cet_build_plugin(GenExtCounterFilter art::module
		      lardataobj::RawData
		      lardataobj::RecoBase
		      dunecore::Geometry
		        larcorealg::Geometry
			larcore::Geometry_Geometry_service
			lardata::Utilities
			nusimdata::SimulationBase
                        ROOT::Core
 			art::Persistency_Common 
			art::Persistency_Provenance
			art::Utilities canvas::canvas
			messagefacility::MF_MessageLogger
			cetlib::cetlib cetlib_except::cetlib_except
              BASENAME_ONLY
)

cet_build_plugin(LArG4ExtCounterFilter art::module
		      lardataobj::RawData
		      lardataobj::RecoBase
		      dunecore::Geometry
		        larcorealg::Geometry
			larcore::Geometry_Geometry_service
			lardata::Utilities
			nusimdata::SimulationBase
                        ROOT::Core
 			art::Persistency_Common canvas::canvas
			art::Persistency_Provenance
			art::Utilities
			messagefacility::MF_MessageLogger
			cetlib::cetlib cetlib_except::cetlib_except
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
