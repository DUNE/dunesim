cet_build_plugin(GenExtCounterFilter art::module
		      lardataobj_RawData
		      lardataobj_RecoBase
		      dunecore_Geometry
		        larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			nusimdata::SimulationBase
                        ROOT::Core
 			art::Persistency_Common canvas
			art::Persistency_Provenance canvas
			art::Utilities canvas
			messagefacility::MF_MessageLogger
			cetlib::cetlib cetlib_except
              BASENAME_ONLY
)

cet_build_plugin(LArG4ExtCounterFilter art::module
		      lardataobj_RawData
		      lardataobj_RecoBase
		      dunecore_Geometry
		        larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			nusimdata::SimulationBase
                        ROOT::Core
 			art::Persistency_Common canvas
			art::Persistency_Provenance canvas
			art::Utilities canvas
			messagefacility::MF_MessageLogger
			cetlib::cetlib cetlib_except
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
