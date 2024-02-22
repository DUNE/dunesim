cet_build_plugin(SpaceChargeServiceProtoDUNE art::service
              dunesim_SpaceCharge
              lardataalg_DetectorInfo
	      art::Framework_Principal
	      art::Persistency_Common canvas
	      art::Persistency_Provenance canvas
	      ROOT::Core
	      art::Utilities canvas
	      messagefacility::MF_MessageLogger
  	      cetlib::cetlib cetlib_except
)              
cet_build_plugin(SpaceChargeServiceProtoDUNEdp art::service
              dunesim_SpaceCharge
              lardataalg_DetectorInfo
	      art::Framework_Principal
              art::Persistency_Common canvas
              art::Persistency_Provenance canvas
              ROOT::Core
              art::Utilities canvas
              messagefacility::MF_MessageLogger
              cetlib::cetlib cetlib_except

)


cet_build_plugin(SpaceChargeServiceDUNE35t art::service
              dunesim_SpaceCharge
	      art::Framework_Principal
	      art::Persistency_Common canvas
	      art::Persistency_Provenance canvas
	      ROOT::Core
	      art::Utilities canvas
	      messagefacility::MF_MessageLogger
  	      cetlib::cetlib cetlib_except
              
)

cet_build_plugin(SpaceChargeService3x1x1dphase art::service
              dunesim_SpaceCharge
	      art::Framework_Principal
	      art::Persistency_Common canvas
	      art::Persistency_Provenance canvas
	      ROOT::Core
	      art::Utilities canvas
	      messagefacility::MF_MessageLogger
  	      cetlib::cetlib cetlib_except
              
)

install_headers()
install_fhicl()
install_source()
