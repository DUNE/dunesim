cet_build_plugin(SpaceChargeServiceProtoDUNE art::service
              dunesim_SpaceCharge
              lardataalg::DetectorInfo
	      art::Framework_Principal
	      art::Persistency_Common 
	      art::Persistency_Provenance
	      ROOT::Core
	      art::Utilities canvas::canvas
	      messagefacility::MF_MessageLogger
  	      cetlib::cetlib cetlib_except::cetlib_except
)              

cet_build_plugin(SpaceChargeServiceProtoDUNEhd art::service
              dunesim_SpaceCharge
              lardataalg::DetectorInfo
              art::Framework_Core
              art::Framework_Principal
              art::Framework_Services_Registry
              art::Persistency_Common
              art::Persistency_Provenance
              art_root_io::tfile_support
              ROOT::Core
              art_root_io::TFileService_service
              art::Framework_Services_System_FileCatalogMetadata_service
              art::Utilities canvas::canvas
              messagefacility::MF_MessageLogger
              cetlib::cetlib cetlib_except::cetlib_except
)


cet_build_plugin(SpaceChargeServiceProtoDUNEvd art::service
              dunesim_SpaceCharge
              lardataalg::DetectorInfo
              art::Framework_Core
              art::Framework_Principal
              art::Framework_Services_Registry
              art::Persistency_Common
              art::Persistency_Provenance
              art_root_io::tfile_support
              ROOT::Core
              art_root_io::TFileService_service
              art::Framework_Services_System_FileCatalogMetadata_service
              art::Utilities
	      canvas::canvas
              messagefacility::MF_MessageLogger
              cetlib::cetlib
	      cetlib_except::cetlib_except
)
              
cet_build_plugin(SpaceChargeServiceProtoDUNEdp art::service
              dunesim_SpaceCharge
              lardataalg::DetectorInfo
              art::Framework_Principal
              art::Persistency_Common 
              art::Persistency_Provenance
              ROOT::Core
              art::Utilities canvas::canvas
              messagefacility::MF_MessageLogger
              cetlib::cetlib cetlib_except::cetlib_except

)


cet_build_plugin(SpaceChargeServiceDUNE35t art::service
              dunesim_SpaceCharge
	      art::Framework_Principal
	      art::Persistency_Common
	      art::Persistency_Provenance
	      ROOT::Core
	      art::Utilities
	      messagefacility::MF_MessageLogger
  	      cetlib::cetlib
              
)

cet_build_plugin(SpaceChargeService3x1x1dphase art::service
              dunesim_SpaceCharge
	      art::Framework_Principal
	      art::Persistency_Common
	      art::Persistency_Provenance
	      ROOT::Core
	      art::Utilities canvas::canvas
	      messagefacility::MF_MessageLogger
  	      cetlib::cetlib
              
)

install_headers()
install_fhicl()
install_source()
