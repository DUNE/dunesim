////////////////////////////////////////////////////////////////////////
// \file CSUDistortionProtoDUNE.h
//
// \brief IDistortion chain link backed by the CSU ProtoDUNE space-charge
//        spatial maps (CSUSpaceChargeMaps). Distort = sim (forward) offsets,
//        Correct = calibration (backward) offsets. Registered as an art tool
//        (see CSUDistortionProtoDUNE_tool.cc); compose with other links (e.g.
//        ElectronDiverterDistortion) in a ChainedDistorter.
////////////////////////////////////////////////////////////////////////
#ifndef DUNESIM_SPACECHARGE_CSUDISTORTIONPROTODUNE_H
#define DUNESIM_SPACECHARGE_CSUDISTORTIONPROTODUNE_H

#include "dunesim/SpaceCharge/CSUSpaceChargeMaps.h"
#include "lardataalg/DetectorInfo/IDistortion.h"

namespace fhicl {
  class ParameterSet;
}

namespace spacecharge {

  class CSUDistortionProtoDUNE : public detinfo::IDistortion {
  public:
    explicit CSUDistortionProtoDUNE(fhicl::ParameterSet const& pset);

    geo::Point_t Distort(geo::Point_t const& point) const override;
    geo::Point_t Correct(geo::Point_t const& point) const override;

  private:
    CSUSpaceChargeMaps fMaps;
  }; // class CSUDistortionProtoDUNE

} // namespace spacecharge

#endif // DUNESIM_SPACECHARGE_CSUDISTORTIONPROTODUNE_H
