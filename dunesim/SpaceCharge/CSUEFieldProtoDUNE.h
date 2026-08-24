////////////////////////////////////////////////////////////////////////
// \file CSUEFieldProtoDUNE.h
//
// \brief IElectricFieldProvider backed by the CSU ProtoDUNE space-charge
//        maps (CSUSpaceChargeMaps). Turns the fractional dE/E offsets into
//        an absolute field vector in kV/cm using a configured nominal field.
////////////////////////////////////////////////////////////////////////
#ifndef DUNESIM_SPACECHARGE_CSUEFIELDPROTODUNE_H
#define DUNESIM_SPACECHARGE_CSUEFIELDPROTODUNE_H

#include "dunesim/SpaceCharge/CSUSpaceChargeMaps.h"
#include "lardataalg/DetectorInfo/IElectricFieldProvider.h"

namespace fhicl {
  class ParameterSet;
}

namespace spacecharge {

  class CSUEFieldProtoDUNE : public detinfo::IElectricFieldProvider {
  public:
    explicit CSUEFieldProtoDUNE(fhicl::ParameterSet const& pset);

    /// Absolute electric field at \p point, in kV/cm.
    geo::Vector_t Efield(geo::Point_t const& point) const override;

  private:
    CSUSpaceChargeMaps fMaps;
    double fE0; ///< nominal drift field [kV/cm]
  }; // class CSUEFieldProtoDUNE

} // namespace spacecharge

#endif // DUNESIM_SPACECHARGE_CSUEFIELDPROTODUNE_H
