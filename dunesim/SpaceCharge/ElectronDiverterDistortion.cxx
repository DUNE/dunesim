////////////////////////////////////////////////////////////////////////
// \file ElectronDiverterDistortion.cxx
////////////////////////////////////////////////////////////////////////

#include "dunesim/SpaceCharge/ElectronDiverterDistortion.h"

#include "lardataalg/DetectorInfo/IgnorableToolConfigKeys.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "TMath.h"

namespace spacecharge {

  //--------------------------------------------------------------------
  ElectronDiverterDistortion::ElectronDiverterDistortion(fhicl::ParameterSet const& pset)
  {
    // Ignore framework + art::make_tool "tool_type" keys during validation.
    fhicl::Table<Config> const config{pset, detinfo::IgnorableToolConfigKeys()};
    
    fZCenter         = config().ZCenter();
    fAXPosOffs       = config().AXPosOffs();
    fBZPosOffs       = config().BZPosOffs();
    fWidth           = config().Width();
    fChargeLossZLow  = config().ChargeLossZLow();
    fChargeLossZHigh = config().ChargeLossZHigh();
  }

  //--------------------------------------------------------------------
  geo::Point_t ElectronDiverterDistortion::Distort(geo::Point_t const& point) const {

    //Diverters only on the beam side
    if (point.X() > 0) return CallNextDistort(point);

    double z = point.Z();

    //'Catastrophic' region
    //Make a note that this is won't pass it along
    if (z > fChargeLossZLow && z < fChargeLossZHigh) {
      return geo::Point_t{point.X() - 2.e9, point.Y() + 2.e9, point.Z() + 2.e9};
    }

    
    double zdiff = z - fZCenter;
    double zexp = TMath::Exp( -TMath::Sq(zdiff/fWidth));
    double zoffset = (fBZPosOffs * zdiff * zexp);

    // the timing offsets need to be computed after the z shift
    double zdiffc = zdiff + zoffset;
    double zexpc = TMath::Exp( -TMath::Sq(zdiffc/fWidth) );
    geo::Vector_t offset{
      -1.*(fAXPosOffs * zexpc),
      0.,
      zoffset
    };

    return CallNextDistort(point + offset);
  }

  //--------------------------------------------------------------------
  geo::Point_t ElectronDiverterDistortion::Correct(geo::Point_t const& point) const
  {
    // Offer no correction
    return CallNextCorrect(point);
  }

} // namespace spacecharge
