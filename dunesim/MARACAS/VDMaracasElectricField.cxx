////////////////////////////////////////////////////////////////////////
// \file VDMaracasElectricField.cxx
////////////////////////////////////////////////////////////////////////

#include "dunesim/MARACAS/VDMaracasElectricField.h"

#include "larcorealg/CoreUtils/ProviderUtil.h" // lar::IgnorableProviderConfigKeys()

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"

#include <vector>

namespace maracas {

  //--------------------------------------------------------------------
  VDMaracasElectricField::VDMaracasElectricField(fhicl::ParameterSet const& pset)
  {
    fhicl::Table<Config> const config{pset, lar::IgnorableProviderConfigKeys()};

    auto const o = config().OriginCm();
    fOriginCm = {o[0], o[1], o[2]};
    fScale = config().FieldScaleToKVCm();

    std::vector<VDMaracasGrid::FieldSpec> const efield{
      {"E", "SCE_3D", {"Ex", "Ey", "Ez"}, true}}; // node-centred

    fTop = std::make_unique<VDMaracasGrid>(config().TopFile(), efield);
    fBot = std::make_unique<VDMaracasGrid>(config().BotFile(), efield);
  }

  //--------------------------------------------------------------------
  geo::Vector_t VDMaracasElectricField::Efield(geo::Point_t const& point) const
  {
    // LArSoft cm -> MARACAS m (x is drift; y,z identical; no flips).
    double const xm = (point.X() - fOriginCm[0]) / 100.;
    double const ym = (point.Y() - fOriginCm[1]) / 100.;
    double const zm = (point.Z() - fOriginCm[2]) / 100.;

    VDMaracasGrid const& g = (xm >= 0.) ? *fTop : *fBot; // top: drift x>=0, bot: x<0
    auto const e = g.interpolate("E", xm, ym, zm);
    return geo::Vector_t(e[0] * fScale, e[1] * fScale, e[2] * fScale);
  }

} // namespace maracas
