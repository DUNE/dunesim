////////////////////////////////////////////////////////////////////////
// \file VDMaracasDistortion.cxx
////////////////////////////////////////////////////////////////////////

#include "dunesim/MARACAS/VDMaracasDistortion.h"

#include "lardataalg/DetectorInfo/IgnorableToolConfigKeys.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"

#include <vector>

namespace maracas {

  //--------------------------------------------------------------------
  VDMaracasDistortion::VDMaracasDistortion(fhicl::ParameterSet const& pset)
  {
    fhicl::Table<Config> const config{pset, detinfo::IgnorableToolConfigKeys()};

    auto const o = config().OriginCm();
    fOriginCm = {o[0], o[1], o[2]};

    std::vector<VDMaracasGrid::FieldSpec> const dist{
      {"fwd", "dist_3D", {"forward_delta_x", "forward_delta_y", "forward_delta_z"}, false},
      {"bwd", "dist_3D", {"backward_delta_x", "backward_delta_y", "backward_delta_z"}, false}};

    fTop = std::make_unique<VDMaracasGrid>(config().TopFile(), dist);
    fBot = std::make_unique<VDMaracasGrid>(config().BotFile(), dist);
  }

  //--------------------------------------------------------------------
  geo::Point_t VDMaracasDistortion::apply(geo::Point_t const& point, bool const forward) const
  {
    // LArSoft cm -> MARACAS m (x is drift; y,z identical; no flips).
    double const xm = (point.X() - fOriginCm[0]) / 100.;
    double const ym = (point.Y() - fOriginCm[1]) / 100.;
    double const zm = (point.Z() - fOriginCm[2]) / 100.;

    VDMaracasGrid const& g = (xm >= 0.) ? *fTop : *fBot; // top: drift x>=0, bot: x<0
    auto const d = g.interpolate(forward ? "fwd" : "bwd", xm, ym, zm); // metres

    // Offsets metres -> cm, added to the input position.
    return geo::Point_t{point.X() + d[0] * 100., point.Y() + d[1] * 100., point.Z() + d[2] * 100.};
  }

  //--------------------------------------------------------------------
  geo::Point_t VDMaracasDistortion::Distort(geo::Point_t const& point) const
  {
    return apply(point, true);
  }

  //--------------------------------------------------------------------
  geo::Point_t VDMaracasDistortion::Correct(geo::Point_t const& point) const
  {
    return apply(point, false);
  }

} // namespace maracas
