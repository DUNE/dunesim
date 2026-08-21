////////////////////////////////////////////////////////////////////////
// \file VDMaracasDistortion.h
//
// \brief IDistortion link backed by the MARACAS /dist_3D offset maps.
//
// A first-match chain link (detinfo::IDistortion) for ChainedDistorter,
// registered as an art tool (see VDMaracasDistortion_tool.cc).
////////////////////////////////////////////////////////////////////////
#ifndef MARACAS_VDMARACASDISTORTION_H
#define MARACAS_VDMARACASDISTORTION_H

#include "dunesim/MARACAS/VDMaracasGrid.h"
#include "lardataalg/DetectorInfo/IDistortion.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include <array>
#include <memory>
#include <string>

namespace fhicl {
  class ParameterSet;
}

namespace maracas {

  class VDMaracasDistortion : public detinfo::IDistortion {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<std::string> TopFile{Name("TopFile"),
                                       Comment("PyTables HDF5 map for the drift x>=0 volume")};
      fhicl::Atom<std::string> BotFile{Name("BotFile"),
                                       Comment("PyTables HDF5 map for the drift x<0 volume")};
      fhicl::Sequence<double, 3> OriginCm{
        Name("OriginCm"),
        Comment("LArSoft-cm position of the MARACAS origin {x,y,z}")};
    };

    explicit VDMaracasDistortion(fhicl::ParameterSet const& pset);

    geo::Point_t Distort(geo::Point_t const& point) const override; ///< forward_delta
    geo::Point_t Correct(geo::Point_t const& point) const override; ///< backward_delta

  private:
    geo::Point_t apply(geo::Point_t const& point, bool forward) const;

    std::array<double, 3> fOriginCm;
    std::unique_ptr<VDMaracasGrid> fTop;
    std::unique_ptr<VDMaracasGrid> fBot;
  }; // class VDMaracasDistortion

} // namespace maracas

#endif // MARACAS_VDMARACASDISTORTION_H
