////////////////////////////////////////////////////////////////////////
// \file VDMaracasDistorter.h
//
// \brief IPositionDistorter backed by the MARACAS /dist_3D offset maps.
////////////////////////////////////////////////////////////////////////
#ifndef MARACAS_VDMARACASDISTORTER_H
#define MARACAS_VDMARACASDISTORTER_H

#include "dunesim/MARACAS/VDMaracasGrid.h"
#include "lardataalg/DetectorInfo/IPositionDistorter.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include <array>
#include <memory>
#include <string>

namespace fhicl {
  class ParameterSet;
}

namespace maracas {

  class VDMaracasDistorter : public detinfo::IPositionDistorter {
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

    explicit VDMaracasDistorter(fhicl::ParameterSet const& pset);

    geo::Point_t Distort(geo::Point_t const& point) const override; ///< forward_delta
    geo::Point_t Correct(geo::Point_t const& point) const override; ///< backward_delta

  private:
    geo::Point_t apply(geo::Point_t const& point, bool forward) const;

    std::array<double, 3> fOriginCm;
    std::unique_ptr<VDMaracasGrid> fTop;
    std::unique_ptr<VDMaracasGrid> fBot;
  }; // class VDMaracasDistorter

} // namespace maracas

#endif // MARACAS_VDMARACASDISTORTER_H
