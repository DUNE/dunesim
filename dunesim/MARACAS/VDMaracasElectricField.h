////////////////////////////////////////////////////////////////////////
// \file VDMaracasElectricField.h
//
// \brief IElectricFieldProvider backed by the MARACAS /SCE_3D E-field map.
////////////////////////////////////////////////////////////////////////
#ifndef MARACAS_VDMARACASELECTRICFIELD_H
#define MARACAS_VDMARACASELECTRICFIELD_H

#include "dunesim/MARACAS/VDMaracasGrid.h"
#include "lardataalg/DetectorInfo/IElectricFieldProvider.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include <array>
#include <memory>
#include <string>

namespace fhicl {
  class ParameterSet;
}

namespace maracas {

  class VDMaracasElectricField : public detinfo::IElectricFieldProvider {
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
      fhicl::Atom<double> FieldScaleToKVCm{
        Name("FieldScaleToKVCm"),
        Comment("factor converting the stored field to kV/cm (V/m -> 1e-5)"),
        1.0e-5};
    };

    explicit VDMaracasElectricField(fhicl::ParameterSet const& pset);

    geo::Vector_t Efield(geo::Point_t const& point) const override;

  private:
    std::array<double, 3> fOriginCm;
    double fScale;
    std::unique_ptr<VDMaracasGrid> fTop;
    std::unique_ptr<VDMaracasGrid> fBot;
  }; // class VDMaracasElectricField

} // namespace maracas

#endif // MARACAS_VDMARACASELECTRICFIELD_H
