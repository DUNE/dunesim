////////////////////////////////////////////////////////////////////////
// \file VDMaracasGrid.h
//
// \brief Reader/interpolator for the MARACAS ProtoDUNE-VD HDF5 maps.
//
// Opens a PyTables-structured HDF5 file read-only, reads the /parameters
// grid header, and loads one or more named 3-component vector fields (each
// three H5T_ARRAY members of a compound dataset such as /SCE_3D or /dist_3D),
// then closes the file. Provides clamped trilinear interpolation in the
// MARACAS frame (metres). Art-free; HDF5 usage is confined to the .cxx.
////////////////////////////////////////////////////////////////////////
#ifndef MARACAS_VDMARACASGRID_H
#define MARACAS_VDMARACASGRID_H

#include <array>
#include <map>
#include <string>
#include <vector>

namespace maracas {

  class VDMaracasGrid {
  public:
    /// One 3-component vector field to load from the file.
    struct FieldSpec {
      std::string key;                    ///< name to retrieve the field later
      std::string dataset;                ///< HDF5 dataset, e.g. "SCE_3D" / "dist_3D"
      std::array<std::string, 3> members; ///< the three component member names
      bool nodeCentered;                  ///< true: samples at min+i*d (N+1 per axis);
                                          ///< false: cell centres at min+(i+0.5)*d (N per axis)
    };

    /// Opens `filename` (resolved on FW_SEARCH_PATH, else used literally),
    /// reads /parameters, loads all `fields`, and closes the file.
    VDMaracasGrid(std::string const& filename, std::vector<FieldSpec> const& fields);

    /// Clamped trilinear interpolation of the named 3-vector field at MARACAS
    /// coordinates (x,y,z) in metres.
    std::array<double, 3> interpolate(std::string const& key,
                                      double x,
                                      double y,
                                      double z) const;

  private:
    struct VecField {
      std::array<int, 3> n;                     ///< samples per axis
      std::array<double, 3> min;                ///< coordinate of sample 0 [m]
      std::array<double, 3> d;                  ///< spacing [m]
      std::array<std::vector<double>, 3> v;     ///< components, flat row-major
    };

    std::map<std::string, VecField> fFields;
  }; // class VDMaracasGrid

} // namespace maracas

#endif // MARACAS_VDMARACASGRID_H
