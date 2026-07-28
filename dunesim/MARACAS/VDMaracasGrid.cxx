////////////////////////////////////////////////////////////////////////
// \file VDMaracasGrid.cxx
////////////////////////////////////////////////////////////////////////

#include "dunesim/MARACAS/VDMaracasGrid.h"

#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"

#include "hdf5.h"

#include <algorithm> // std::min
#include <cmath>     // std::floor
#include <utility>   // std::move


namespace maracas {

  // Partial compound reads: describe a memory type holding just the one named
  // member (at offset 0) and let HDF5 extract it from the file's compound row.
  double readScalarD(hid_t dset, char const* name)
  {
    hid_t const mt = H5Tcreate(H5T_COMPOUND, sizeof(double));
    H5Tinsert(mt, name, 0, H5T_NATIVE_DOUBLE);
    double v = 0.;
    herr_t const st = H5Dread(dset, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
    H5Tclose(mt);
    if (st < 0)
      throw cet::exception("VDMaracasGrid") << "failed reading scalar '" << name << "'\n";
    return v;
  }

  unsigned readScalarU(hid_t dset, char const* name)
  {
    hid_t const mt = H5Tcreate(H5T_COMPOUND, sizeof(unsigned));
    H5Tinsert(mt, name, 0, H5T_NATIVE_UINT);
    unsigned v = 0;
    herr_t const st = H5Dread(dset, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
    H5Tclose(mt);
    if (st < 0)
      throw cet::exception("VDMaracasGrid") << "failed reading scalar '" << name << "'\n";
    return v;
  }

  std::vector<double>
  readArray(hid_t dset, char const* name, hsize_t d0, hsize_t d1, hsize_t d2)
  {
    hsize_t const dims[3] = {d0, d1, d2};
    hid_t const arrT = H5Tarray_create2(H5T_NATIVE_DOUBLE, 3, dims);
    hsize_t const n = d0 * d1 * d2;
    hid_t const mt = H5Tcreate(H5T_COMPOUND, sizeof(double) * n);
    H5Tinsert(mt, name, 0, arrT);
    std::vector<double> buf(n);
    herr_t const st = H5Dread(dset, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    H5Tclose(mt);
    H5Tclose(arrT);
    if (st < 0)
      throw cet::exception("VDMaracasGrid") << "failed reading array member '" << name << "'\n";
    return buf;
  }

  std::string resolvePath(std::string const& name)
  {
    cet::search_path sp{"FW_SEARCH_PATH"};
    std::string full;
    if (sp.find_file(name, full)) return full;
    return name; // fall back to a literal (e.g. absolute) path
  }




  //--------------------------------------------------------------------
  VDMaracasGrid::VDMaracasGrid(std::string const& filename,
                               std::vector<FieldSpec> const& fields)
  {
    std::string const path = resolvePath(filename);

    hid_t const fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
    hid_t const file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, fapl);
    H5Pclose(fapl);
    if (file < 0)
      throw cet::exception("VDMaracasGrid") << "cannot open HDF5 file '" << path << "'\n";

    // Grid header (metres / cell counts).
    hid_t const par = H5Dopen2(file, "parameters", H5P_DEFAULT);
    if (par < 0) {
      H5Fclose(file);
      throw cet::exception("VDMaracasGrid") << "no /parameters in '" << path << "'\n";
    }
    double const gmin[3] = {readScalarD(par, "xmin"), readScalarD(par, "ymin"),
                            readScalarD(par, "zmin")};
    double const gd[3] = {readScalarD(par, "dx"), readScalarD(par, "dy"), readScalarD(par, "dz")};
    unsigned const gN[3] = {readScalarU(par, "Nx"), readScalarU(par, "Ny"), readScalarU(par, "Nz")};
    H5Dclose(par);

    for (auto const& fs : fields) {
      hid_t const dset = H5Dopen2(file, fs.dataset.c_str(), H5P_DEFAULT);
      if (dset < 0) {
        H5Fclose(file);
        throw cet::exception("VDMaracasGrid")
          << "no /" << fs.dataset << " in '" << path << "'\n";
      }

      VecField vf;
      double const off = fs.nodeCentered ? 0.0 : 0.5;
      for (int a = 0; a < 3; ++a) {
        vf.n[a] = static_cast<int>(fs.nodeCentered ? gN[a] + 1 : gN[a]);
        vf.min[a] = gmin[a] + off * gd[a];
        vf.d[a] = gd[a];
      }
      for (int c = 0; c < 3; ++c)
        vf.v[c] = readArray(dset,
                            fs.members[c].c_str(),
                            static_cast<hsize_t>(vf.n[0]),
                            static_cast<hsize_t>(vf.n[1]),
                            static_cast<hsize_t>(vf.n[2]));

      H5Dclose(dset);
      fFields.emplace(fs.key, std::move(vf));
    }

    H5Fclose(file);
  }

  //--------------------------------------------------------------------
  std::array<double, 3>
  VDMaracasGrid::interpolate(std::string const& key, double x, double y, double z) const
  {
    auto const it = fFields.find(key);
    if (it == fFields.end())
      throw cet::exception("VDMaracasGrid") << "unknown field '" << key << "'\n";
    VecField const& f = it->second;

    double const p[3] = {x, y, z};
    int lo[3];
    double t[3];
    for (int a = 0; a < 3; ++a) {
      int const nmax = f.n[a] - 1; // largest valid index
      double fi = (p[a] - f.min[a]) / f.d[a];
      if (fi < 0.) fi = 0.;
      if (fi > nmax) fi = nmax;
      int l = static_cast<int>(std::floor(fi));
      if (l > nmax - 1) l = (nmax >= 1) ? nmax - 1 : 0; // keep l+1 in range
      lo[a] = l;
      t[a] = fi - l;
    }

    std::array<double, 3> out{0., 0., 0.};
    for (int c = 0; c < 3; ++c) {
      double res = 0.;
      for (int di = 0; di < 2; ++di)
        for (int dj = 0; dj < 2; ++dj)
          for (int dk = 0; dk < 2; ++dk) {
            int const i = std::min(lo[0] + di, f.n[0] - 1);
            int const j = std::min(lo[1] + dj, f.n[1] - 1);
            int const k = std::min(lo[2] + dk, f.n[2] - 1);
            double const w = (di ? t[0] : 1. - t[0]) * (dj ? t[1] : 1. - t[1]) *
                             (dk ? t[2] : 1. - t[2]);
            res += w * f.v[c][(i * f.n[1] + j) * f.n[2] + k];
          }
      out[c] = res;
    }
    return out;
  }

} // namespace maracas
