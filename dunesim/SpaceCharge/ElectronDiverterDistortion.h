////////////////////////////////////////////////////////////////////////
// \file ElectronDiverterDistortion.h
//
// \brief IDistortion link modeling an electron-diverter position distortion.
//
// A chain link (detinfo::IDistortion) for ChainedDistorter, registered as an
// art tool (see ElectronDiverterDistortion_tool.cc). Each link decides
// internally whether to apply its transform and/or forward to the next link
// via CallNextDistort()/CallNextCorrect().
////////////////////////////////////////////////////////////////////////
#ifndef DUNESIM_SPACECHARGE_ELECTRONDIVERTERDISTORTION_H
#define DUNESIM_SPACECHARGE_ELECTRONDIVERTERDISTORTION_H

#include "lardataalg/DetectorInfo/IDistortion.h"

#include "fhiclcpp/types/Atom.h"

namespace fhicl {
  class ParameterSet;
}

namespace spacecharge {

  class ElectronDiverterDistortion : public detinfo::IDistortion {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<double> ZCenter{
          Name("ZCenter"),
          Comment("Location of this electron diverter region.")
      };
      fhicl::Atom<double> AXPosOffs{
          Name("AXPosOffs"),
          Comment("X distortion shift scale in cm.")
      };
      fhicl::Atom<double> BZPosOffs{
          Name("BZPosOffs"),
          Comment("Z distortion shift scale in cm.")
      };
      fhicl::Atom<double> Width{
          Name("Width"),
          Comment("Width of distortion function in cm.")
      };
      fhicl::Atom<double> ChargeLossZLow{
          Name("ChargeLossZLow"),
          Comment("Lower bound of this electron diverter region.")
      };
      fhicl::Atom<double> ChargeLossZHigh{
          Name("ChargeLossZHigh"),
          Comment("Upper bound of this electron diverter region.")
      };
    };

    explicit ElectronDiverterDistortion(fhicl::ParameterSet const& pset);

    geo::Point_t Distort(geo::Point_t const& point) const override;
    geo::Point_t Correct(geo::Point_t const& point) const override;
  private:
    double fZCenter;
    double fAXPosOffs;
    double fBZPosOffs;
    double fWidth;
    double fChargeLossZLow;
    double fChargeLossZHigh;
  }; // class ElectronDiverterDistortion

} // namespace spacecharge

#endif // DUNESIM_SPACECHARGE_ELECTRONDIVERTERDISTORTION_H
