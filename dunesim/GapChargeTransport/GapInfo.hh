#ifndef GAPINFO_HH
#define GAPINFO_HH

#include <limits>
#include <string>

//all the info about the gap
class GapInfo {
public:
  enum class GapType { X, Y, Z, NotAGap, NotInDetectorVolume };
  enum class MovingDirection { LEFT, RIGHT };

  // Constructor
  GapInfo(GapType type = GapType::NotInDetectorVolume,
          double size = 0,
          double left = std::numeric_limits<double>::quiet_NaN(),
          double right = std::numeric_limits<double>::quiet_NaN(),
          MovingDirection mdir = MovingDirection::LEFT)
    : fGapType(type), fGapSize(size), fLeftBorder(left), fRightBorder(right), fMovingDirection(mdir)
  {}

  // Getters
  GapType GetType() const { return fGapType; }
  double GetSize() const { return fGapSize; }
  double GetLeftBorder() const { return fLeftBorder; }
  double GetRightBorder() const { return fRightBorder; }
  MovingDirection GetDirection() const { return fMovingDirection; }
  //Methods to quickly print the type and direction
  static std::string GapToString(GapType type);
  static std::string GapToString(MovingDirection direction);

private:
  GapType fGapType;
  double fGapSize;
  double fLeftBorder;
  double fRightBorder;
  MovingDirection fMovingDirection;
};

#endif
