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
    : gapType(type), gapSize(size), leftBorder(left), rightBorder(right), movingDirection(mdir)
  {}

  // Getters
  GapType GetType() const { return gapType; }
  double GetSize() const { return gapSize; }
  double GetLeftBorder() const { return leftBorder; }
  double GetRightBorder() const { return rightBorder; }
  MovingDirection GetDirection() const { return movingDirection; }
  //Methods to quickly print the type and direction
  static std::string GapToString(GapType type);
  static std::string GapToString(MovingDirection direction);

private:
  GapType gapType;
  double gapSize;
  double leftBorder;
  double rightBorder;
  MovingDirection movingDirection;
};

#endif
