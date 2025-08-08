#include "GapInfo.hh"
#include <unordered_map>

//to make couts easier
std::string GapInfo::GapToString(GapType type)
{
  static const std::unordered_map<GapType, std::string> gapTypeToStringMap = {
    {GapType::X, "X"},
    {GapType::Y, "Y"},
    {GapType::Z, "Z"},
    {GapType::NotAGap, "NotAGap"},
    {GapType::NotInDetectorVolume, "NotInDetectorVolume"}};

  return gapTypeToStringMap.at(type);
}

std::string GapInfo::GapToString(MovingDirection direction)
{
  static const std::unordered_map<MovingDirection, std::string> movingDirectionToStringMap = {
    {MovingDirection::LEFT, "LEFT"}, {MovingDirection::RIGHT, "RIGHT"}};

  return movingDirectionToStringMap.at(direction);
}
