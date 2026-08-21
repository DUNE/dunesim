////////////////////////////////////////////////////////////////////////
// \file ElectronDiverterDistortion_tool.cc
//
// Registers spacecharge::ElectronDiverterDistortion as an art tool so it can
// be loaded as a link of a ChainedDistorter via
// art::make_tool<detinfo::IDistortion>.
////////////////////////////////////////////////////////////////////////

#include "dunesim/SpaceCharge/ElectronDiverterDistortion.h"

#include "art/Utilities/ToolMacros.h"

DEFINE_ART_CLASS_TOOL(spacecharge::ElectronDiverterDistortion)
