////////////////////////////////////////////////////////////////////////
// \file CSUDistortionProtoDUNE_tool.cc
//
// Registers spacecharge::CSUDistortionProtoDUNE as an art tool so it can be
// loaded as a link of a ChainedDistorter via
// art::make_tool<detinfo::IDistortion>.
////////////////////////////////////////////////////////////////////////

#include "dunesim/SpaceCharge/CSUDistortionProtoDUNE.h"

#include "art/Utilities/ToolMacros.h"

DEFINE_ART_CLASS_TOOL(spacecharge::CSUDistortionProtoDUNE)
