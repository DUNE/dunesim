////////////////////////////////////////////////////////////////////////
// \file VDMaracasDistortion_tool.cc
//
// Registers maracas::VDMaracasDistortion as an art tool so it can be loaded
// as a link of a ChainedDistorter via art::make_tool<detinfo::IDistortion>.
////////////////////////////////////////////////////////////////////////

#include "dunesim/MARACAS/VDMaracasDistortion.h"

#include "art/Utilities/ToolMacros.h"

DEFINE_ART_CLASS_TOOL(maracas::VDMaracasDistortion)
