#include "tilebeam2013.h"
#include "tilebeam2014.h"
#include "tilebeam2016.h"
#include "lnaimpedance.h"
#include "tileimpedance.h"
#include "tilebeambase.h"

#ifndef TILE_BEAM_H

// Define the default implementation
// Other options:
// typedef TileBeamBase<TileBeam2013> TileBeam;
// typedef TileBeamBase<TileBeam2016> TileBeam;
typedef TileBeamBase<TileBeam2014> TileBeam;

#endif
