//
// Created by lurker on 3/23/17.
//

#ifndef LEVELSET_ELECTRIC_H
#define LEVELSET_ELECTRIC_H

#include "bbfmm.h"
#include "gmres.h"
#include "levelset.h"
#include "Config.h"

void electric(Grid& g, levelset& ls, Surface& surf, Molecule& mol, scalar_t rescale, Config& cfg);

#endif //LEVELSET_ELECTRIC_H
