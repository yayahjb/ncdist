#pragma once

#ifndef FOLLOWERTOOLS_H
#define FOLLOWERTOOLS_H

#include <algorithm>
#include <sstream>
#include <string>
#include <utility>

#include "ProjectorTools.h"
#include "triple.h"

std::pair<int,std::string>  IdentifyNearbyBoundaries(const double v[6], const double cutoff);

namespace FollowerTools {
   void OpenOutputFile( std::ofstream& svgOut, const std::string& sFileName );
};
#endif
