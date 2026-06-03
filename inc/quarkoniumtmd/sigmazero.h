//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "quarkoniumtmd/constantsquarkonium.h"
#include <cmath>

namespace apfel
{
  double sigmazero(double const& alphas, double const& M, std::string const& fit, QuarkoniumSpecies const& Species);
}
