//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "quarkoniumtmd/constantsquarkonium.h"

#include <apfel/tmdbuilder.h>

namespace apfel
{
  std::map<int, TmdObjects> InitializeTmdObjectsQuarkonium(Grid                const& g,
                                                           std::vector<double> const& Thresholds,
                                                           QuarkoniumSpecies   const& Species,
                                                           double              const& IntEps = 1e-5);
}
