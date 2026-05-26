//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "quarkoniumtmd/constantsquarkonium.h"

#include <apfel/expression.h>

namespace apfel
{  
  class C1ggpdfQuarkonium: public Expression
  {
  public:
    C1ggpdfQuarkonium(QuarkoniumSpecies const& Species);
    double Local(double const&) const;
  private:
    double _BQc;
  };
}
