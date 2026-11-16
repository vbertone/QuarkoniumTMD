//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "quarkoniumtmd/matchingfunctionspdfquarkonium.h"
#include "quarkoniumtmd/hardfactorquarkonium.h"

#include <apfel/constants.h>
#include <apfel/messages.h>

namespace apfel
{
  //_________________________________________________________________________________
  C1ggpdfQuarkonium::C1ggpdfQuarkonium(QuarkoniumSpecies const& Species):
    Expression()
  {
    // Only check that the species is known (throws otherwise). The
    // NLO virtual constant B no longer enters here.
    BQuarkonium(Species);
  }
  double C1ggpdfQuarkonium::Local(double const&) const
  {
    // Universal C_gg^(1) (coefficient of alpha_s / 4 pi). The NLO
    // virtual constant B now enters the hard factor
    // (see HardFactorQuarkonium in hardfactorquarkonium.h).
    return - CA * zeta2;
  }
}
