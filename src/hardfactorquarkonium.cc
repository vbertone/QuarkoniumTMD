//
// QuarkoniumTMD
//

#include "quarkoniumtmd/hardfactorquarkonium.h"

#include <apfel/constants.h>
#include <apfel/messages.h>

#include <cmath>
#include <stdexcept>

namespace apfel
{
  //_________________________________________________________________________________
  double BQuarkonium(QuarkoniumSpecies const& Species)
  {
    // Reference: https://arxiv.org/pdf/hep-ph/9707223v2
    if (Species == Q_1S0_8 || Species == Q_approx)
      return CF * ( - 5 + Pi2 / 4 ) + CA * ( 2 + 7 * Pi2 / 24 ); // Eq. (129)
    else if (Species == Q_3P0_8)
      return CF * ( - 7. / 3 + Pi2 / 4 ) + CA * ( 17. / 54 + 35. / 27 * log(2) + 7 * Pi2 / 24 ); // Eq. (130)
    else if (Species == Q_3P2_8)
      return - 4 * CF + CA * ( 23. / 36 + 7 * log(2) / 9 + 5 * Pi2 / 12 ); // Eq. (131)
    else if (Species == Q_3P2_1)
      return - 4 * CF + CA * ( 1. / 3 + 5 * log(2) / 3. + Pi2 / 6 ); // Eq. (127)
    else
      throw std::runtime_error(error("BQuarkonium", "Unknown quarkonium species"));
  }

  //_________________________________________________________________________________
  double HardFactorQuarkonium(double const& alphas, QuarkoniumSpecies const& Species)
  {
    return 1 + alphas / M_PI * ( BQuarkonium(Species) - Pi2 * CA / 12 );
  }
}
