//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "quarkoniumtmd/sigmazero.h"

#include <apfel/messages.h>

namespace apfel
{
  //_________________________________________________________________________
  double sigmazero(double const& alphas, double const& M, std::string const& fit, QuarkoniumSpecies const& Species)
  {
    if (Species == Q_1S0_8)
      return 5. / 12. * pow(alphas, 2) * pow(M_PI, 3) / pow(M, 5) * LDME.at(fit).at(Species);
    else if(Species == Q_3P0_8)
      return       5. * pow(alphas, 2) * pow(M_PI, 3) / pow(M, 5) * LDME.at(fit).at(Species);
    else if(Species == Q_3P2_8)
      return   4. / 3. * pow(alphas, 2) * pow(M_PI, 3) / pow(M, 7) * LDME.at(fit).at(Species);
    else if(Species == Q_3P2_1)
      return 64. / 15. * pow(alphas, 2) * pow(M_PI, 3) / pow(M, 7) * LDME.at(fit).at(Species);
    else
      throw std::runtime_error(error("gammaFg0Quarkonium", "Unknown quarkonium species"));

  }
}
