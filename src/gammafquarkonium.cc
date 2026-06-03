//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "quarkoniumtmd/gammafquarkonium.h"

#include <apfel/constants.h>
#include <apfel/messages.h>

namespace apfel
{
  //_________________________________________________________________________
  double gammaFg0Quarkonium(int const& nf, QuarkoniumSpecies const& Species)
  {
    bool octet = false;
    if (Species == Q_1S0_8)
      octet = true;
    else if(Species == Q_3P0_8)
      octet = true;
    else if(Species == Q_3P2_8)
      octet = true;
    else if(Species == Q_approx)
      octet = true;
    else if(Species == Q_3P2_1)
      octet = false;
    else
      throw std::runtime_error(error("gammaFg0Quarkonium", "Unknown quarkonium species"));

    return - ( - 22 * CA / 3 + 8 * TR * nf / 3 - (octet ? 2 * CA : 0) );
  }
}
