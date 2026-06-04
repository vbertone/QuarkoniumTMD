//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "quarkoniumtmd/matchingfunctionspdfquarkonium.h"

#include <apfel/constants.h>
#include <apfel/messages.h>

namespace apfel
{
  //_________________________________________________________________________________
  C1ggpdfQuarkonium::C1ggpdfQuarkonium(QuarkoniumSpecies const& Species):
    Expression()
  {
    // Reference: https://arxiv.org/pdf/hep-ph/9707223v2
    if (Species == Q_1S0_8){
      _BQc = CF * ( - 5 + Pi2 / 4 ) + CA * ( 2 + 7 * Pi2 / 24 ); // Eq. (129)
      std::cout << "BQc: " << _BQc << std::endl;
    }
    else if(Species == Q_3P0_8){
      _BQc = CF * ( - 7. / 3 + Pi2 / 4 ) + CA * ( 17. / 54 + 35.  / 27 * log(2)+ 7 * Pi2 / 24 ); // Eq. (130)
      std::cout << "BQc: " << _BQc << std::endl;
    }
    else if(Species == Q_3P2_8){
      _BQc = - 4 * CF + CA * ( 23. / 36 + 7 * log(2) / 9 + 5 * Pi2 / 12 ); // Eq. (131)
      std::cout << "BQc: " << _BQc << std::endl;
    }
    else if(Species == Q_3P2_1){
      _BQc = - 4 * CF + CA * ( 1. / 3 + 5 * log(2) / 3. + Pi2 / 6 ); // Eq. (127)
      std::cout << "BQc: " << _BQc << std::endl;
    }
    else
      throw std::runtime_error(error("C1ggpdfQuarkonium", "Unknown quarkonium species"));
  }
  double C1ggpdfQuarkonium::Local(double const&) const
  {
    // The factor of 2 in front of CA * zeta2 is suspicious
    // double Cgg = - 2 * CA * zeta2 + 2 * _BQc;
    // std::cout << "Cgg: " << Cgg / 4. << std::endl;
    return - 2 * CA * zeta2 + 2 * _BQc;
  }
}
