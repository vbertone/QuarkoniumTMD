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
  double sigmazero(double const& alphas, double const& M, std::string const& fit, QuarkoniumSpecies const& Species, QuarkoniumFamily const& Family)
  {
    // Retrieve the LDME of the requested family, determination and
    // Fock state, reporting explicitly what is missing if the
    // combination is not available.
    if (LDME.count(Family) == 0)
      throw std::runtime_error(error("sigmazero", "Unknown quarkonium family"));

    const std::map<std::string, std::map<QuarkoniumSpecies, double>>& LDMEFamily = LDME.at(Family);
    if (LDMEFamily.count(fit) == 0)
      throw std::runtime_error(error("sigmazero", "No '" + fit + "' LDME determination available for the "
                                     + std::string(Family == Charmonium ? "charmonium" : "bottomonium")));

    const std::map<QuarkoniumSpecies, double>& LDMEFit = LDMEFamily.at(fit);
    if (LDMEFit.count(Species) == 0)
      throw std::runtime_error(error("sigmazero", "The '" + fit + "' determination does not provide the LDME of the requested Fock state"));

    const double ldme = LDMEFit.at(Species);

    if (Species == Q_1S0_8)
      return 5. / 12. * pow(alphas, 2) * pow(M_PI, 3) / pow(M, 5) * ldme;
    else if(Species == Q_3P0_8)
      return       5. * pow(alphas, 2) * pow(M_PI, 3) / pow(M, 7) * ldme;
    else if(Species == Q_3P2_8)
      return   4. / 3. * pow(alphas, 2) * pow(M_PI, 3) / pow(M, 7) * ldme;
    else if(Species == Q_3P2_1)
      return 64. / 15. * pow(alphas, 2) * pow(M_PI, 3) / pow(M, 7) * ldme;
    else if(Species == Q_approx)
      // @note No B_Qc factor here: the NLO coefficients of
      // Nucl. Phys. B 514 (1998) 245 enter the matching function
      // C1ggpdfQuarkonium.
      return 5. / 12. * pow(alphas, 2) * pow(M_PI, 3) / pow(M, 5) * ldme;
    else
      throw std::runtime_error(error("sigmazero", "Unknown quarkonium species"));

  }
}
