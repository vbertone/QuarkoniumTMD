//
// QuarkoniumTMD
//

#pragma once

#include "quarkoniumtmd/constantsquarkonium.h"

namespace apfel
{
  /**
   * @brief NLO virtual constant B of Petrelli et al., NPB 514 (1998)
   * 245 [hep-ph/9707223], Eqs. (127)-(131), normalised as
   *
   *   Gamma_V = Gamma_Born (alpha_s / pi) f_eps ( ... + B ).
   *
   * @param Species: quarkonium Fock state
   * @return B
   */
  double BQuarkonium(QuarkoniumSpecies const& Species);

  /**
   * @brief One-loop hard factor of g g -> Q[n]:
   *
   *   H = 1 + (alpha_s / pi) ( B - pi^2 C_A / 12 ).
   *
   * Used together with the universal gluon matching coefficient
   * C_gg^(1) = - C_A zeta2 delta(1-x) (coefficient of alpha_s / 4 pi),
   * it reproduces the O(alpha_s) constant of Sun, Yuan and Yuan,
   * PRD 88 (2013) 054008, Eq. (10).
   *
   * @param alphas: strong coupling alpha_s at the hard scale M
   * @param Species: quarkonium Fock state
   * @return H
   */
  double HardFactorQuarkonium(double const& alphas, QuarkoniumSpecies const& Species);
}
