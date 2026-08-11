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
  /**
   * @brief Leading-order partonic cross section for the production of
   * a heavy-quark pair in a given Fock state, multiplied by the
   * corresponding long-distance matrix element.
   * @param alphas: the strong coupling evaluated at the relevant scale
   * @param M: the mass of the state produced in the hard process
   * @param fit: the tag of the LDME determination (e.g. "SYY")
   * @param Species: the Fock state
   * @param Family: the quarkonium family (default: Charmonium)
   * @return sigma_0 x LDME
   */
  double sigmazero(double            const& alphas,
                   double            const& M,
                   std::string       const& fit,
                   QuarkoniumSpecies const& Species,
                   QuarkoniumFamily  const& Family = Charmonium);
}
