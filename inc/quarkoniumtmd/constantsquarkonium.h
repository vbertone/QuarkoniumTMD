//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <map>
#include <string>

namespace apfel
{
  enum QuarkoniumSpecies: int {Q_1S0_8, Q_3P0_8, Q_3P2_8, Q_3P2_1, Q_approx};

  /**
   * @name Quarkonium mass
   * @brief Mass of the quarkonia. To Do: choose PDG values or not
   */
  ///@{
  const double MJpsi    = 3.0;
  const double MUpsilon = 9.5;
  const double Mpsi2S   = 3.7;
  const double MChic2   = 3.556;
  const double MChib2   = 9.912;
  ///@}

  /**
   * @name LDMEs
   * @brief Quarkonia Long Distance Matrix Elements. 
   */
  ///@{
  const std::map<std::string, std::map<QuarkoniumSpecies, double>> LDME = {
    {"SV",   {{Q_1S0_8, 0.018},  {Q_3P0_8, 0.0405},    {Q_3P2_8, 0.02025},   {Q_3P2_1, 0.}}},
    {"BK11", {{Q_1S0_8, 0.0304}, {Q_3P0_8, -0.00908},  {Q_3P2_8, -0.0454},   {Q_3P2_1, 0.}}},
    {"C12",  {{Q_1S0_8, 0.089},  {Q_3P0_8, 0.0126},    {Q_3P2_8, 0.063},     {Q_3P2_1, 0.}}},
    {"SYY",  {{Q_1S0_8, 0.1423}, {Q_3P0_8, -0.039375}, {Q_3P2_8, -0.196875}, {Q_3P2_1, 0.018}, {Q_approx, 0.0197}}},
    };
  ///@}

    /**
   * @name Branching mass
   * @brief Mass of the quarkonia. To Do: choose PDG values or not
   */
  ///@{
  const std::map<std::string, std::map<std::string, double> > Br = {
    {"Jpsi",    {{"mumu",    5.93e-2}, {"ee", 5.94e-2}}},
    {"Upsilon", {{"mumu",    2.48e-2}, {"ee", 2.39e-2}}},
    {"Chic2",   {{"Jpsi",    0.19}}},
    {"Chib2",   {{"Upsilon", 0.18}}}
  };
  ///@}

/**
    * @name BQc coefficients
    * @brief: NLO Corrections for onium production.
    * Taken from Nucl. Phys. B 514 (1998) 245-309, see Eqs. (127), (129)-(131)
   */
   ///@{
   const std::map <std::string, double> BQc = {{"1S08", 3.16}, {"3P08", 3.76}, {"3P28", 2.8}, {"3P21", 4.07}};
   ///@}
   ///@}


}
