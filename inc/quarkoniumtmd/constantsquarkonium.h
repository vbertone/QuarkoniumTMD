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
  /**
   * @brief Fock state of the heavy-quark pair produced in the hard
   * process. "Q_approx" is the approximated combination of the 1S0[8]
   * and 3P0[8] channels a la Sun-Yuan-Yuan.
   */
  enum QuarkoniumSpecies: int {Q_1S0_8, Q_3P0_8, Q_3P2_8, Q_3P2_1, Q_approx};

  /**
   * @brief Quarkonium family, i.e. the observed quarkonium state:
   * charmonium (J/psi, with the chi_c2 feed-down) or bottomonium
   * (Upsilon(1S), with the chi_b2 feed-down).
   */
  enum QuarkoniumFamily: int {Charmonium, Bottomonium};

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
   * @name Family lookups
   * @brief Mass of the observed quarkonium, mass of the chi_2 state
   * that feeds down into it, and names of the two states as used in
   * the "Br" map below, for each quarkonium family.
   */
  ///@{
  const std::map<QuarkoniumFamily, double> MOnium = {{Charmonium, MJpsi},   {Bottomonium, MUpsilon}};
  const std::map<QuarkoniumFamily, double> MChi2  = {{Charmonium, MChic2},  {Bottomonium, MChib2}};

  const std::map<QuarkoniumFamily, std::string> OniumName = {{Charmonium, "Jpsi"},  {Bottomonium, "Upsilon"}};
  const std::map<QuarkoniumFamily, std::string> Chi2Name  = {{Charmonium, "Chic2"}, {Bottomonium, "Chib2"}};
  ///@}

  /**
   * @name LDMEs
   * @brief Quarkonia Long Distance Matrix Elements, in GeV^3 for the
   * S-wave states and in GeV^5 for the P-wave ones. The map is
   * indexed by quarkonium family, then by the tag of the
   * determination, then by Fock state.
   */
  ///@{
  const std::map<QuarkoniumFamily, std::map<std::string, std::map<QuarkoniumSpecies, double>>> LDME =
  {
    {
      Charmonium,
      {
        {"SV",   {{Q_1S0_8, 0.018},  {Q_3P0_8, 0.0405},    {Q_3P2_8, 0.02025},   {Q_3P2_1, 0.}}},
        {"BK11", {{Q_1S0_8, 0.0304}, {Q_3P0_8, -0.00908},  {Q_3P2_8, -0.0454},   {Q_3P2_1, 0.}}},
        {"C12",  {{Q_1S0_8, 0.089},  {Q_3P0_8, 0.0126},    {Q_3P2_8, 0.063},     {Q_3P2_1, 0.}}},
        {"SYY",  {{Q_1S0_8, 0.1423}, {Q_3P0_8, -0.039375}, {Q_3P2_8, -0.196875}, {Q_3P2_1, 0.0896}, {Q_approx, 0.0197}}}
      }
    },
    {
      Bottomonium,
      {
        {"SYY",  {{Q_1S0_8, 0.119},  {Q_3P0_8, -0.27075},  {Q_3P2_8, -1.35375},  {Q_3P2_1, 1.973}, {Q_approx, 0.0321}}},
      }
    }
  };
  ///@}

  /**
   * @name Branching ratios
   * @brief Leptonic branching ratios of the observed quarkonia and
   * branching ratios of the chi_2 states into them.
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
