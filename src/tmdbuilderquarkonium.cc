//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "quarkoniumtmd/matchingfunctionspdfquarkonium.h"
#include "quarkoniumtmd/gammafquarkonium.h"
#include "quarkoniumtmd/constantsquarkonium.h"

#include <apfel/tmdbuilder.h>
#include <apfel/timer.h>
#include <apfel/matchingfunctionspdf.h>
#include <apfel/evolutionbasisqcd.h>
#include <apfel/betaqcd.h>
#include <apfel/gammak.h>
#include <apfel/gammaf.h>
#include <apfel/kcs.h>
#include <apfel/hardfactors.h>
#include <apfel/tools.h>

namespace apfel
{
  //_____________________________________________________________________________
  std::map<int, TmdObjects> InitializeTmdObjectsQuarkonium(Grid                const& g,
                                                           std::vector<double> const& Thresholds,
                                                           QuarkoniumSpecies   const& Species,
                                                           double              const& IntEps)
  {
    // Initialise space-like splitting functions on the grid required
    // to compute the log terms of the matching functions.
    const std::map<int, DglapObjects> DglapObjpdf = InitializeDglapObjectsQCD(g, Thresholds, false, IntEps);

    report("Initializing TMD quarkonium objects for matching and evolution... ");
    Timer t;

    // Compute initial and final number of active flavours according
    // to the vector of thresholds (it assumes that the threshold
    // vector entries are ordered).
    int nfi = 0;
    int nff = Thresholds.size();
    for (auto const& v : Thresholds)
      if (v <= 0)
        nfi++;

    // ===============================================================
    // LO matching functions operators
    std::map<int, std::map<int, Operator>> C00;
    const Operator Id  {g, Identity{}, IntEps};
    const Operator Zero{g, Null{},     IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, Id});
        OM.insert({EvolutionBasisQCD::PNSM, Id});
        OM.insert({EvolutionBasisQCD::PNSV, Id});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * Id});
        OM.insert({EvolutionBasisQCD::PQG,  Zero});
        OM.insert({EvolutionBasisQCD::PGQ,  Zero});
        OM.insert({EvolutionBasisQCD::PGG,  Id});
        C00.insert({nf, OM});
      }

    // ===============================================================
    // NLO matching functions operators
    // PDFs
    std::map<int, std::map<int, Operator>> C10pdf;
    const Operator O1nspdf{g, C1nspdf{},                  IntEps};
    const Operator O1qgpdf{g, C1qgpdf{},                  IntEps};
    const Operator O1gqpdf{g, C1gqpdf{},                  IntEps};
    const Operator O1ggpdf{g, C1ggpdfQuarkonium{Species}, IntEps};
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O1nspdf});
        OM.insert({EvolutionBasisQCD::PNSM, O1nspdf});
        OM.insert({EvolutionBasisQCD::PNSV, O1nspdf});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O1nspdf});
        OM.insert({EvolutionBasisQCD::PQG,           nf * O1qgpdf});
        OM.insert({EvolutionBasisQCD::PGQ,  ( nf / 6. ) * O1gqpdf});
        OM.insert({EvolutionBasisQCD::PGG,                O1ggpdf});
        C10pdf.insert({nf, OM});
      }

    // Terms proportion to one power of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C11pdf;
    for (int nf = nfi; nf <= nff; nf++)
      {
        const Operator O11gmVq = gammaFq0() * Id;
        const Operator O11gmVg = gammaFg0Quarkonium(nf, Species) * Id;
        const auto P0 = DglapObjpdf.at(nf).SplittingFunctions.at(0);
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O11gmVq - 2 * P0.at(0)});
        OM.insert({EvolutionBasisQCD::PNSM, O11gmVq - 2 * P0.at(1)});
        OM.insert({EvolutionBasisQCD::PNSV, O11gmVq - 2 * P0.at(2)});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * ( O11gmVq - 2 * P0.at(3) )});
        OM.insert({EvolutionBasisQCD::PQG,                          - 2 * P0.at(4)});
        OM.insert({EvolutionBasisQCD::PGQ,                - ( nf / 3. ) * P0.at(5)});
        OM.insert({EvolutionBasisQCD::PGG,                  O11gmVg - 2 * P0.at(6)});
        C11pdf.insert({nf, OM});
      }

    // Terms proportion to two powers of log(mu0/mub)
    std::map<int, std::map<int, Operator>> C12;
    const Operator O12gmKq = - CF * gammaK0() / 2 * Id;
    const Operator O12gmKg = - CA * gammaK0() / 2 * Id;
    for (int nf = nfi; nf <= nff; nf++)
      {
        std::map<int, Operator> OM;
        OM.insert({EvolutionBasisQCD::PNSP, O12gmKq});
        OM.insert({EvolutionBasisQCD::PNSM, O12gmKq});
        OM.insert({EvolutionBasisQCD::PNSV, O12gmKq});
        OM.insert({EvolutionBasisQCD::PQQ,  ( nf / 6. ) * O12gmKq});
        OM.insert({EvolutionBasisQCD::PQG,                Zero});
        OM.insert({EvolutionBasisQCD::PGQ,                Zero});
        OM.insert({EvolutionBasisQCD::PGG,                O12gmKg});
        C12.insert({nf, OM});
      }

    // Define map containing the TmdObjects for each nf
    std::map<int, TmdObjects> TmdObj;

    // Construct a set of operators for each perturbative order for
    // the matching functions. Initialize also coefficients of: beta
    // function, gammaK, gammaF, and Collins-Soper anomalous
    // dimensions.
    for (int nf = nfi; nf <= nff; nf++)
      {
        TmdObjects obj;

        // Threshold
        obj.Threshold = Thresholds[nf-1];

        // Beta function
        obj.Beta.insert({0, beta0qcd(nf)});
        obj.Beta.insert({1, beta1qcd(nf)});

        // GammaF quark
        obj.GammaFq.insert({0, gammaFq0()});

        // GammaF gluon
        obj.GammaFg.insert({0, gammaFg0Quarkonium(nf, Species)});

        // gammaK (multiply by CF for quarks and by CA for gluons)
        obj.GammaK.insert({0, gammaK0()});
        obj.GammaK.insert({1, gammaK1(nf)});

        // Collins-Soper anomalous dimensions (multiply by CF for
        // quarks and by CA for gluons).
        obj.KCS.insert({0, {KCS00(), KCS01()}});

        // Matching functions
        const EvolutionBasisQCD evb{nf};

        // PDFs
        obj.MatchingFunctionsPDFs.insert({0, {{evb, C00.at(nf)}}});
        obj.MatchingFunctionsPDFs.insert({1, {{evb, C10pdf.at(nf)}, {evb, C11pdf.at(nf)}, {evb, C12.at(nf)}}});

        // Insert full object
        TmdObj.insert({nf, obj});
      }
    t.stop();

    return TmdObj;
  }
}
