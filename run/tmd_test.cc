//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#include <quarkoniumtmd/tmdbuilderquarkonium.h>

// Main program
int main()
{
  // Vectors of masses and thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};
  //const std::vector<double> Thresholds = {0, 0, 0, 0, 0};

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, apfel::FixedOrderAccuracy::NLO};
  const apfel::TabulateObject<double> tabas{a, 500, 0.9, 1001, 3};
  const auto Alphas = [=] (double const& mu) -> double{ return tabas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{100, 1e-1, 3}, apfel::SubGrid{100, 6e-1, 3}, apfel::SubGrid{80, 8.5e-1, 5}}};

  // Construct the DGLAP objects
  const auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Thresholds), apfel::LHToyPDFs, sqrt(2), apfel::FixedOrderAccuracy::NLO, Alphas);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabPDFs{*EvolvedPDFs, 200, 1, 20000, 3};
  const auto CollPDFs = [=] (double const& mu) -> apfel::Set<apfel::Distribution> { return TabPDFs.Evaluate(mu); };

  // Get timer
  apfel::Timer t;

  // Initialize TMD objects
  const auto TmdObj = apfel::InitializeTmdObjectsQuarkonium(g, Thresholds, apfel::QuarkoniumSpecies::Q_1S0_8);

  // Kinematics
  const double Vs = 13000;
  const double y  = 0;
  const double Q  = 3;

  // Compute 'x1' and 'x2'
  const double x1 = Q * exp(y) / Vs;
  const double x2 = Q / exp(y) / Vs;

  // Build evolved TMD PDFs
  const auto EvTMDPDFs = BuildTmdPDFs(TmdObj, CollPDFs, Alphas, apfel::LogAccuracy::NLLp);

  // Get Drell-Yan hard-factor function
  const double hcs = 1;//HardFactor("DY", TmdObj, Alphas, apfel::LogAccuracy::NNNLL)(Q);

  // Get Evolved TMD PDFs at bmax. This is subtracted off from the
  // b-dependent luminosity to improve the convergence of the Hankel
  // transform, particularly at small qT.
  const std::map<int,apfel::Distribution> xF = QCDEvToPhys(EvTMDPDFs(2 * exp( - apfel::emc), Q, Q * Q).GetObjects());
  const double lumimax = xF.at(0).Evaluate(x1) * xF.at(0).Evaluate(x2);

  // Construct the TMD luminosity in b space to be trasformed in qT
  // space.
  const std::function<double(double const&)> TMDLumib = [=] (double const& b) -> double
  {
    // Compute b*
    const double bstar = std::max(b / sqrt( 1 + pow(b / 2 / exp( - apfel::emc), 2) ), 1e-4);

    // Get Evolved TMD PDFs and rotate them into the physical
    // basis
    const std::map<int,apfel::Distribution> xF = QCDEvToPhys(EvTMDPDFs(bstar, Q, Q * Q).GetObjects());

    // Combine TMDs through the EW charges
    double lumi = xF.at(0).Evaluate(x1) * xF.at(0).Evaluate(x2);

    // Combine all pieces and return
    return b * ( lumi - lumimax );
  };

  // Double exponential quadrature
  //apfel::DoubleExponentialQuadrature DEObj{};
  apfel::OgataQuadrature DEObj{};

  // Phase-space reduction factor
  apfel::TwoBodyPhaseSpace ps{20, -1, 2.4};

  // Compute predictions
  const int nqT = 100;
  const double qTmin = 0.01;
  const double qTmax = 3;
  const double qTstp = ( qTmax - qTmin ) / ( nqT - 1 );
  std::cout << std::scientific;
  std::cout << "\n     Q         "
            << "       y        "
            << "      qT        "
            << " cross section  "
            << " PS red. fact.  "
            << "PV PS red. fact."
            << std::endl;
  for (double qT = qTmin; qT <= qTmax; qT += qTstp)
    std::cout << Q << "\t" << y << "\t" << qT << "\t"
              << apfel::ConvFact * qT * 8 * M_PI * pow(apfel::alphaem, 2) * hcs / pow(Q, 3) / 9 * DEObj.transform(TMDLumib, qT) << "\t"
              << ps.PhaseSpaceReduction(Q, y, qT) << "\t"
              << ps.ParityViolatingPhaseSpaceReduction(Q, y, qT) << "\t"
              << std::endl;
  t.stop();

  return 0;
}
